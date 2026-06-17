#ifndef HWLOC_TOPO_H
#define HWLOC_TOPO_H

#include <hwloc.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef MAX_L3_POOLS
#    define MAX_L3_POOLS 16
#endif
#define MAX_CORES_PER_POOL 128
#define MAX_SMT_PER_CORE 8

typedef struct {
    int pool_idx;
    int worker_idx;
    int core_idx;
    int pu_idx;
    int os_pu_index;
} pu_binding_t;

typedef struct {
    int domain_idx;
    int num_cores;
    int num_pus;
    int smt_per_core;
    pu_binding_t phase1[MAX_CORES_PER_POOL];
    int phase1_count;
    pu_binding_t phase2[MAX_CORES_PER_POOL * (MAX_SMT_PER_CORE - 1)];
    int phase2_count;
} l3_domain_t;

typedef struct topo_config {
    int num_pools;
    int threads_per_pool[MAX_L3_POOLS];
    int total_workers;
    bool smt_available;
    int smt_per_core;
    l3_domain_t domains[MAX_L3_POOLS];
    pu_binding_t* all_bindings[MAX_L3_POOLS];
    int all_bindings_count[MAX_L3_POOLS];
    int total_cores;
    int total_physical_pus;
    int chunk_size;
    hwloc_topology_t topo;
} topo_config_t;

static void topo_build_domains(hwloc_topology_t topology, topo_config_t* cfg) {
    int l3_depth = hwloc_get_type_depth(topology, HWLOC_OBJ_L3CACHE);
    int num_l3 = (l3_depth != HWLOC_TYPE_DEPTH_UNKNOWN) ? hwloc_get_nbobjs_by_depth(topology, l3_depth) : 0;
    bool synthetic_single = (num_l3 <= 0);
    if (synthetic_single) num_l3 = 1;

    cfg->num_pools = num_l3 > MAX_L3_POOLS ? MAX_L3_POOLS : num_l3;
    cfg->total_cores = 0;
    cfg->total_physical_pus = 0;
    cfg->smt_per_core = 0;

    for (int d = 0; d < cfg->num_pools; d++) {
        l3_domain_t* dom = &cfg->domains[d];
        dom->domain_idx = d;

        hwloc_cpuset_t search_set;
        if (!synthetic_single) {
            hwloc_obj_t l3_obj = hwloc_get_obj_by_depth(topology, l3_depth, d);
            search_set = l3_obj->cpuset;
        } else {
            hwloc_obj_t machine = hwloc_get_obj_by_type(topology, HWLOC_OBJ_MACHINE, 0);
            search_set = machine->cpuset;
        }

        dom->num_cores = hwloc_get_nbobjs_inside_cpuset_by_type(topology, search_set, HWLOC_OBJ_CORE);
        dom->num_pus = hwloc_get_nbobjs_inside_cpuset_by_type(topology, search_set, HWLOC_OBJ_PU);
        dom->smt_per_core = (dom->num_cores > 0) ? (dom->num_pus / dom->num_cores) : 1;
        if (dom->smt_per_core < 1) dom->smt_per_core = 1;
        if (dom->smt_per_core > cfg->smt_per_core) cfg->smt_per_core = dom->smt_per_core;
        cfg->total_cores += dom->num_cores;
        cfg->total_physical_pus += dom->num_pus;

        // Phase 1: one thread per core (first PU of each core)
        dom->phase1_count = 0;
        hwloc_obj_t core = NULL;
        int core_idx = 0;
        while ((core = hwloc_get_next_obj_inside_cpuset_by_type(topology, search_set, HWLOC_OBJ_CORE, core)) != NULL) {
            hwloc_obj_t pu = hwloc_get_next_obj_inside_cpuset_by_type(topology, core->cpuset, HWLOC_OBJ_PU, NULL);
            if (pu && dom->phase1_count < MAX_CORES_PER_POOL) {
                pu_binding_t* b = &dom->phase1[dom->phase1_count];
                b->pool_idx = d;
                b->worker_idx = dom->phase1_count;
                b->core_idx = core_idx;
                b->pu_idx = 0;
                b->os_pu_index = (int)pu->os_index;
                dom->phase1_count++;
            }
            core_idx++;
        }

        // Phase 2: additional SMT siblings (PU1, PU2, ... per core)
        dom->phase2_count = 0;
        core = NULL;
        core_idx = 0;
        while ((core = hwloc_get_next_obj_inside_cpuset_by_type(topology, search_set, HWLOC_OBJ_CORE, core)) != NULL) {
            hwloc_obj_t pu = NULL;
            int pu_rank = 0;
            while ((pu = hwloc_get_next_obj_inside_cpuset_by_type(topology, core->cpuset, HWLOC_OBJ_PU, pu)) != NULL) {
                if (pu_rank > 0 && dom->phase2_count < MAX_CORES_PER_POOL * (MAX_SMT_PER_CORE - 1)) {
                    pu_binding_t* b = &dom->phase2[dom->phase2_count];
                    b->pool_idx = d;
                    b->core_idx = core_idx;
                    b->pu_idx = pu_rank;
                    b->os_pu_index = (int)pu->os_index;
                    dom->phase2_count++;
                }
                pu_rank++;
            }
            core_idx++;
        }
    }

    cfg->smt_available = (cfg->smt_per_core > 1);
}

static void topo_assign_workers(topo_config_t* cfg, int max_total_threads) {
    int budget = max_total_threads;
    if (budget < 1) budget = cfg->total_physical_pus;

    // Determine cores_per_pool (max physical cores available in any single domain)
    int cores_per_pool = 0;
    for (int d = 0; d < cfg->num_pools; d++) {
        if (cfg->domains[d].phase1_count > cores_per_pool) cores_per_pool = cfg->domains[d].phase1_count;
    }
    if (cores_per_pool < 1) cores_per_pool = 1;

    // Determine how many pools to activate.
    // If budget >= total_physical_pus (or default), use all pools.
    // Otherwise create just enough pools to hold the threads.
    int required_pools;
    if (budget >= cfg->total_physical_pus) {
        required_pools = cfg->num_pools;
    } else {
        required_pools = (budget + cores_per_pool - 1) / cores_per_pool;
        if (required_pools > cfg->num_pools) required_pools = cfg->num_pools;
    }
    if (required_pools < 1) required_pools = 1;

    // Distribute threads evenly across required pools.
    // First 'remainder' pools get base+1 threads, the rest get base.
    int base = budget / required_pools;
    int remainder = budget % required_pools;

    for (int d = 0; d < cfg->num_pools; d++) {
        if (d < required_pools) {
            cfg->threads_per_pool[d] = base + (d < remainder ? 1 : 0);
        } else {
            cfg->threads_per_pool[d] = 0;
        }
    }

    // Build per-pool flat binding arrays: phase1 entries first, then phase2 (SMT siblings)
    int min_threads = INT32_MAX;
    cfg->total_workers = 0;
    for (int d = 0; d < cfg->num_pools; d++) {
        l3_domain_t* dom = &cfg->domains[d];
        int n = cfg->threads_per_pool[d];
        cfg->all_bindings[d] = (pu_binding_t*)calloc(n > 0 ? (size_t)n : 1U, sizeof(pu_binding_t));
        cfg->all_bindings_count[d] = n;

        int filled = 0;
        for (int i = 0; i < dom->phase1_count && filled < n; i++) {
            cfg->all_bindings[d][filled] = dom->phase1[i];
            cfg->all_bindings[d][filled].worker_idx = filled;
            filled++;
        }
        for (int i = 0; i < dom->phase2_count && filled < n; i++) {
            cfg->all_bindings[d][filled] = dom->phase2[i];
            cfg->all_bindings[d][filled].worker_idx = filled;
            filled++;
        }
        cfg->total_workers += n;
        if (n > 0 && n < min_threads) min_threads = n;
    }

    // chunk_size: smallest multiple of min_threads_per_pool that is strictly greater than the number of active pools
    int active_pools = required_pools;
    if (min_threads > 0 && min_threads < INT32_MAX) {
        int k = (active_pools / min_threads) + 1;
        cfg->chunk_size = k * min_threads;
    } else {
        cfg->chunk_size = 1;
    }
}

static topo_config_t* probe_topology(int max_total_threads) {
    hwloc_topology_t topology;
    if (hwloc_topology_init(&topology) != 0) {
        fprintf(stderr, "hwloc_topology_init failed\n");
        return NULL;
    }
    if (hwloc_topology_load(topology) != 0) {
        fprintf(stderr, "hwloc_topology_load failed\n");
        hwloc_topology_destroy(topology);
        return NULL;
    }

    topo_config_t* cfg = (topo_config_t*)calloc(1, sizeof(topo_config_t));
    if (!cfg) {
        hwloc_topology_destroy(topology);
        return NULL;
    }
    cfg->topo = topology;
    topo_build_domains(topology, cfg);
    topo_assign_workers(cfg, max_total_threads);
    return cfg;
}

static void print_topology_summary(const topo_config_t* cfg) {
    printf("Hardware topology:\n");
    printf("  L3 cache domains (threadpools): %d\n", cfg->num_pools);
    printf("  Total physical cores:           %d\n", cfg->total_cores);
    printf("  Total PUs (incl. SMT):          %d\n", cfg->total_physical_pus);
    printf("  SMT per core:                   %d (%s)\n", cfg->smt_per_core, cfg->smt_available ? "enabled" : "disabled");
    printf("  Total workers allocated:        %d\n", cfg->total_workers);
    printf("  Dispatch chunk size:            %d\n", cfg->chunk_size);

    for (int d = 0; d < cfg->num_pools; d++) {
        const l3_domain_t* dom = &cfg->domains[d];
        printf("  Pool %d (L3#%d): %d cores, %d PUs, %d workers\n", d, dom->domain_idx, dom->num_cores, dom->num_pus, cfg->threads_per_pool[d]);
        if (cfg->threads_per_pool[d] > 0 && cfg->all_bindings[d]) {
            int show = cfg->threads_per_pool[d] < 8 ? cfg->threads_per_pool[d] : 8;
            printf("    PUs: ");
            for (int w = 0; w < show; w++) {
                printf("%d%s", cfg->all_bindings[d][w].os_pu_index, w < show - 1 ? "," : "");
            }
            if (cfg->threads_per_pool[d] > 8) printf(",...");
            printf("\n");
        }
    }
}

static int save_topology_config(const topo_config_t* cfg) {
    const char* home = getenv("HOME");
    if (!home) home = "/tmp";

    char dirpath[512];
    snprintf(dirpath, sizeof(dirpath), "%s/.ihsnpeaks", home);
    (void)mkdir(dirpath, 0755);

    char xmlpath[512];
    snprintf(xmlpath, sizeof(xmlpath), "%s/hwloc_config", dirpath);

    if (hwloc_topology_export_xml(cfg->topo, xmlpath, 0) != 0) {
        fprintf(stderr, "Failed to export hwloc topology to %s\n", xmlpath);
        return -1;
    }
    printf("Hardware topology saved to %s\n", xmlpath);
    print_topology_summary(cfg);
    return 0;
}

static topo_config_t* load_or_probe_topology(bool do_generate, int max_total_threads) {
    topo_config_t* cfg = probe_topology(max_total_threads);
    if (!cfg) return NULL;

    if (do_generate) {
        save_topology_config(cfg);
        return cfg;
    }

    // If no cached config exists yet, auto-save one
    const char* home = getenv("HOME");
    if (home) {
        char xmlpath[512];
        snprintf(xmlpath, sizeof(xmlpath), "%s/.ihsnpeaks/hwloc_config", home);
        struct stat st;
        if (stat(xmlpath, &st) != 0) {
            save_topology_config(cfg);
        }
    }
    return cfg;
}

static void bind_thread_pu(const topo_config_t* cfg, int pool_idx, int worker_idx) {
    if (pool_idx < 0 || pool_idx >= cfg->num_pools) return;
    if (worker_idx < 0 || worker_idx >= cfg->all_bindings_count[pool_idx]) return;

    const pu_binding_t* b = &cfg->all_bindings[pool_idx][worker_idx];
    hwloc_cpuset_t cpuset = hwloc_bitmap_alloc();
    hwloc_bitmap_only(cpuset, (unsigned)b->os_pu_index);

    if (hwloc_set_cpubind(cfg->topo, cpuset, HWLOC_CPUBIND_THREAD) != 0) {
        fprintf(stderr, "Warning: cpubind failed for pool %d worker %d (PU %d)\n", pool_idx, worker_idx, b->os_pu_index);
    }
    hwloc_bitmap_free(cpuset);
}

static void free_topology(topo_config_t* cfg) {
    if (!cfg) return;
    for (int d = 0; d < cfg->num_pools; d++) {
        free(cfg->all_bindings[d]);
        cfg->all_bindings[d] = NULL;
    }
    if (cfg->topo) hwloc_topology_destroy(cfg->topo);
    free(cfg);
}

#endif  // HWLOC_TOPO_H
