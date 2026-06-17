#ifndef POOL_H
#define POOL_H

#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "hwloc_topo.h"

// Forward declarations resolved by inclusion order in main.c:
//   process_target() is declared before pool.h (via process.h)

typedef struct l3_dispatcher l3_dispatcher_t;

typedef struct {
    int pool_idx;
    int n_workers;
    parameters* params;
    topo_config_t* topo;
    kt_forpool_t* pool;
    buffer_t** buffers;  // points into params->buffers[buffer_offset..]
    int buffer_offset;
    bool* worker_bound;  // [n_workers] lazy binding flags
    long chunk_offset;   // current chunk start (set by pool manager before each kt_forpool call)
    l3_dispatcher_t* disp;
} pool_context_t;

struct l3_dispatcher {
    int num_pools;
    pool_context_t pools[MAX_L3_POOLS];
    volatile long next_chunk;
    long chunk_size;
    long total_files;
    pthread_t managers[MAX_L3_POOLS];
    topo_config_t* topo;
};

// Pool worker callback: lazy-bind to PU, allocate buffer, process file.
static void pool_process_targets(void* data, long i, int thread_id) {
    pool_context_t* ctx = (pool_context_t*)data;

    if (!ctx->worker_bound[thread_id]) {
        bind_thread_pu(ctx->topo, ctx->pool_idx, thread_id);
        if (!ctx->buffers[thread_id]->allocated) {
            alloc_buffer(ctx->buffers[thread_id], ctx->params);
        }
        ctx->worker_bound[thread_id] = true;
    }

    long global_i = ctx->chunk_offset + i;
    process_target(kv_A(ctx->params->targets, global_i).path, ctx->buffers[thread_id], ctx->params, true, NULL);

    int total = (int)ctx->disp->total_files;
    int permile = total / 1000;
    if (permile == 0) permile = 1;
    int current = __sync_add_and_fetch(&ctx->params->iter_count, 1);
    if (current % permile == 0 || current >= total) {
        float progress = (float)current * 100.0f / (float)total;
        double elapsed = elapsed_seconds_since(&ctx->params->batch_start_time);
        char time_buf[64];
        if (elapsed > 0.001 && current > 0) {
            double remaining_sec = (double)(total - current) * elapsed / (double)current;
            format_time_remaining(remaining_sec, time_buf, sizeof(time_buf));
        } else {
            snprintf(time_buf, sizeof(time_buf), "calculating...");
        }
        printf("\rComputation in progress: %.1f%% complete | %s%s", progress, time_buf, current >= total ? "\n" : "");
        fflush(stdout);
    }
}

static void* pool_manager_thread(void* arg) {
    pool_context_t* ctx = (pool_context_t*)arg;
    l3_dispatcher_t* d = ctx->disp;

    for (;;) {
        long chunk_idx = __sync_fetch_and_add(&d->next_chunk, 1L);
        long start = chunk_idx * d->chunk_size;
        if (start >= d->total_files) break;

        long count = d->chunk_size;
        if (start + count > d->total_files) count = d->total_files - start;

        ctx->chunk_offset = start;
        kt_forpool(ctx->pool, pool_process_targets, ctx, count);
    }
    return NULL;
}

static l3_dispatcher_t* l3_dispatcher_create(parameters* params, topo_config_t* topo) {
    l3_dispatcher_t* d = (l3_dispatcher_t*)calloc(1, sizeof(l3_dispatcher_t));
    if (!d) return NULL;
    d->total_files = (long)kv_size(params->targets);
    d->chunk_size = topo->chunk_size;
    d->next_chunk = 0;
    d->topo = topo;

    // Allocate a flat buffer array covering all pool workers.
    // params->buffers is populated so the existing fprint_buffer loop works unchanged.
    params->nbuffers = topo->total_workers;
    params->buffers = (buffer_t**)calloc((size_t)topo->total_workers, sizeof(buffer_t*));
    for (int i = 0; i < topo->total_workers; i++) {
        params->buffers[i] = (buffer_t*)calloc(1, sizeof(buffer_t));
    }

    int running_offset = 0;
    d->num_pools = 0;
    for (int p = 0; p < topo->num_pools; p++) {
        if (topo->threads_per_pool[p] == 0) continue;
        pool_context_t* ctx = &d->pools[d->num_pools];
        ctx->pool_idx = p;
        ctx->n_workers = topo->threads_per_pool[p];
        ctx->params = params;
        ctx->topo = topo;
        ctx->disp = d;
        ctx->buffer_offset = running_offset;
        ctx->buffers = &params->buffers[running_offset];
        ctx->worker_bound = (bool*)calloc((size_t)ctx->n_workers, sizeof(bool));
        ctx->chunk_offset = 0;

        if (ctx->n_workers > 1) {
            ctx->pool = (kt_forpool_t*)kt_forpool_init(ctx->n_workers, params->idle);
        } else {
            ctx->pool = NULL;
        }
        running_offset += ctx->n_workers;
        d->num_pools++;
    }
    return d;
}

static void l3_dispatcher_run(l3_dispatcher_t* d) {
    // Launch one pool manager pthread per L3 domain
    for (int p = 0; p < d->num_pools; p++) {
        if (pthread_create(&d->managers[p], NULL, pool_manager_thread, &d->pools[p]) != 0) {
            perror("pthread_create (pool manager)");
            exit(EXIT_FAILURE);
        }
    }
    for (int p = 0; p < d->num_pools; p++) {
        pthread_join(d->managers[p], NULL);
    }
}

static void l3_dispatcher_destroy(l3_dispatcher_t* d) {
    if (!d) return;
    for (int p = 0; p < d->num_pools; p++) {
        pool_context_t* ctx = &d->pools[p];
        if (ctx->pool) kt_forpool_destroy(ctx->pool);
        free(ctx->worker_bound);
        ctx->worker_bound = NULL;
    }
    free(d);
}

#endif  // POOL_H
