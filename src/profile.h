#ifndef PROFILE_H
#define PROFILE_H

#ifdef IHSNPEAKS_PROFILE

#    include <stdint.h>
#    include <stdio.h>
#    include <stdlib.h>
#    include <string.h>
#    include <time.h>

typedef enum {
    PROFILE_PHASE_READ = 0,
    PROFILE_PHASE_PREPROCESS,
    PROFILE_PHASE_PRECOMPUTE,
    PROFILE_PHASE_FILL_INPUT,
    PROFILE_PHASE_LADDER_SETUP,
    PROFILE_PHASE_RESET_LADDER,
    PROFILE_PHASE_NUFFT_EXECUTE,
    PROFILE_PHASE_POWER_ACCUMULATE,
    PROFILE_PHASE_LADDER_ADVANCE,
    PROFILE_PHASE_PEAK_SCAN,
    PROFILE_PHASE_SORT,
    PROFILE_PHASE_OUTPUT,
    PROFILE_PHASE_COUNT
} profile_phase_t;

static unsigned long long profile_ns[PROFILE_PHASE_COUNT];
static unsigned long long profile_calls[PROFILE_PHASE_COUNT];
static unsigned long long profile_targets;

#    if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#        define IHSNPEAKS_PROFILE_THREAD_LOCAL _Thread_local
#    else
#        define IHSNPEAKS_PROFILE_THREAD_LOCAL __thread
#    endif

static IHSNPEAKS_PROFILE_THREAD_LOCAL unsigned long long profile_local_ns[PROFILE_PHASE_COUNT];
static IHSNPEAKS_PROFILE_THREAD_LOCAL unsigned long long profile_local_calls[PROFILE_PHASE_COUNT];
static IHSNPEAKS_PROFILE_THREAD_LOCAL unsigned long long profile_local_targets;

static inline uint64_t profile_now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ((uint64_t)ts.tv_sec * UINT64_C(1000000000)) + (uint64_t)ts.tv_nsec;
}

static inline void profile_add(profile_phase_t phase, uint64_t start_ns) {
    uint64_t elapsed = profile_now_ns() - start_ns;
    profile_local_ns[phase] += (unsigned long long)elapsed;
    profile_local_calls[phase] += 1ULL;
}

static inline void profile_count_target(void) { ++profile_local_targets; }

static inline void profile_flush_thread(void) {
    if (profile_local_targets) {
        __sync_fetch_and_add(&profile_targets, profile_local_targets);
        profile_local_targets = 0;
    }
    for (int i = 0; i < PROFILE_PHASE_COUNT; ++i) {
        if (profile_local_ns[i]) {
            __sync_fetch_and_add(&profile_ns[i], profile_local_ns[i]);
            profile_local_ns[i] = 0;
        }
        if (profile_local_calls[i]) {
            __sync_fetch_and_add(&profile_calls[i], profile_local_calls[i]);
            profile_local_calls[i] = 0;
        }
    }
}

#    define PROFILE_START(name) uint64_t profile_start_##name = profile_now_ns()
#    define PROFILE_END(name, phase) profile_add((phase), profile_start_##name)
#    define PROFILE_COUNT_TARGET() profile_count_target()
#    define PROFILE_FLUSH_THREAD() profile_flush_thread()

static inline double profile_seconds(unsigned long long ns) { return (double)ns / 1.0e9; }

static inline double profile_mib(size_t bytes) { return (double)bytes / (1024.0 * 1024.0); }

static inline size_t profile_parse_cache_size(const char *text) {
    char *end = NULL;
    unsigned long value = strtoul(text, &end, 10);
    while (end && (*end == ' ' || *end == '\t')) ++end;
    if (!end || value == 0) return 0;
    if (*end == 'K' || *end == 'k') return (size_t)value * 1024U;
    if (*end == 'M' || *end == 'm') return (size_t)value * 1024U * 1024U;
    return (size_t)value;
}

static inline size_t profile_read_cache_size(int wanted_level, const char *wanted_type) {
    for (int idx = 0; idx < 10; ++idx) {
        char path[128];
        char buf[64];
        snprintf(path, sizeof(path), "/sys/devices/system/cpu/cpu0/cache/index%d/level", idx);
        FILE *f = fopen(path, "r");
        if (!f) continue;
        int level = 0;
        if (fscanf(f, "%d", &level) != 1) level = 0;
        fclose(f);
        if (level != wanted_level) continue;

        snprintf(path, sizeof(path), "/sys/devices/system/cpu/cpu0/cache/index%d/type", idx);
        f = fopen(path, "r");
        if (!f) continue;
        if (!fgets(buf, sizeof(buf), f)) buf[0] = '\0';
        fclose(f);
        buf[strcspn(buf, "\r\n")] = '\0';
        if (strcmp(buf, wanted_type) != 0) continue;

        snprintf(path, sizeof(path), "/sys/devices/system/cpu/cpu0/cache/index%d/size", idx);
        f = fopen(path, "r");
        if (!f) continue;
        if (!fgets(buf, sizeof(buf), f)) buf[0] = '\0';
        fclose(f);
        return profile_parse_cache_size(buf);
    }
    return 0;
}

static inline uint32_t profile_padded_len(uint32_t len) {
    uint32_t lanes = (uint32_t)VEC_LEN;
    return (uint32_t)(((size_t)len + (size_t)lanes - 1U) & ~((size_t)lanes - 1U));
}

static inline void profile_report_cache_line(const char *label, size_t bytes, double l1_mib, double l2_mib, double l3_mib) {
    fprintf(stderr, "  %-24s %8.2f MiB", label, profile_mib(bytes));
    if (l1_mib > 0.0 && profile_mib(bytes) <= l1_mib) {
        fprintf(stderr, "  <= L1d/thread");
    } else if (l2_mib > 0.0 && profile_mib(bytes) <= l2_mib) {
        fprintf(stderr, "  <= L2/thread");
    } else if (l3_mib > 0.0 && profile_mib(bytes) <= l3_mib) {
        fprintf(stderr, "  <= shared L3");
    } else {
        fprintf(stderr, "  > shared L3");
    }
    fputc('\n', stderr);
}

static inline void profile_report(const parameters *params) {
    static const char *phase_names[PROFILE_PHASE_COUNT] = {
        "read_dat",      "preprocess",       "nufft_precompute", "fill_input",         "ladder_setup", "reset_ladder",
        "nufft_execute", "power_accumulate", "ladder_advance",   "peak_scan/spectrum", "sort/refine",  "output",
    };

    fflush(stdout);
    fprintf(stderr, "\n[ihsnpeaks profile]\n");
    fprintf(stderr, "targets: %llu\n", profile_targets);
    fprintf(stderr, "threads: %d\n", params ? (params->nbuffers > 0 ? params->nbuffers : 1) : 0);
    fprintf(stderr, "maxLen: %u  maxFreqCount: %u  gridLen: %u  outputLen: %u  ladderLevels: %u\n", params ? params->maxLen : 0U,
            params ? params->maxFreqCount : 0U, params ? params->gridLen : 0U, params ? params->outputLen : 0U, params ? params->ladderLevels : 0U);

    fprintf(stderr, "\nphase totals, summed across threads:\n");
    for (int i = 0; i < PROFILE_PHASE_COUNT; ++i) {
        fprintf(stderr, "  %-20s %10.6f s  calls=%llu\n", phase_names[i], profile_seconds(profile_ns[i]), profile_calls[i]);
    }

    if (!params) return;

    uint32_t padded_len = profile_padded_len(params->maxLen);
    size_t time_series = (size_t)params->maxLen * (sizeof(double) + 3U * sizeof(float));
    size_t power_grid = ((size_t)params->maxFreqCount + 2U) * sizeof(float);
    size_t block_io = (size_t)params->outputLen * 2U * sizeof(float);
    size_t ladder = (size_t)NUFFT_LADDER_LEVEL_CAP * (size_t)padded_len * 4U * sizeof(float);
    size_t fft = params->nufftExternalSizes.fft_len * 2U * sizeof(float);
    size_t cobra = params->nufftExternalSizes.cobra_len * 2U * sizeof(float);
    size_t input = (size_t)padded_len * 2U * sizeof(float);
    size_t read_buf = params->maxSize;
    size_t nufft_workspace = ((size_t)params->maxLen * (size_t)params->nterms * sizeof(int)) +
                             ((size_t)params->maxLen * (size_t)params->nterms * 2U * sizeof(float)) +
                             ((size_t)params->maxLen * (size_t)params->nterms * (size_t)NUFFT1_PSWF_WIDTH * sizeof(float));
    size_t per_thread_total = time_series + power_grid + block_io + ladder + fft + cobra + input + read_buf + nufft_workspace;
    size_t shared_twiddles = params->nufftExternalSizes.twiddle_len * 2U * sizeof(float);

    size_t l1_bytes = profile_read_cache_size(1, "Data");
    size_t l2_bytes = profile_read_cache_size(2, "Unified");
    size_t l3_bytes = profile_read_cache_size(3, "Unified");
    double l1_mib = l1_bytes ? profile_mib(l1_bytes) : 32.0 / 1024.0;
    double l2_mib = l2_bytes ? profile_mib(l2_bytes) : 1.0;
    double l3_mib = l3_bytes ? profile_mib(l3_bytes) : 64.0;

    fprintf(stderr, "\nworking-set estimate per worker:\n");
    profile_report_cache_line("read buffer", read_buf, l1_mib, l2_mib, l3_mib);
    profile_report_cache_line("x/y/dy/wy", time_series, l1_mib, l2_mib, l3_mib);
    profile_report_cache_line("power grid", power_grid, l1_mib, l2_mib, l3_mib);
    profile_report_cache_line("block real/imag", block_io, l1_mib, l2_mib, l3_mib);
    profile_report_cache_line("input real/imag", input, l1_mib, l2_mib, l3_mib);
    profile_report_cache_line("twiddle ladder arrays", ladder, l1_mib, l2_mib, l3_mib);
    profile_report_cache_line("fft/cobra buffers", fft + cobra, l1_mib, l2_mib, l3_mib);
    profile_report_cache_line("nufft workspace", nufft_workspace, l1_mib, l2_mib, l3_mib);
    profile_report_cache_line("estimated total", per_thread_total, l1_mib, l2_mib, l3_mib);
    profile_report_cache_line("shared twiddles", shared_twiddles, l1_mib, l2_mib, l3_mib);

    fprintf(stderr, "\ncache reference: detected cpu0 cache instance L1d=%.2f MiB L2=%.2f MiB L3=%.2f MiB.\n", l1_mib, l2_mib, l3_mib);
    fprintf(stderr, "diagnosis hint: compare physical-core scaling with SMT scaling and v3/v4 phase mix; no hardware counters are used here.\n");
}

#else

#    define PROFILE_START(name) ((void)0)
#    define PROFILE_END(name, phase) ((void)0)
#    define PROFILE_COUNT_TARGET() ((void)0)
#    define PROFILE_FLUSH_THREAD() ((void)0)
#    define profile_report(params) ((void)0)

#endif

#endif
