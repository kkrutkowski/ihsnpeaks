#include <pthread.h>
#include <sched.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>

/***************************
 * kt_for with thread pool *
 ***************************/

struct kt_forpool_t;

typedef struct {
    struct kt_forpool_t *t;
    long i;
    int action;
} kto_worker_t;

typedef struct kt_forpool_t {
    int n_threads, n_pending;
    long n;
    pthread_t *tid;
    kto_worker_t *w;
    void (*func)(void*, long, int);
    void *data;
    pthread_mutex_t mutex;
    pthread_cond_t cv_m, cv_s;
} kt_forpool_t;

static inline long kt_fp_steal_work(kt_forpool_t *t) {
    int i, min_i = -1;
    long k, min = LONG_MAX;
    for (i = 0; i < t->n_threads; ++i)
        if (min > t->w[i].i) min = t->w[i].i, min_i = i;
    k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
    return k >= t->n ? -1 : k;
}

static void *kt_fp_worker(void *data) {
    kto_worker_t *w = (kto_worker_t*)data;
    kt_forpool_t *fp = w->t;
    for (;;) {
        long i;
        int action;
        pthread_mutex_lock(&fp->mutex);
        if (--fp->n_pending == 0)
            pthread_cond_signal(&fp->cv_m);
        w->action = 0;
        while (w->action == 0) pthread_cond_wait(&fp->cv_s, &fp->mutex);
        action = w->action;
        pthread_mutex_unlock(&fp->mutex);
        if (action < 0) break;
        for (;;) { // process jobs allocated to this worker
            i = __sync_fetch_and_add(&w->i, fp->n_threads);
            if (i >= fp->n) break;
            fp->func(fp->data, i, w - fp->w);
        }
        while ((i = kt_fp_steal_work(fp)) >= 0) // steal jobs allocated to other workers
            fp->func(fp->data, i, w - fp->w);
    }
    pthread_exit(0);
}

void *kt_forpool_init(int n_threads) {
    kt_forpool_t *fp;
    int i;
    fp = (kt_forpool_t*)calloc(1, sizeof(kt_forpool_t));
    fp->n_threads = fp->n_pending = n_threads;
    fp->tid = (pthread_t*)calloc(fp->n_threads, sizeof(pthread_t));
    fp->w = (kto_worker_t*)calloc(fp->n_threads, sizeof(kto_worker_t));

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setschedpolicy(&attr, SCHED_BATCH); // Set to SCHED_BATCH

    // Optionally, you can also set a lower priority for batch jobs, but it's not required
    struct sched_param param;
    param.sched_priority = 0; // Default priority for SCHED_BATCH
    pthread_attr_setschedparam(&attr, &param);

    for (i = 0; i < fp->n_threads; ++i) {
        fp->w[i].t = fp;
        // Create each worker thread with custom attributes
        pthread_create(&fp->tid[i], &attr, kt_fp_worker, &fp->w[i]);
    }
    pthread_attr_destroy(&attr); // Clean up the thread attributes

    pthread_mutex_init(&fp->mutex, 0);
    pthread_cond_init(&fp->cv_m, 0);
    pthread_cond_init(&fp->cv_s, 0);

    pthread_mutex_lock(&fp->mutex);
    while (fp->n_pending) pthread_cond_wait(&fp->cv_m, &fp->mutex);
    pthread_mutex_unlock(&fp->mutex);

    return fp;
}

void kt_forpool_destroy(void *_fp) {
    kt_forpool_t *fp = (kt_forpool_t*)_fp;
    int i;
    for (i = 0; i < fp->n_threads; ++i) fp->w[i].action = -1;
    pthread_cond_broadcast(&fp->cv_s);
    for (i = 0; i < fp->n_threads; ++i) pthread_join(fp->tid[i], 0);
    pthread_cond_destroy(&fp->cv_s);
    pthread_cond_destroy(&fp->cv_m);
    pthread_mutex_destroy(&fp->mutex);
    free(fp->w);
    free(fp->tid);
    free(fp);
}

void kt_forpool(void *_fp, void (*func)(void*, long, int), void *data, long n) {
    kt_forpool_t *fp = (kt_forpool_t*)_fp;
    long i;
    if (fp && fp->n_threads > 1) {
        fp->n = n, fp->func = func, fp->data = data, fp->n_pending = fp->n_threads;
        for (i = 0; i < fp->n_threads; ++i) fp->w[i].i = i, fp->w[i].action = 1;
        pthread_mutex_lock(&fp->mutex);
        pthread_cond_broadcast(&fp->cv_s);
        while (fp->n_pending) pthread_cond_wait(&fp->cv_m, &fp->mutex);
        pthread_mutex_unlock(&fp->mutex);
    } else {
        for (i = 0; i < n; ++i) func(data, i, 0);
    }
}
