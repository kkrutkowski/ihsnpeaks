/*
 * lc_readout.c — Single C translation unit bridging upstream ihsnpeaks
 * readout functions (read_dat, read_fits, linregw_buffer) into lc-qt.
 *
 * All upstream header-only libraries with non-static definitions
 * (fast_convert.h, qfits.h, sds.h, fdist.h) are compiled here exactly once.
 */

#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "lc_readout.h"
#include "utils/readout.h"

#define LC_MAX_READ_BUF (64U * 1024U * 1024U) /* 64 MB hard cap */

/*
 * Allocate a minimal buffer_t suitable for read_dat / read_fits.
 * Only x, y, dy, readBuf, and len are populated; everything else is zeroed.
 */
static buffer_t *alloc_read_buffer(uint32_t maxLen, size_t readBufSize) {
    buffer_t *buf = calloc(1, sizeof(buffer_t));
    if (!buf) return NULL;

    buf->len = maxLen;
    buf->spectrum = sdsempty();
    buf->outBuf = sdsempty();

    buf->x = malloc((size_t)maxLen * sizeof(double));
    buf->y = malloc((size_t)maxLen * sizeof(float));
    buf->dy = malloc((size_t)maxLen * sizeof(float));

    if (!buf->x || !buf->y || !buf->dy) goto fail;

    if (readBufSize > 0) {
        buf->readBuf = malloc(readBufSize);
        if (!buf->readBuf) goto fail;
    }

    return buf;

fail:
    free(buf->x);
    free(buf->y);
    free(buf->dy);
    free(buf->readBuf);
    sdsfree(buf->spectrum);
    sdsfree(buf->outBuf);
    free(buf);
    return NULL;
}

static void free_read_buffer(buffer_t *buf) {
    if (!buf) return;
    free(buf->x);
    free(buf->y);
    free(buf->dy);
    free(buf->readBuf);
    sdsfree(buf->spectrum);
    sdsfree(buf->outBuf);
    free(buf);
}

/*
 * Transfer ownership of buffer arrays into lc_data_t.
 * The buffer_t shell is freed but the x/y/dy arrays live on inside lc_data_t.
 */
static int finalize(buffer_t *buf, lc_data_t *out) {
    if (buf->n == 0) {
        free_read_buffer(buf);
        return -1;
    }
    out->n = buf->n;
    out->x = buf->x;
    out->y = buf->y;
    out->dy = buf->dy;
    out->magnitude = buf->magnitude;
    out->_buf = buf; /* keep shell for detrend */

    /* Detach arrays from buf so free_read_buffer won't double-free */
    buf->x = NULL;
    buf->y = NULL;
    buf->dy = NULL;
    return 0;
}

int lc_load_dat(const char *path, lc_data_t *out) {
    memset(out, 0, sizeof(*out));

    /* Stat the file to determine allocation sizes */
    struct stat st;
    if (stat(path, &st) != 0 || st.st_size <= 0) return -1;

    size_t fileSize = (size_t)st.st_size;
    size_t readBufSize = fileSize + 1; /* +1 for null terminator */
    if (readBufSize > LC_MAX_READ_BUF) readBufSize = LC_MAX_READ_BUF;

    /*
     * Count newlines to get a tight upper-bound on measurement rows.
     * Each data row is terminated by '\n', so newline count >= row count.
     */
    uint32_t maxLen = 0;
    FILE *fp = fopen(path, "rb");
    if (!fp) return -1;
    {
        char chunk[8192];
        size_t rd;
        while ((rd = fread(chunk, 1, sizeof(chunk), fp)) > 0) {
            for (size_t i = 0; i < rd; i++) {
                if (chunk[i] == '\n') maxLen++;
            }
        }
    }
    fclose(fp);
    if (maxLen == 0) maxLen = 1;

    buffer_t *buf = alloc_read_buffer(maxLen, readBufSize);
    if (!buf) return -1;

    read_dat(path, buf);
    return finalize(buf, out);
}

int lc_load_fits(const char *path, lc_data_t *out) {
    memset(out, 0, sizeof(*out));

    /*
     * For FITS, peek at the table to get the row count.
     * We open the table just to read nr, then close and let read_fits
     * re-open it (keeps the wrapper simple and reuses upstream logic).
     */
    qfits_table *table = qfits_table_open(path, 1);
    if (!table) return -1;
    int nr = table->nr;
    qfits_table_close(table);

    if (nr <= 0) return -1;
    uint32_t maxLen = (uint32_t)nr;

    buffer_t *buf = alloc_read_buffer(maxLen, 0);
    if (!buf) return -1;

    read_fits(path, buf);
    return finalize(buf, out);
}

void lc_detrend(lc_data_t *data) {
    if (!data || data->n == 0 || !data->_buf) return;
    buffer_t *buf = (buffer_t *)data->_buf;

    /* Re-attach arrays so linregw_buffer can operate on them */
    buf->x = data->x;
    buf->y = data->y;
    buf->dy = data->dy;
    buf->n = data->n;

    /* Save the original time offset (linregw_buffer zeroes x[0]) */
    double t0 = data->x[0];

    /* Compute mean before detrending */
    double mean = 0.0;
    for (unsigned int i = 0; i < data->n; i++) {
        mean += (double)data->y[i];
    }
    mean /= (double)data->n;

    /* Weighted linear regression detrend (stores intercept in buf->magnitude) */
    linregw_buffer(buf);

    /* Restore original timescale */
    for (unsigned int i = 0; i < data->n; i++) {
        data->x[i] += t0;
    }

    /* Add mean value back so the plot shows physically meaningful levels */
    for (unsigned int i = 0; i < data->n; i++) {
        data->y[i] += (float)mean;
    }

    data->magnitude = buf->magnitude;

    /* Detach again */
    buf->x = NULL;
    buf->y = NULL;
    buf->dy = NULL;
}

void lc_free(lc_data_t *data) {
    if (!data) return;
    if (data->_buf) {
        free_read_buffer((buffer_t *)data->_buf);
        data->_buf = NULL;
    }
    free(data->x);
    free(data->y);
    free(data->dy);
    data->x = NULL;
    data->y = NULL;
    data->dy = NULL;
    data->n = 0;
}
