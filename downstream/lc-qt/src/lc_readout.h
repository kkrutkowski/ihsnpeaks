#ifndef LC_READOUT_H
#define LC_READOUT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    unsigned int n;
    double *x;
    float *y;
    float *dy;
    float magnitude;
    /* internal */
    void *_buf;
} lc_data_t;

int  lc_load_dat(const char *path, lc_data_t *out);
int  lc_load_fits(const char *path, lc_data_t *out);
void lc_detrend(lc_data_t *data);
void lc_free(lc_data_t *data);

#ifdef __cplusplus
}
#endif

#endif /* LC_READOUT_H */
