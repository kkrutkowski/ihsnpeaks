#ifndef TRIG_H
#define TRIG_H

static inline float sin2pif_tls(float x) {
    float f = x - (float)((int)x);
    if (f < 0.0f) f += 1.0f;
    float sign = 1.0f;
    if (f >= 0.5f) {
        sign = -1.0f;
        f -= 0.5f;
    }
    if (f > 0.25f) f = 0.5f - f;
    float f2 = f * f;
    float p = f2 * 39.536706065730207835108712734262f - 76.549782293595742666226937116116f;
    p = p * f2 + 81.601004073261773523492199897936f;
    p = p * f2 - 41.341655031416278077153126232486f;
    p = p * f2 + 6.2831851600894774430188071795666f;
    p *= f;
    return p * sign;
}

static inline float cos2pif_tls(float x) {
    float f = x - (float)((int)x);
    if (f < 0.0f) f += 1.0f;
    if (f > 0.5f) f = 1.0f - f;
    float sign = 1.0f;
    if (f > 0.25f) {
        sign = -1.0f;
        f = 0.5f - f;
    }
    float f2 = f * f;
    float p = f2 * 56.242380464873243259663276802701f - 85.240330322699427859509454517828f;
    p = p * f2 + 64.934590626780991246193352727536f;
    p = p * f2 - 19.739171434702393618770795066531f;
    p = p * f2 + 0.99999995346667013630639784578184f;
    return p * sign;
}

#endif  // TRIG_H
