#ifndef RHOSORT_H
#define RHOSORT_H
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint64_t T; // Changed to 64-bit type
typedef size_t U;

#define LIKELY(X) __builtin_expect(X, 1)
#define RARE(X) __builtin_expect(X, 0)

#ifndef PROF_START
  #define PROF_START(n) (void) 0
  #define PROF_CONT(n)  (void) 0
  #define PROF_END(n)   (void) 0
#endif

// Minimum size to steal from buffer
static const U BLOCK = 16;

// Extract key from the 64-bit value
static inline uint32_t extract_key(T value) {
  return (uint32_t)(value >> 32); // Upper 32 bits
}

// Merge arrays of length l and n-l starting at a, using buffer aux.
static void merge(T *a, U l, U n, T *aux) {
  if (extract_key(a[l - 1]) <= extract_key(a[l])) return;
  if (extract_key(a[n - 1]) < extract_key(a[0]) && l + l == n) {
    T *b = a + l;
    for (U i = 0; i < l; i++) { T t = a[i]; a[i] = b[i]; b[i] = t; }
    return;
  }
  memcpy(aux, a, l * sizeof(T));
  for (U ai = 0, bi = l, i = 0; i < bi; i++) {
    if (bi >= n || extract_key(aux[ai]) <= extract_key(a[bi]))
      a[i] = aux[ai++];
    else
      a[i] = a[bi++];
  }
}

// Merge array x of size n, if units of length block are pre-sorted
static void mergefrom(T *x, U n, U block, T *aux) {
  for (U w = block; w < n; w *= 2)
    for (U i = 0, ww = 2 * w; i < n - w; i += ww)
      merge(x + i, w, n - i < ww ? n - i : ww, aux);
}

// Counting sort of the n values starting at x
static void count(T *x, U n, uint32_t min, U range) {
  U* count = (U*) calloc(range, sizeof(U));
  if (range < n / 8) {
    for (U i = 0; i < n; i++) count[extract_key(x[i]) - min]++;
    for (U i = 0; i < range; i++)
      for (U j = 0; j < count[i]; j++)
        *x++ = (((uint64_t)(min + i)) << 32) | (x[0] & 0xFFFFFFFF);
  } else {
    for (U i = 0; i < n; i++) { count[extract_key(x[i]) - min]++; x[i] = 0; }
    x[0] = ((uint64_t)min) << 32;
    for (U i = 0, s = count[i]; s < n; s += count[++i]) x[s]++;
    { U i = 0;
      for (; i + 4 < n; i += 4) { x[i + 4] += x[i + 3] += x[i + 2] += x[i + 1] += x[i]; }
      for (; i + 1 < n; i++) { x[i + 1] += x[i]; }
    }
  }
  free(count);
}

// The main attraction. Sort array of 64-bit key-value pairs with length n.
void rhsort64(T *array, U n) {
  T *x = array, *xb = x;
  PROF_START(0);
  uint32_t min = extract_key(x[0]), max = min;
  for (U i = 1; i < n; i++) {
    uint32_t e = extract_key(x[i]);
    if (e < min) min = e;
    if (e > max) max = e;
  }
  U r = (U)(max - min) + 1;
  PROF_END(0);
  if (RARE(r / 4 < n)) {
    PROF_START(5); count(x, n, min, r); PROF_END(5); return;
  }

  PROF_START(1);
  uint32_t s = max;
  U sh = 0;
  while (r > 5 * n) { sh++; r >>= 1; }
  U threshold = 2 * BLOCK;
  U sz = r + threshold;

  T* aux = (T*) malloc((sz > n ? sz : n) * sizeof(T));
  for (U i = 0; i < sz; i++) aux[i] = ((uint64_t)s) << 32;
  PROF_END(1);

  PROF_START(2);
  #define POS(E) ((U)((E) - min) >> sh)
  for (U i = 0; i < n; i++) {
    T e = x[i];
    uint32_t key = extract_key(e);
    U j = POS(key);
    T h = aux[j];
    if (LIKELY(extract_key(h) == s)) { aux[j] = e; continue; }

    U j0 = j, f = j;
    do {
      T n = aux[++f];
      int c = key >= extract_key(h);
      j += c;
      aux[f - c] = h;
      h = n;
    } while (extract_key(h) != s);
    aux[j] = e;
    f += 1;

    if (RARE(f - j0 >= threshold)) {
      threshold = BLOCK;
      while (j0 && extract_key(aux[j0 - 1]) != s) j0--;
      T *hj = aux + j0, *hf = aux + f;
      while (hj <= hf - BLOCK) {
        for (U i = 0; i < BLOCK; i++) { xb[i] = hj[i]; hj[i] = ((uint64_t)s) << 32; }
        hj += BLOCK; xb += BLOCK;
      }
      U pr = j0;
      while (hj < hf) {
        e = *hj; *hj++ = ((uint64_t)s) << 32;
        U pp = POS(extract_key(e));
        pr = pp > pr ? pp : pr;
        aux[pr++] = e;
      }
    }
  }
  #undef POS
  PROF_END(2);

  PROF_START(3);
  while (extract_key(aux[--sz]) == s); sz++;
  T *xt = xb;
  {
    static const U u = 8;
    #define WR(I) xt += extract_key(aux[i + I]) != s ? (*xt = aux[i + I], 1) : 0
    U i = 0;
    for (; i < (sz & ~(u - 1)); i += u) { WR(0); WR(1); WR(2); WR(3); WR(4); WR(5); WR(6); WR(7); }
    for (; i < sz; i++) WR(0);
    #undef WR
  }
  while (xt < x + n) *xt++ = ((uint64_t)s) << 32;
  PROF_END(3);

  U l = xb - x;
  if (l) {
    PROF_START(4);
    mergefrom(x, l, BLOCK, aux);
    merge(x, l, n, aux);
    PROF_END(4);
  }

  PROF_CONT(1);
  free(aux);
  PROF_END(1);
}

void rhmergesort64(T *x, U n) {
  static const U size = 1 << 16;
  for (U i = 0; i < n; i += size) rhsort64(x + i, n > i + size ? size : n - i);
  PROF_START(5);
  T* aux = (T*) malloc(n * sizeof(T));
  mergefrom(x, n, size, aux);
  free(aux);
  PROF_END(5);
}

#endif
