#ifndef RHOSORT_H
#define RHOSORT_H
#include <stdlib.h>
#include <string.h>

typedef unsigned int T;  // Changed to unsigned type
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

// Merge arrays of length l and n-l starting at a, using buffer aux.
static void merge(T *a, U l, U n, T *aux) {
  if (a[l-1] <= a[l]) return;
  if (a[n-1] < a[0] && l+l == n) {
    T *b = a + l;
    for (U i = 0; i < l; i++) { T t = a[i]; a[i] = b[i]; b[i] = t; }
    return;
  }
  memcpy(aux, a, l * sizeof(T));
  for (U ai = 0, bi = l, i = 0; i < bi; i++) {
    if (bi >= n || aux[ai] <= a[bi])
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
static void count(T *x, U n, T min, U range) {
  U* count = (U*) calloc(range, sizeof(U));
  if (range < n / 8) {
    for (U i = 0; i < n; i++) count[x[i] - min]++;
    for (U i = 0; i < range; i++)
      for (U j = 0; j < count[i]; j++)
        *x++ = min + i;
  } else {
    for (U i = 0; i < n; i++) { count[x[i] - min]++; x[i] = 0; }
    x[0] = min;
    for (U i = 0, s = count[i]; s < n; s += count[++i]) x[s]++;
    { U i = 0;
      for (; i + 4 < n; i += 4) { x[i + 4] += x[i + 3] += x[i + 2] += x[i + 1] += x[i]; }
      for (; i + 1 < n; i++) { x[i + 1] += x[i]; }
    }
  }
  free(count);
}

// The main attraction. Sort array of unsigned ints with length n.
void rhsort32(T *array, U n) {
  T *x = array, *xb = x;
  PROF_START(0);
  T min = x[0], max = min;
  for (U i = 1; i < n; i++) {
    T e = x[i]; if (e < min) min = e; if (e > max) max = e;
  }
  U r = (U)(max - min) + 1;
  PROF_END(0);
  if (RARE(r / 4 < n)) {
    PROF_START(5); count(x, n, min, r); PROF_END(5); return;
  }

  PROF_START(1);
  T s = max;
  U sh = 0;
  while (r > 5 * n) { sh++; r >>= 1; }
  U threshold = 2 * BLOCK;
  U sz = r + threshold;

  T* aux = (T*) malloc((sz > n ? sz : n) * sizeof(T));
  for (U i = 0; i < sz; i++) aux[i] = s;
  PROF_END(1);

  PROF_START(2);
  #define POS(E) ((U)((E) - min) >> sh)
  for (U i = 0; i < n; i++) {
    T e = x[i];
    U j = POS(e);
    T h = aux[j];
    if (LIKELY(h == s)) { aux[j] = e; continue; }

    U j0 = j, f = j;
    do {
      T n = aux[++f];
      int c = e >= h;
      j += c;
      aux[f - c] = h;
      h = n;
    } while (h != s);
    aux[j] = e;
    f += 1;

    if (RARE(f - j0 >= threshold)) {
      threshold = BLOCK;
      while (j0 && aux[j0 - 1] != s) j0--;
      T *hj = aux + j0, *hf = aux + f;
      while (hj <= hf - BLOCK) {
        for (U i = 0; i < BLOCK; i++) { xb[i] = hj[i]; hj[i] = s; }
        hj += BLOCK; xb += BLOCK;
      }
      U pr = j0;
      while (hj < hf) {
        e = *hj; *hj++ = s;
        U pp = POS(e);
        pr = pp > pr ? pp : pr;
        aux[pr++] = e;
      }
    }
  }
  #undef POS
  PROF_END(2);

  PROF_START(3);
  while (aux[--sz] == s); sz++;
  T *xt = xb;
  {
    static const U u = 8;
    #define WR(I) xt += s != (*xt = aux[i + I])
    U i = 0;
    for (; i < (sz & ~(u - 1)); i += u) { WR(0); WR(1); WR(2); WR(3); WR(4); WR(5); WR(6); WR(7); }
    for (; i < sz; i++) WR(0);
    #undef WR
  }
  while (xt < x + n) *xt++ = s;
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

void rhmergesort32(T *x, U n) {
  static const U size = 1 << 16;
  for (U i = 0; i < n; i += size) rhsort32(x + i, n > i + size ? size : n - i);
  PROF_START(5);
  T* aux = (T*) malloc(n * sizeof(T));
  mergefrom(x, n, size, aux);
  free(aux);
  PROF_END(5);
}

#endif
