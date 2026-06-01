#ifndef IHSNPEAKS_COMPAT_H
#define IHSNPEAKS_COMPAT_H

#include <stddef.h>
#include <stdlib.h>

#if !defined(__STDC_VERSION__) || __STDC_VERSION__ < 202311L
#    ifndef constexpr
#        define constexpr const
#    endif
#endif

#if (!defined(__STDC_VERSION__) || __STDC_VERSION__ < 201112L) && !defined(aligned_alloc)
static inline size_t ihsnpeaks_bitceil_size(size_t n) {
    size_t out = sizeof(void *);
    while (out < n) out <<= 1U;
    return out;
}

static inline void *ihsnpeaks_aligned_alloc_fallback(size_t alignment, size_t size) {
    void *ptr = NULL;
    if (alignment == 0) return NULL;
    if (alignment < sizeof(void *)) alignment = sizeof(void *);
    if ((alignment & (alignment - 1U)) != 0) alignment = ihsnpeaks_bitceil_size(alignment);
    if (size % alignment != 0) size = ((size + alignment - 1U) / alignment) * alignment;
    if (posix_memalign(&ptr, alignment, size) != 0) return NULL;
    return ptr;
}

#    define aligned_alloc ihsnpeaks_aligned_alloc_fallback
#endif

#endif
