/* $Id: qfits_byteswap.c,v 1.5 2006/02/17 10:24:52 yjung Exp $
 *
 * This file is part of the ESO QFITS Library
 * Copyright (C) 2001-2004 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Author: yjung $
 * $Date: 2006/02/17 10:24:52 $
 * $Revision: 1.5 $
 * $Name: qfits-6_2_0 $
 */

#ifndef QFITS_DEFS
#define QFITS_DEFS
#define _POSIX_C_SOURCE 200809L

#include <stddef.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <regex.h>
#include <ctype.h>
#include <stdarg.h>
#include <pwd.h>
#include <time.h>
#include <sys/time.h>

#include "config.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/* The following defines the maximum acceptable size for a FITS value */
#define FITSVALSZ                    60

#define QFITS_INVALIDTABLE            0
#define QFITS_BINTABLE                1
#define QFITS_ASCIITABLE            2

/*-----------------------------------------------------------------------------
                                   New types
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Column data type
 */
/*----------------------------------------------------------------------------*/
typedef enum _TFITS_DATA_TYPE_ {
    TFITS_ASCII_TYPE_A,
    TFITS_ASCII_TYPE_D,
    TFITS_ASCII_TYPE_E,
    TFITS_ASCII_TYPE_F,
    TFITS_ASCII_TYPE_I,
    TFITS_BIN_TYPE_A,
    TFITS_BIN_TYPE_B,
    TFITS_BIN_TYPE_C,
    TFITS_BIN_TYPE_D,
    TFITS_BIN_TYPE_E,
    TFITS_BIN_TYPE_I,
    TFITS_BIN_TYPE_J,
    TFITS_BIN_TYPE_L,
    TFITS_BIN_TYPE_M,
    TFITS_BIN_TYPE_P,
    TFITS_BIN_TYPE_X,
    TFITS_BIN_TYPE_UNKNOWN
} tfits_type ;

/*----------------------------------------------------------------------------*/
/**
  @brief    Column object

  This structure contains all information needed to read a column in a table.
  These informations come from the header.
  The qfits_table object contains a list of qfits_col objects.

  This structure has to be created from scratch and filled if one want to
  generate a FITS table.
 */
/*----------------------------------------------------------------------------*/
typedef struct qfits_col
{
    /**
      Number of atoms in one field.
     In ASCII tables, it is the number of characters in the field as defined
     in TFORM%d keyword.
     In BIN tables, it is the number of atoms in each field. For type 'A',
     it is the number of characters. A field with two complex object will
     have atom_nb = 4.
    */
    int            atom_nb ;

    /**
     Number of decimals in a ASCII field.
     This value is always 0 for BIN tables
    */
    int         atom_dec_nb ;

    /**
      Size of one element in bytes. In ASCII tables, atom_size is the size
      of the element once it has been converted in its 'destination' type.
      For example, if "123" is contained in an ASCII table in a column
      defined as I type, atom_nb=3, atom_size=4.
      In ASCII tables:
       - type 'A' : atom_size = atom_nb = number of chars
       - type 'I', 'F' or 'E' : atom_size = 4
       - type 'D' : atom_size = 8
      In BIN tables :
       - type 'A', 'L', 'X', 'B': atom_size = 1
       - type 'I' : atom_size = 2
       - type 'E', 'J', 'C', 'P' : atom_size = 4
       - type 'D', 'M' : atom_size = 8
      In ASCII table, there is one element per field. The size in bytes and
      in number of characters is atom_nb, and the size in bytes after
      conversion of the field is atom_size.
      In BIN tables, the size in bytes of a field is always atom_nb*atom_size.
     */
    int            atom_size ;

    /**
      Type of data in the column as specified in TFORM keyword
      In ASCII tables : TFITS_ASCII_TYPE_* with *=A, I, F, E or D
      In BIN tables : TFITS_BIN_TYPE_* with *=L, X, B, I, J, A, E, D, C, M or P
    */
    tfits_type    atom_type ;

    /** Label of the column */
    char        tlabel[FITSVALSZ] ;

    /** Unit of the data */
    char        tunit[FITSVALSZ] ;

    /** Null value */
    char        nullval[FITSVALSZ] ;

    /** Display format */
    char        tdisp[FITSVALSZ] ;

    /**
      zero and scale are used when the quantity in the field does not
      represent a true physical quantity. Basically, thez should be used
      when they are present: physical_value = zero + scale * field_value
      They are read from TZERO and TSCAL in the header
     */
    int            zero_present ;
    float        zero ;
    int            scale_present ;
    float        scale ;

    /** Offset between the beg. of the table and the beg. of the column.  */
    int            off_beg ;

    /** Flag to know if the column is readable. An empty col is not readable */
    int            readable ;
} qfits_col ;


/*----------------------------------------------------------------------------*/
/**
  @brief    Table object

  This structure contains all information needed to read a FITS table.
  These information come from the header. The object is created by
  qfits_open().

  To read a FITS table, here is a code example:
  @code
  int main(int argc, char* argv[])
  {
      qfits_table     *   table ;
     int                    n_ext ;
    int                    i ;

    // Query the number of extensions
    n_ext = qfits_query_n_ext(argv[1]) ;

    // For each extension
    for (i=0 ; i<n_ext ; i++) {
        // Read all the infos about the current table
        table = qfits_table_open(argv[1], i+1) ;
        // Display the current table
        dump_extension(table, stdout, '|', 1, 1) ;
    }
    return ;
  }
  @endcode
 */
/*----------------------------------------------------------------------------*/
typedef struct qfits_table
{
    /**
        Name of the file the table comes from or it is intended to end to
     */
    char            filename[512] ;
    /**
        Table type.
        Possible values: QFITS_INVALIDTABLE, QFITS_BINTABLE, QFITS_ASCIITABLE
     */
    int                tab_t ;
    /** Width in bytes of the table */
    int                tab_w ;
    /** Number of columns */
    int                nc ;
    /** Number of raws */
    int                nr ;
    /** Array of qfits_col objects */
    qfits_col    *    col ;
} qfits_table ;

/*-----------------------------------------------------------------------------
                        Function ANSI C prototypes
 -----------------------------------------------------------------------------*/
typedef unsigned char byte ;
typedef struct qfits_header {
    void    *   first ;         /* Pointer to list start */
    void    *   last ;          /* Pointer to list end */
    int         n ;             /* Number of cards in list */
    /* For efficient looping internally */
    void    *   current ;
    int         current_idx ;
} qfits_header ;

typedef struct qfitsloader {

    /** Private field to see if structure has been initialized */
    int            _init ;
    /** input: Name of the file you want to read pixels from */
    char    *    filename ;
    /** input: xtension number you want to read */
    int            xtnum ;
    /** input: Index of the plane you want, from 0 to np-1 */
    int            pnum ;
    /** input: Pixel type you want (PTYPE_FLOAT, PTYPE_INT or PTYPE_DOUBLE) */
    int            ptype ;
    /** input: Guarantee file copy or allow file mapping */
    int         map ;

    /** output: Total number of extensions found in file */
    int            exts ;
    /** output: Size in X of the requested plane */
    int            lx ;
    /** output: Size in Y of the requested plane */
    int            ly ;
    /** output: Number of planes present in this extension */
    int            np ;
    /** output: BITPIX for this extension */
    int            bitpix ;
    /** output: Start of the data segment (in bytes) for your request */
    int            seg_start ;
    /** output: Size of the data segment (in bytes) for your request */
    int         seg_size ;
    /** output: BSCALE found for this extension */
    double        bscale ;
    /** output: BZERO found for this extension */
    double        bzero ;

    /** output: Pointer to pixel buffer loaded as integer values */
    int        *    ibuf ;
    /** output: Pointer to pixel buffer loaded as float values */
    float    *    fbuf ;
    /** output: Pointer to pixel buffer loaded as double values */
    double    *    dbuf ;
} qfitsloader ;

typedef struct qfitsdumper {

    /** Name of the file to dump to, "STDOUT" to dump to stdout */
    char     *    filename ;
    /** Number of pixels in the buffer to dump */
    int            npix ;
    /** Buffer type: PTYPE_FLOAT, PTYPE_INT or PTYPE_DOUBLE */
    int            ptype ;

    /** Pointer to input integer pixel buffer */
    int        *    ibuf ;
    /** Pointer to input float pixel buffer */
    float    *    fbuf ;
    /** Pointer to input double pixel buffer */
    double    *    dbuf ;

    /** Requested BITPIX in output FITS file */
    int            out_ptype ;
} qfitsdumper ;

static void byteReverse(unsigned char *buf, unsigned longs) ;
const char * qfits_datamd5(const char *) ;

unsigned short qfits_swap_bytes_16(unsigned short w);
unsigned int qfits_swap_bytes_32(unsigned int dw);
void qfits_swap_bytes(void * p, int s);

static void qfits_cache_activate(void);
static int qfits_is_cached(const char * filename);
static int qfits_cache_add(const char * name);
static void qfits_cache_dump(void) ;
void qfits_card_build(char *, const char *, const char *, const char *) ;

char * qfits_getkey(const char *) ;
char * qfits_getvalue(const char *) ;
char * qfits_getcomment(const char *) ;
char * qfits_expand_keyword(const char *) ;
void qfits_warning(const char *fmt, ...);
void qfits_error(const char *fmt, ...);
char * qfits_get_dir_name(const char *) ;
char * qfits_get_base_name(const char *) ;
char * qfits_get_root_name(const char *) ;
char * qfits_get_ext_name(const char *) ;

int _qfits_isnanf(float) ;
int _qfits_isinff(float) ;
int _qfits_isnand(double) ;
int _qfits_isinfd(double) ;

qfits_header * qfits_header_new(void) ;
qfits_header * qfits_header_default(void) ;
void qfits_header_add(qfits_header *, const char *, const char *, const char *,
        const char *) ;
void qfits_header_add_after(qfits_header *, const char *, const char *,
        const char *, const char *, const char *) ;
void qfits_header_append(qfits_header *, const char *, const char *,
        const char *, const char *) ;
void qfits_header_del(qfits_header *, const char *) ;
int qfits_header_sort(qfits_header **) ;
qfits_header * qfits_header_copy(const qfits_header *) ;
void qfits_header_mod(qfits_header *, const char *, const char *, const char *);
void qfits_header_destroy(qfits_header *) ;
char * qfits_header_getstr(const qfits_header *, const char *) ;
int qfits_header_getitem(const qfits_header *, int, char *, char *, char *,
        char *) ;
char * qfits_header_getcom(const qfits_header *, const char *) ;
int qfits_header_getint(const qfits_header *, const char *, int) ;
double qfits_header_getdouble(const qfits_header *, const char *, double) ;
int qfits_header_getboolean(const qfits_header *, const char *, int) ;
int qfits_header_dump(const qfits_header *, FILE *) ;
int qfitsloader_init(qfitsloader *) ;
int qfits_loadpix(qfitsloader *) ;
int qfits_loadpix_window(qfitsloader *, int, int, int, int) ;
int qfits_pixdump(qfitsdumper *) ;

void * qfits_memory_malloc(size_t, const char *, int) ;
void * qfits_memory_calloc(size_t, size_t, const char *, int) ;
void * qfits_memory_realloc(void *, size_t, const char *, int) ;
void   qfits_memory_free(void *, const char *, int) ;
char * qfits_memory_strdup(const char *, const char *, int) ;
char * qfits_memory_falloc(char *, size_t, size_t *, const char *, int) ;
void qfits_memory_fdealloc(void *, size_t, size_t, const char *, int) ;

static float * qfits_pixin_float(byte *, int, int, double, double) ;
static int * qfits_pixin_int(byte *, int, int, double, double) ;
static double * qfits_pixin_double(byte *, int, int, double, double) ;
static byte * qfits_pixdump_float(float *, int, int) ;
static byte * qfits_pixdump_int(int *, int, int) ;
static byte * qfits_pixdump_double(double *, int, int) ;

void qfits_memory_status(void) ;
int qfits_memory_is_empty(void) ;

qfits_header * qfits_header_read(const char *) ;
qfits_header * qfits_header_read_hdr(const char *) ;
qfits_header * qfits_header_read_hdr_string(const unsigned char *, int) ;
qfits_header * qfits_header_readext(const char *, int) ;
void qfits_zeropad(const char *) ;
int qfits_is_fits(const char *) ;
int qfits_get_hdrinfo(const char *, int, int *, int *) ;
int qfits_get_datinfo(const char *, int, int *, int *) ;

char * qfits_get_datetime_iso8601(void) ;

int qfits_is_table(const char * filename, int xtnum) ;
qfits_header * qfits_table_prim_header_default(void) ;
qfits_header * qfits_table_ext_header_default(const qfits_table *) ;
qfits_table * qfits_table_new(const char *, int, int, int, int) ;
int qfits_col_fill(qfits_col *, int, int, int, tfits_type, const char *,
        const char *, const char *, const char *, int, float, int, float, int) ;
qfits_table * qfits_table_open(const char *, int) ;
void qfits_table_close(qfits_table *) ;
unsigned char * qfits_query_column(const qfits_table *, int, const int *) ;
unsigned char * qfits_query_column_seq(const qfits_table *, int, int, int) ;
void * qfits_query_column_data(const qfits_table *, int, const int *,
        const void *) ;
void * qfits_query_column_seq_data(const qfits_table *, int, int, int,
        const void *) ;
int * qfits_query_column_nulls(const qfits_table *, int, const int *, int *,
        int *);
int qfits_save_table_hdrdump(const void **, const qfits_table *,
        const qfits_header *) ;
int qfits_table_append_xtension(FILE *, const qfits_table *, const void **) ;
int qfits_table_append_xtension_hdr(FILE *, const qfits_table *, const void **,
        const qfits_header *) ;
char * qfits_table_field_to_string(const qfits_table *, int, int, int) ;

static char * qfits_bintable_field_to_string(const qfits_table *, int, int,int);
static char * qfits_asciitable_field_to_string(const qfits_table *, int, int,
        int) ;
static char * qfits_build_format(const qfits_col *) ;
static int qfits_table_append_bin_xtension(FILE *, const qfits_table *,
        const void **) ;
static int qfits_table_append_ascii_xtension(FILE *, const qfits_table *,
        const void **) ;
static int qfits_table_append_data(FILE *, const qfits_table *, const void **) ;
static int qfits_table_get_field_size(int, const qfits_col *) ;
static int qfits_table_interpret_type(const char *, int *, int*, tfits_type *,
        int) ;
static char * qfits_strstrip(const char *);
static double qfits_str2dec(const char *, int) ;
static int qfits_compute_table_width(const qfits_table *) ;

char * qfits_query_hdr(const char *, const char *) ;
char * qfits_query_ext(const char *, const char *, int) ;
int qfits_query_n_ext(const char *) ;
int qfits_query_nplanes(const char *, int) ;
char * qfits_pretty_string(const char *) ;
int qfits_is_boolean(const char *) ;
int qfits_is_int(const char *) ;
int qfits_is_float(const char *) ;
int qfits_is_complex(const char *) ;
int qfits_is_string(const char *) ;
int qfits_get_type(const char *) ;
char * qfits_query_card(const char *, const char *) ;
int qfits_replace_card(const char *, const char *, const char *) ;
const char * qfits_version(void) ;


#endif

#ifndef QFITS_STD_H
#define QFITS_STD_H

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/* FITS header constants */

/** FITS block size */
#define FITS_BLOCK_SIZE     (2880)
/** FITS number of cards per block */
#define FITS_NCARDS         (36)
/** FITS size of each line in bytes */
#define FITS_LINESZ         (80)

#endif

#ifndef MD5_H
#define MD5_H

typedef unsigned int word32 ;

struct MD5Context {
    word32 buf[4];
    word32 bits[2];
    unsigned char in[64];
};

void MD5Init(struct MD5Context *context);
void MD5Update(struct MD5Context *context, unsigned char const *buf,
           unsigned len);
void MD5Final(unsigned char digest[16], struct MD5Context *context);
void MD5Transform(word32 buf[4], word32 const in[16]);

/*
 * This is needed to make RSAREF happy on some MS-DOS compilers.
 */
typedef struct MD5Context MD5_CTX;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup    md5     MD5 message-digest algorithm
 *
 *  The algorithm is due to Ron Rivest.  This code was written by Colin Plumb
 *  in 1993, no copyright is claimed. This code is in the public domain; do
 *  with it what you wish.
 *  Equivalent code is available from RSA Data Security, Inc. This code has
 *  been tested against that, and is equivalent, except that you don't need to
 *  include two pages of legalese with every copy.
 *  To compute the message digest of a chunk of bytes, declare an MD5Context
 *  structure, pass it to MD5Init, call MD5Update as needed on buffers full of
 *  bytes, and then call MD5Final, which will fill a supplied 16-byte array with
 *  the digest.
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

/*
 * Start MD5 accumulation.  Set bit count to 0 and buffer to mysterious
 * initialization constants.
 */
void MD5Init(struct MD5Context *ctx)
{
    ctx->buf[0] = 0x67452301;
    ctx->buf[1] = 0xefcdab89;
    ctx->buf[2] = 0x98badcfe;
    ctx->buf[3] = 0x10325476;

    ctx->bits[0] = 0;
    ctx->bits[1] = 0;
}

/*
 * Update context to reflect the concatenation of another buffer full
 * of bytes.
 */
void MD5Update(struct MD5Context *ctx, unsigned char const *buf, unsigned len)
{
    register word32 t;

    /* Update bitcount */

    t = ctx->bits[0];
    if ((ctx->bits[0] = t + ((word32) len << 3)) < t)
    ctx->bits[1]++;        /* Carry from low to high */
    ctx->bits[1] += len >> 29;

    t = (t >> 3) & 0x3f;    /* Bytes already in shsInfo->data */

    /* Handle any leading odd-sized chunks */

    if (t) {
    unsigned char *p = (unsigned char *) ctx->in + t;

    t = 64 - t;
    if (len < t) {
        memmove(p, buf, len);
        return;
    }
    memmove(p, buf, t);
    byteReverse(ctx->in, 16);
    MD5Transform(ctx->buf, (word32 *) ctx->in);
    buf += t;
    len -= t;
    }
    /* Process data in 64-byte chunks */

    while (len >= 64) {
    memmove(ctx->in, buf, 64);
    byteReverse(ctx->in, 16);
    MD5Transform(ctx->buf, (word32 *) ctx->in);
    buf += 64;
    len -= 64;
    }

    /* Handle any remaining bytes of data. */

    memmove(ctx->in, buf, len);
}

/*
 * Final wrapup - pad to 64-byte boundary with the bit pattern
 * 1 0* (64-bit count of bits processed, MSB-first)
 */
void MD5Final(unsigned char digest[16], struct MD5Context *ctx)
{
    unsigned int count;
    unsigned char *p;

    /* Compute number of bytes mod 64 */
    count = (ctx->bits[0] >> 3) & 0x3F;

    /* Set the first char of padding to 0x80.  This is safe since there is
       always at least one byte free */
    p = ctx->in + count;
    *p++ = 0x80;

    /* Bytes of padding needed to make 64 bytes */
    count = 64 - 1 - count;

    /* Pad out to 56 mod 64 */
    if (count < 8) {
    /* Two lots of padding:  Pad the first block to 64 bytes */
    memset(p, 0, count);
    byteReverse(ctx->in, 16);
    MD5Transform(ctx->buf, (word32 *) ctx->in);

    /* Now fill the next block with 56 bytes */
    memset(ctx->in, 0, 56);
    } else {
    /* Pad block to 56 bytes */
    memset(p, 0, count - 8);
    }
    byteReverse(ctx->in, 14);

    /* Append length in bits and transform */
    ((word32 *) ctx->in)[14] = ctx->bits[0];
    ((word32 *) ctx->in)[15] = ctx->bits[1];

    MD5Transform(ctx->buf, (word32 *) ctx->in);
    byteReverse((unsigned char *) ctx->buf, 4);
    memmove(digest, ctx->buf, 16);
    memset(ctx, 0, sizeof(ctx));    /* In case it's sensitive */
}

/* The four core functions - F1 is optimized somewhat */

/* #define F1(x, y, z) (x & y | ~x & z) */
#define F1(x, y, z) (z ^ (x & (y ^ z)))
#define F2(x, y, z) F1(z, x, y)
#define F3(x, y, z) (x ^ y ^ z)
#define F4(x, y, z) (y ^ (x | ~z))

/* This is the central step in the MD5 algorithm. */
#define MD5STEP(f, w, x, y, z, data, s) \
    ( w += f(x, y, z) + data,  w = w<<s | w>>(32-s),  w += x )

/*
 * The core of the MD5 algorithm, this alters an existing MD5 hash to
 * reflect the addition of 16 longwords of new data.  MD5Update blocks
 * the data and converts bytes into longwords for this routine.
 */
void MD5Transform(word32 buf[4], word32 const in[16])
{
    register word32 a, b, c, d;

    a = buf[0];
    b = buf[1];
    c = buf[2];
    d = buf[3];

    MD5STEP(F1, a, b, c, d, in[0] + 0xd76aa478, 7);
    MD5STEP(F1, d, a, b, c, in[1] + 0xe8c7b756, 12);
    MD5STEP(F1, c, d, a, b, in[2] + 0x242070db, 17);
    MD5STEP(F1, b, c, d, a, in[3] + 0xc1bdceee, 22);
    MD5STEP(F1, a, b, c, d, in[4] + 0xf57c0faf, 7);
    MD5STEP(F1, d, a, b, c, in[5] + 0x4787c62a, 12);
    MD5STEP(F1, c, d, a, b, in[6] + 0xa8304613, 17);
    MD5STEP(F1, b, c, d, a, in[7] + 0xfd469501, 22);
    MD5STEP(F1, a, b, c, d, in[8] + 0x698098d8, 7);
    MD5STEP(F1, d, a, b, c, in[9] + 0x8b44f7af, 12);
    MD5STEP(F1, c, d, a, b, in[10] + 0xffff5bb1, 17);
    MD5STEP(F1, b, c, d, a, in[11] + 0x895cd7be, 22);
    MD5STEP(F1, a, b, c, d, in[12] + 0x6b901122, 7);
    MD5STEP(F1, d, a, b, c, in[13] + 0xfd987193, 12);
    MD5STEP(F1, c, d, a, b, in[14] + 0xa679438e, 17);
    MD5STEP(F1, b, c, d, a, in[15] + 0x49b40821, 22);

    MD5STEP(F2, a, b, c, d, in[1] + 0xf61e2562, 5);
    MD5STEP(F2, d, a, b, c, in[6] + 0xc040b340, 9);
    MD5STEP(F2, c, d, a, b, in[11] + 0x265e5a51, 14);
    MD5STEP(F2, b, c, d, a, in[0] + 0xe9b6c7aa, 20);
    MD5STEP(F2, a, b, c, d, in[5] + 0xd62f105d, 5);
    MD5STEP(F2, d, a, b, c, in[10] + 0x02441453, 9);
    MD5STEP(F2, c, d, a, b, in[15] + 0xd8a1e681, 14);
    MD5STEP(F2, b, c, d, a, in[4] + 0xe7d3fbc8, 20);
    MD5STEP(F2, a, b, c, d, in[9] + 0x21e1cde6, 5);
    MD5STEP(F2, d, a, b, c, in[14] + 0xc33707d6, 9);
    MD5STEP(F2, c, d, a, b, in[3] + 0xf4d50d87, 14);
    MD5STEP(F2, b, c, d, a, in[8] + 0x455a14ed, 20);
    MD5STEP(F2, a, b, c, d, in[13] + 0xa9e3e905, 5);
    MD5STEP(F2, d, a, b, c, in[2] + 0xfcefa3f8, 9);
    MD5STEP(F2, c, d, a, b, in[7] + 0x676f02d9, 14);
    MD5STEP(F2, b, c, d, a, in[12] + 0x8d2a4c8a, 20);

    MD5STEP(F3, a, b, c, d, in[5] + 0xfffa3942, 4);
    MD5STEP(F3, d, a, b, c, in[8] + 0x8771f681, 11);
    MD5STEP(F3, c, d, a, b, in[11] + 0x6d9d6122, 16);
    MD5STEP(F3, b, c, d, a, in[14] + 0xfde5380c, 23);
    MD5STEP(F3, a, b, c, d, in[1] + 0xa4beea44, 4);
    MD5STEP(F3, d, a, b, c, in[4] + 0x4bdecfa9, 11);
    MD5STEP(F3, c, d, a, b, in[7] + 0xf6bb4b60, 16);
    MD5STEP(F3, b, c, d, a, in[10] + 0xbebfbc70, 23);
    MD5STEP(F3, a, b, c, d, in[13] + 0x289b7ec6, 4);
    MD5STEP(F3, d, a, b, c, in[0] + 0xeaa127fa, 11);
    MD5STEP(F3, c, d, a, b, in[3] + 0xd4ef3085, 16);
    MD5STEP(F3, b, c, d, a, in[6] + 0x04881d05, 23);
    MD5STEP(F3, a, b, c, d, in[9] + 0xd9d4d039, 4);
    MD5STEP(F3, d, a, b, c, in[12] + 0xe6db99e5, 11);
    MD5STEP(F3, c, d, a, b, in[15] + 0x1fa27cf8, 16);
    MD5STEP(F3, b, c, d, a, in[2] + 0xc4ac5665, 23);

    MD5STEP(F4, a, b, c, d, in[0] + 0xf4292244, 6);
    MD5STEP(F4, d, a, b, c, in[7] + 0x432aff97, 10);
    MD5STEP(F4, c, d, a, b, in[14] + 0xab9423a7, 15);
    MD5STEP(F4, b, c, d, a, in[5] + 0xfc93a039, 21);
    MD5STEP(F4, a, b, c, d, in[12] + 0x655b59c3, 6);
    MD5STEP(F4, d, a, b, c, in[3] + 0x8f0ccc92, 10);
    MD5STEP(F4, c, d, a, b, in[10] + 0xffeff47d, 15);
    MD5STEP(F4, b, c, d, a, in[1] + 0x85845dd1, 21);
    MD5STEP(F4, a, b, c, d, in[8] + 0x6fa87e4f, 6);
    MD5STEP(F4, d, a, b, c, in[15] + 0xfe2ce6e0, 10);
    MD5STEP(F4, c, d, a, b, in[6] + 0xa3014314, 15);
    MD5STEP(F4, b, c, d, a, in[13] + 0x4e0811a1, 21);
    MD5STEP(F4, a, b, c, d, in[4] + 0xf7537e82, 6);
    MD5STEP(F4, d, a, b, c, in[11] + 0xbd3af235, 10);
    MD5STEP(F4, c, d, a, b, in[2] + 0x2ad7d2bb, 15);
    MD5STEP(F4, b, c, d, a, in[9] + 0xeb86d391, 21);

    buf[0] += a;
    buf[1] += b;
    buf[2] += c;
    buf[3] += d;
}

/**@}*/

/*
 * Note: this code is harmless on little-endian machines.
 */
static void byteReverse(unsigned char *buf, unsigned longs)
{
    word32 t;
    do {
    t = (word32) ((unsigned) buf[3] << 8 | buf[2]) << 16 |
        ((unsigned) buf[1] << 8 | buf[0]);
    *(word32 *) buf = t;
    buf += 4;
    } while (--longs);
}
#endif
#ifndef QFITS_MD5_H
#define QFITS_MD5_H
/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/






/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/* Size of an MD5 hash in bytes (32 bytes are 128 bits) */
#define MD5HASHSZ    32

/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_md5   FITS data block MD5 computation routine
 *
 * This module offers MD5 computation over all data areas of a FITS file.
 *
*/
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                              Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the MD5 hash of data zones in a FITS file.
  @param    filename    Name of the FITS file to examine.
  @return    1 statically allocated character string, or NULL.

  This function expects the name of a FITS file.
  It will compute the MD5 hash on all data blocks in the main data section
  and possibly extensions (including zero-padding blocks if necessary) and
  return it as a string suitable for inclusion into a FITS keyword.

  The returned string is statically allocated inside this function,
  so do not free it or modify it. This function returns NULL in case
  of error.
 */
/*----------------------------------------------------------------------------*/
const char * qfits_datamd5(const char * filename)
{
    static char         datamd5[MD5HASHSZ+1] ;
    struct MD5Context    ctx ;
    unsigned char         digest[16] ;
    FILE             *    in ;
    char                 buf[FITS_BLOCK_SIZE];
    char            *    buf_c ;
    int                    i ;
    int                    in_header ;
    int                    check_fits ;

    /* Check entries */
    if (filename==NULL) return NULL ;
    /* Open input file */
    if ((in=fopen(filename, "r"))==NULL) {
        qfits_error("cannot open file %s", filename);
        return NULL ;
    }
    /* Initialize all variables */
    MD5Init(&ctx);
    in_header=1 ;
    check_fits=0 ;
    /* Loop over input file */
    while (fread(buf, 1, FITS_BLOCK_SIZE, in)==FITS_BLOCK_SIZE) {
        /* First time in the loop: check the file is FITS */
        if (check_fits==0) {
            /* Examine first characters in block */
            if (buf[0]!='S' ||
                buf[1]!='I' ||
                buf[2]!='M' ||
                buf[3]!='P' ||
                buf[4]!='L' ||
                buf[5]!='E' ||
                buf[6]!=' ' ||
                buf[7]!=' ' ||
                buf[8]!='=') {
                qfits_error("file [%s] is not FITS\n", filename);
                fclose(in);
                return NULL ;
            } else {
                check_fits=1 ;
            }
        }
        if (in_header) {
            buf_c = buf ;
            for (i=0 ; i<FITS_NCARDS ; i++) {
                if (buf_c[0]=='E' &&
                    buf_c[1]=='N' &&
                    buf_c[2]=='D' &&
                    buf_c[3]==' ') {
                    in_header=0 ;
                    break ;
                }
                buf_c += FITS_LINESZ ;
            }
        } else {
            /* If current block is a data block */
            /* Try to locate an extension header */
            if (buf[0]=='X' &&
                buf[1]=='T' &&
                buf[2]=='E' &&
                buf[3]=='N' &&
                buf[4]=='S' &&
                buf[5]=='I' &&
                buf[6]=='O' &&
                buf[7]=='N' &&
                buf[8]=='=') {
                in_header=1 ;
                buf_c = buf ;
                for (i=0 ; i<FITS_NCARDS ; i++) {
                    /* Try to find an END marker in this block */
                    if (buf_c[0]=='E' &&
                        buf_c[1]=='N' &&
                        buf_c[2]=='D' &&
                        buf_c[3]==' ') {
                        /* Found END marker in same block as XTENSION */
                        in_header=0;
                        break ;
                    }
                    buf_c += FITS_LINESZ ;
                }
            } else {
                MD5Update(&ctx, (unsigned char *)buf, FITS_BLOCK_SIZE);
            }
        }
    }
    fclose(in);
    if (check_fits==0) {
        /* Never went through the read loop: file is not FITS */
        qfits_error("file [%s] is not FITS", filename);
        return NULL ;
    }
    /* Got to the end of file: summarize */
    MD5Final(digest, &ctx);
    /* Write digest into a string */
    sprintf(datamd5,
    "%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x",
    digest[ 0],
    digest[ 1],
    digest[ 2],
    digest[ 3],
    digest[ 4],
    digest[ 5],
    digest[ 6],
    digest[ 7],
    digest[ 8],
    digest[ 9],
    digest[10],
    digest[11],
    digest[12],
    digest[13],
    digest[14],
    digest[15]);
    return datamd5 ;
}
#endif

#ifndef QFITS_MEMORY_H
#define QFITS_MEMORY_H

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

#define qfits_malloc(s)         qfits_memory_malloc(s,      __FILE__,__LINE__)
#define qfits_calloc(n,s)       qfits_memory_calloc(n,s,    __FILE__,__LINE__)
#define qfits_realloc(p,s)      qfits_memory_realloc(p,s,   __FILE__,__LINE__)
#define qfits_free(p)           qfits_memory_free(p,        __FILE__,__LINE__)
#define qfits_strdup(s)         qfits_memory_strdup(s,      __FILE__,__LINE__)
#define qfits_falloc(f,o,s)     qfits_memory_falloc(f,o,s,  __FILE__,__LINE__)
#define qfits_fdealloc(f,o,s)   qfits_memory_fdealloc(f,o,s,__FILE__,__LINE__)

/*-----------------------------------------------------------------------------
                                Defines
 -----------------------------------------------------------------------------*/

/*
  This symbol sets the debug level for the memory module. Debug levels are
  defined as follows:

  0   no debugging.
  1   add tracing for memory leaks and diagnostics in qfits_memory_status
  2   add lots of debug messages
*/
#ifndef QFITS_MEMORY_DEBUG
#define QFITS_MEMORY_DEBUG       0
#endif

/*
  This symbol defines the level of usage of the memory module.

  0   Use the memory system calls.
  1   Use the memory system calls, but exit if they are not succesfull
  2   Fully use the memory functions
*/
#ifndef QFITS_MEMORY_MODE
#define QFITS_MEMORY_MODE        2
#endif

/* Initial number of entries in memory table */
/* If this number is big, the size of the memory table can become
   problematic.
 */
#ifndef QFITS_MEMORY_MAXPTRS
#define QFITS_MEMORY_MAXPTRS     200003
#endif

/* Identify true RAM memory */
#define MEMTYPE_RAM         'R'
/* Identify swap memory */
#define MEMTYPE_SWAP        'S'
/* Identify memory-mapped file */
#define MEMTYPE_MMAP        'M'

/* Minimal page size in bytes */
#define MEMPAGESZ           2048

/* Size of temporary dir name */
#define TMPDIRNAMESZ        1024

/* Size of temporary file names */
#define TMPFILENAMESZ       1024

/* Size of source file names */
#define SRCFILENAMESZ       64

/* Size of mapped file names */
#define MAPFILENAMESZ       256

/*-----------------------------------------------------------------------------
                                Macros
 -----------------------------------------------------------------------------*/

/* Debug */
#if (QFITS_MEMORY_DEBUG>=2)
#define qfits_mem_debug( code ) { code }
#else
#define qfits_mem_debug( code )
#endif

/* A very simple hash */
#define PTR_HASH(ptr) (((unsigned long int) ptr) % QFITS_MEMORY_MAXPTRS)

/*-----------------------------------------------------------------------------
                        Private variables
 -----------------------------------------------------------------------------*/

/* Initialization flag */
static int  qfits_memory_initialized=0 ;

/* Path to temporary directory */
static char qfits_memory_tmpdirname[TMPDIRNAMESZ] = "." ;

/*----------------------------------------------------------------------------*/
/*
  This table holds a list pointer cells (all the ones allocated so far).
  It is strictly internal to this source file.
 */
/*----------------------------------------------------------------------------*/
static struct {
    /* Number of active cells */
    int                 ncells ;
    /* Total allocated memory in bytes */
    size_t              alloc_total ;
    /* Total allocated RAM in bytes */
    size_t              alloc_ram ;
    /* Total allocated VM in bytes */
    size_t              alloc_swap ;
    /* Peak allocation ever seen for diagnostics */
    size_t              alloc_max ;
    /* Peak number of pointers ever seen for diagnostics */
    int                 max_cells ;

    /* Current number of swap files */
    int                 nswapfiles ;
    /* Registration counter for swap files */
    int                 file_reg ;

    /* Current number of memory-mapped files */
    int                 n_mm_files ;
    /* Current number of mappings derived from files */
    int                 n_mm_mappings ;

#ifdef __linux__
    /* Page size in bytes (Linux only) */
    int                 pagesize ;
    /* Value found for RLIMIT_DATA (Linux only) */
    int                 rlimit_data ;
#endif
} qfits_memory_table ;

/* Various infos about the pointers */
/* List of pointers (outside of cells for efficiency reason) */
static void *   qfits_memory_p_val[QFITS_MEMORY_MAXPTRS] ;
/* Pointed size in bytes */
static size_t   qfits_memory_p_size[QFITS_MEMORY_MAXPTRS] ;
#if (QFITS_MEMORY_DEBUG>=1)
/* Name of the source file where the alloc was requested */
static char *   qfits_memory_p_filename[QFITS_MEMORY_MAXPTRS] ;
/* Line number where the alloc was requested */
static int      qfits_memory_p_lineno[QFITS_MEMORY_MAXPTRS] ;
#endif
/* Memory type: RAM, swap, or mapped file */
static char     qfits_memory_p_memtype[QFITS_MEMORY_MAXPTRS] ;
/* Swap memory only */
/* Swap file ID */
static int      qfits_memory_p_swapfileid[QFITS_MEMORY_MAXPTRS] ;
/* Swap file descriptor */
static int      qfits_memory_p_swapfd[QFITS_MEMORY_MAXPTRS] ;
/* Mapped files only */
/* Name of mapped file */
static char     qfits_memory_p_mm_filename[QFITS_MEMORY_MAXPTRS][MAPFILENAMESZ];
/* Hash of mapped file name for quick search */
static unsigned qfits_memory_p_mm_hash[QFITS_MEMORY_MAXPTRS] ;
/* Reference counter for this pointer */
static int      qfits_memory_p_mm_refcount[QFITS_MEMORY_MAXPTRS] ;

/*-----------------------------------------------------------------------------
                    Private function prototypes
 -----------------------------------------------------------------------------*/

static unsigned qfits_memory_hash(char *) ;
static void qfits_memory_init(void) ;
static void qfits_memory_cleanup(void);
static int qfits_memory_addcell(void*, size_t, const char*, int, char, int,
        int, char*) ;
static int qfits_memory_remcell(int) ;
static void qfits_memory_dumpcell(int, FILE*) ;
static char * qfits_memory_tmpfilename(int) ;
static char * strdup_(const char * str) ;
void qfits_memory_status_(const char *, int) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_memory     POSIX-compatible extended memory handling
 *
 * qfits_memory is a small and efficient module offering memory extension
 * capabitilies to ANSI C programs running on POSIX-compliant systems. It
 * offers several useful features such as memory leak detection, protection for
 * free on NULL or unallocated pointers, and virtually unlimited memory space.
 * qfits_memory requires the @c mmap() system call to be implemented in the
 * local C library to function. This module has been tested on a number of
 * current Unix * flavours and is reported to work fine.
 * The current limitation is the limited number of pointers it can handle at
 * the same time.
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                            Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Allocate memory.
  @param    size        Size (in bytes) to allocate.
  @param    filename    Name of the file where the alloc took place.
  @param    lineno      Line number in the file.
  @return   1 newly allocated pointer.

  This function is a replacement call for malloc. It should never be called
  directly but through a macro instead, as:

  @code
  qfits_memory_malloc(size, __FILE__, __LINE__)
  @endcode
 */
/*----------------------------------------------------------------------------*/
void * qfits_memory_malloc(
        size_t          size,
        const char  *   filename,
        int             lineno)
{
    void    *   ptr ;
    char    *   fname ;
    int         swapfileid ;
    int         swapfd ;
    char        wbuf[MEMPAGESZ] ;
    int         nbufs ;
    int         memtype ;
    int         i ;
    int         pos ;
#ifdef __linux__
    int         p ;
#endif

    /* If QFITS_MEMORY_MODE is 0 or 1, do not use the qfits_memory model  */
    if (QFITS_MEMORY_MODE == 0) return malloc(size);
    else if (QFITS_MEMORY_MODE == 1) {
        ptr = malloc(size);
        if (ptr == NULL) exit(1) ;
        else return ptr ;
    }

    /* Initialize table if needed */
    if (qfits_memory_initialized==0) {
        qfits_memory_init() ;
        qfits_memory_initialized++ ;
    }

    /* Protect the call */
    if (size==0) {
        qfits_mem_debug(
            fprintf(stderr, "qfits_mem: malloc called with 0 size - %s (%d)\n",
                    filename, lineno);
        );
        return NULL ;
    }

    /* Try to allocate in memory */
#ifdef __linux__
    /* Linux does not honor the RLIMIT_DATA limit.
     * The only way to limit the amount of memory taken by
     * a process is to set RLIMIT_AS, which unfortunately also
     * limits down the maximal amount of memory addressable with
     * mmap() calls, making on-the-fly swap space creation useless
     * in this module. To avoid this, the RLIMIT_DATA value
     * is honored here with this test.
     */
    ptr = NULL ;
    if (qfits_memory_table.rlimit_data<1) {
        /* No limit set on RLIMIT_DATA: proceed with malloc */
        ptr = malloc(size);
    } else if (qfits_memory_table.alloc_total+size <=
            (size_t)qfits_memory_table.rlimit_data) {
        /* Next allocation will still be within limits: proceed */
        ptr = malloc(size);
    }
#else
    ptr = malloc(size);
#endif
    if (ptr==NULL) {
        /* No more RAM available: try to allocate private swap */
        qfits_mem_debug(
            fprintf(stderr, "qfits_mem: hit a NULL pointer -- swapping\n");
        );

        /* Create swap file with rights: rw-rw-rw- */
        swapfileid = ++ qfits_memory_table.file_reg ;
        fname = qfits_memory_tmpfilename(swapfileid);
        swapfd = open(fname, O_RDWR | O_CREAT);
        if (swapfd==-1) {
            fprintf(stderr, "qfits_mem: cannot create swap file\n");
            exit(-1);
        }
        fchmod(swapfd, S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH);

        /* Compute number of passes to insert buffer */
        nbufs = size / MEMPAGESZ ;
        if (size % MEMPAGESZ != 0) nbufs ++ ;

        /* Dump empty buffers into file */
        memset(wbuf, 0, MEMPAGESZ);
        for (i=0 ; i<nbufs ; i++) {
            if (write(swapfd, wbuf, MEMPAGESZ)==-1) {
                perror("write");
                fprintf(stderr,
                        "qfits_mem: fatal error: cannot create swapfile\n");
                close(swapfd);
                remove(fname);
                exit(-1);
            }
        }

        /* mmap() the swap file */
        ptr = (void*)mmap(0,
                          nbufs * MEMPAGESZ,
                          PROT_READ | PROT_WRITE,
                          MAP_PRIVATE,
                          swapfd,
                          0);
        if ((char*)ptr == (char*)-1) {
            perror("mmap");
            fprintf(stderr,
                    "qfits_mem: fatal error: mmap failed for swap file\n");
            close(swapfd);
            remove(fname);
            exit(-1);
        }

        qfits_mem_debug(
            fprintf(stderr, "qfits_mem: swap [%s] created for %ld bytes\n",
                fname, (long)size);
        );

        memtype = MEMTYPE_SWAP ;
        qfits_memory_table.alloc_swap += size ;
        qfits_memory_table.nswapfiles ++ ;
    } else {
        /* Memory allocation succeeded */
#ifdef __linux__
        /*
         * On Linux, the returned pointer might not be honored later.
         * To make sure the returned memory is actually usable, it has to
         * be touched. The following will touch one byte every 'pagesize'
         * bytes to make sure all blocks are visited and properly allocated
         * by the OS.
         */
        qfits_mem_debug(
            fprintf(stderr, "qfits_mem: touching memory (Linux)\n");
        );
        for (p=0 ; p<(int)size ; p+=qfits_memory_table.pagesize)
            ((char*)ptr)[p] = 0;
#endif
        swapfd = -1 ;
        swapfileid = -1 ;
        memtype = MEMTYPE_RAM ;
        qfits_memory_table.alloc_ram   += size ;
    }

    /* Print out message in debug mode */
    qfits_mem_debug(
        fprintf(stderr, "qfits_mem: %p alloc(%ld) in %s (%d)\n",
            ptr, (long)size, filename, lineno) ;
    );

    /* Add cell into general table */
    pos = qfits_memory_addcell(  ptr,
                            size,
                            filename,
                            lineno,
                            memtype,
                            swapfileid,
                            swapfd,
                            NULL);
    /* Adjust size */
    qfits_memory_table.alloc_total += size ;
    /* Remember biggest allocated block */
    if (qfits_memory_table.alloc_total > qfits_memory_table.alloc_max)
        qfits_memory_table.alloc_max = qfits_memory_table.alloc_total ;

    /* Insert memory stamp */
    return (void*)ptr ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Allocate memory.
  @param    nmemb       Number of elements to allocate.
  @param    size        Size (in bytes) of each element.
  @param    filename    Name of the file where the alloc took place.
  @param    lineno      Line number in the file.
  @return   1 newly allocated pointer.

  This function is a replacement call for calloc. It should never be called
  directly but through a macro instead, as:

  @code
  qfits_memory_calloc(nmemb, size, __FILE__, __LINE__)
  @endcode
 */
/*----------------------------------------------------------------------------*/
void * qfits_memory_calloc(
        size_t          nmemb,
        size_t          size,
        const char  *   filename,
        int             lineno)
{
    void    *   ptr ;

    /* If QFITS_MEMORY_MODE is 0 or 1, do not use the qfits_memory model  */
    if (QFITS_MEMORY_MODE == 0) return calloc(nmemb, size) ;
    else if (QFITS_MEMORY_MODE == 1) {
        ptr = calloc(nmemb, size) ;
        if (ptr == NULL) exit(1) ;
        else return ptr ;
    }

    ptr = qfits_memory_malloc(nmemb * size, filename, lineno) ;
    return memset(ptr, 0, nmemb * size) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Map a file's contents to memory as a char pointer.
  @param    name        Name of the file to map
  @param    offs        Offset to the first mapped byte in file.
  @param    size        Returned size of the mapped file in bytes.
  @param    srcname     Name of the source file making the call.
  @param    srclin      Line # where the call was made.
  @return   A pointer to char, to be freed using qfits_memory_free().

  This function takes in input the name of a file. It tries to map the file
  into memory and if it succeeds, returns the file's contents as a char pointer.
  It also modifies the input size variable to be the size of the mapped file in
  bytes. This function is normally never directly called but through the
  falloc() macro.

  The offset indicates the starting point for the mapping, i.e. if you are not
  interested in mapping the whole file but only from a given place.

  The returned pointer ptr must be deallocated with qfits_memory_fdealloc(ptr)
 */
/*----------------------------------------------------------------------------*/
char * qfits_memory_falloc(
        char        *   name,
        size_t          offs,
        size_t      *   size,
        const char  *   srcname,
        int             srclin)
{
    unsigned        mm_hash ;
    char        *   ptr ;
    struct stat     sta ;
    int             fd ;
    int             nptrs ;
    int             i ;

    /* If QFITS_MEMORY_MODE is 0 or 1, do not use the qfits_memory model  */
    if ((QFITS_MEMORY_MODE == 0) || (QFITS_MEMORY_MODE == 1)) {

        if (size!=NULL) *size = 0 ;

        /* Check file's existence and compute its size */
        if (stat(name, &sta)==-1) {
            qfits_mem_debug(
                fprintf(stderr, "qfits_mem: cannot stat file %s - %s (%d)\n",
                        name, srcname, srclin);
            );
            if (QFITS_MEMORY_MODE == 0) return NULL ;
            else exit(1) ;
        }
        /* Check offset request does not go past end of file */
        if (offs>=(size_t)sta.st_size) {
            qfits_mem_debug(
                fprintf(stderr,
                    "qfits_mem: falloc offsets larger than file size");
            );
            if (QFITS_MEMORY_MODE == 0) return NULL ;
            else exit(1) ;
        }

        /* Open file */
        if ((fd=open(name, O_RDONLY))==-1) {
            qfits_mem_debug(
                fprintf(stderr, "qfits_mem: cannot open file %s - %s (%d)\n",
                        name, srcname, srclin);
            );
            if (QFITS_MEMORY_MODE == 0) return NULL ;
            else exit(1) ;
        }

        /* Memory-map input file */
        ptr = (char*)mmap(0, sta.st_size,
                PROT_READ | PROT_WRITE, MAP_PRIVATE,fd,0);

        /* Close file */
        close(fd);
        if (ptr == (char*)-1 || ptr==NULL) {
            qfits_mem_debug(
                perror("mmap");
                fprintf(stderr, "qfits_mem: falloc cannot mmap file %s", name);
            );
            if (QFITS_MEMORY_MODE == 0) return NULL ;
            else exit(1) ;
        }

        qfits_mem_debug(
            fprintf(stderr,
                    "qfits_mem: falloc mmap succeeded for [%s] - %s (%d)\n",
                    name, srcname, srclin);
        );

        if (size!=NULL) (*size) = sta.st_size ;

        return ptr + offs ;
    }

    /* Protect the call */
    if (size!=NULL) *size = 0 ;

    /* Initialize table if needed */
    if (qfits_memory_initialized==0) {
        qfits_memory_init() ;
        qfits_memory_initialized++ ;
    }

    if (qfits_memory_table.ncells>0) {
        /* Check if file has already been mapped */
        /* Compute hash for this name */
        mm_hash = qfits_memory_hash(name);
        /* Loop over all memory cells */
        nptrs=0 ;
        for (i=0 ; i<QFITS_MEMORY_MAXPTRS ; i++) {
            if (qfits_memory_p_val[i]!=NULL)
                nptrs++ ;
            if ((qfits_memory_p_val[i]!=NULL) &&
                (qfits_memory_p_mm_filename[i] != NULL) &&
                (qfits_memory_p_mm_hash[i] == mm_hash)) {
                if (!strncmp(qfits_memory_p_mm_filename[i], name,
                             MAPFILENAMESZ)) {
                    /* File already mapped */
                    /* Check offset consistency wrt file size */
                    if (offs >= qfits_memory_p_size[i]) {
                        qfits_mem_debug(
                            fprintf(stderr,
                                "qfits_mem: falloc offset larger than file sz");
                        );
                        return NULL ;
                    }
                    /* Increase reference counter */
                    qfits_memory_p_mm_refcount[i] ++ ;
                    qfits_mem_debug(
                        fprintf(stderr,
                                "qfits_mem: incref on %s (%d mappings)\n",
                                name,
                                qfits_memory_p_mm_refcount[i]);
                    );
                    /* Increase number of mappings */
                    qfits_memory_table.n_mm_mappings ++ ;
                    /* Build up return pointer */
                    ptr = (char*)qfits_memory_p_val[i] + offs ;
                    /* Available size is filesize minus offset */
                    if (size!=NULL) {
                        *size = qfits_memory_p_size[i] - offs ;
                    }
                    /* Return constructed pointer as void * */
                    return (void*)ptr ;
                }
            }
            if (nptrs>=qfits_memory_table.ncells) break ;
        }
    }

    /* First mapping attempt for this file */
    /* Check file's existence and compute its size */
    if (stat(name, &sta)==-1) {
        qfits_mem_debug(
            fprintf(stderr, "qfits_mem: cannot stat file %s - %s (%d)\n",
                    name, srcname, srclin);
        );
        return NULL ;
    }
    /* Check offset request does not go past end of file */
    if (offs>=(size_t)sta.st_size) {
        qfits_mem_debug(
            fprintf(stderr,
                "qfits_mem: falloc offsets larger than file size");
        );
        return NULL ;
    }

    /* Open file */
    if ((fd=open(name, O_RDONLY))==-1) {
        qfits_mem_debug(
            fprintf(stderr, "qfits_mem: cannot open file %s - %s (%d)\n",
                    name, srcname, srclin);
        );
        return NULL ;
    }

    /* Memory-map input file */
    ptr = (char*)mmap(0, sta.st_size, PROT_READ | PROT_WRITE, MAP_PRIVATE,fd,0);

    /* Close file */
    close(fd);
    if (ptr == (char*)-1 || ptr==NULL) {
        qfits_mem_debug(
            perror("mmap");
            fprintf(stderr, "qfits_mem: falloc cannot mmap file %s", name);
        );
        return NULL ;
    }

    qfits_memory_table.n_mm_files ++ ;
    qfits_memory_table.n_mm_mappings ++ ;
    qfits_mem_debug(
        fprintf(stderr,
                "qfits_mem: falloc mmap succeeded for [%s] - %s (%d)\n",
                name, srcname, srclin);
    );

    /* Add cell into general table */
    (void) qfits_memory_addcell((void*)ptr, sta.st_size, srcname, srclin,
                           MEMTYPE_MMAP, -1, -1, name) ;

    if (size!=NULL) (*size) = sta.st_size ;

    return ptr + offs ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Free memory that has been allocated with falloc
  @param    ptr         Pointer to free.
  @param    offs        Offset to the first mapped byte in file.
  @param    size        size to unmap
  @param    filename    Name of the file where the dealloc took place.
  @param    lineno      Line number in the file.
  @return   void
 */
/*----------------------------------------------------------------------------*/
void qfits_memory_fdealloc(
        void        *   ptr,
        size_t          offs,
        size_t          size,
        const char  *   filename,
        int             lineno)
{
    int     i ;
    int     pos ;
    char *  swapname ;
    int     nptrs ;
    int     ii;

    /* Do nothing for a NULL pointer */
    if (ptr==NULL) {
        /* Output a warning */
        fprintf(stderr, "qfits_mem: free requested on NULL ptr -- %s (%d)\n",
                filename, lineno);
        return ;
    }

    /* If QFITS_MEMORY_MODE is 0 or 1, do not use the qfits_memory model  */
    if ((QFITS_MEMORY_MODE == 0) || (QFITS_MEMORY_MODE == 1)) {
        munmap((char*)(ptr)-offs, size) ;
        return ;
    }

    /* Locate pointer in main table */
    nptrs = 0 ;
    pos = -1 ;
    i = PTR_HASH(ptr);
    for (ii=0 ; ii<QFITS_MEMORY_MAXPTRS ; ii++) {
        if (++i == QFITS_MEMORY_MAXPTRS) i = 0;
        if (qfits_memory_p_val[i] == NULL) continue ;
        nptrs++ ;
        if (qfits_memory_p_val[i] == ptr) {
            pos=i ;
            break ;
        }
        if (qfits_memory_p_memtype[i]==MEMTYPE_MMAP) {
            if (((char*)qfits_memory_p_val[i]<=(char*)ptr) &&
                (((char*)qfits_memory_p_val[i] +
                  qfits_memory_p_size[i]) >= (char*)ptr)) {
                pos = i ;
                break ;
            }
        }
        if (nptrs>=qfits_memory_table.ncells) break ;
    }
    if (pos==-1) {
        fprintf(stderr,
                "qfits_mem: %s (%d) free req. on unallocated pointer (%p)\n",
                filename, lineno, ptr);
        /* Pointer sent to system's free() function, maybe it should not? */
        free(ptr);
        return ;
    }

    /* Deallocate pointer */
    switch (qfits_memory_p_memtype[pos]) {
        case MEMTYPE_RAM:
            /* --- RAM pointer */
            /* Free normal memory pointer */
            free(ptr);
            qfits_memory_table.alloc_ram -= qfits_memory_p_size[pos] ;
            break ;
        case MEMTYPE_SWAP:
            /* --- SWAP pointer */
            swapname = qfits_memory_tmpfilename(qfits_memory_p_swapfileid[pos]);
            qfits_mem_debug(
                    fprintf(stderr, "qfits_mem: deallocating swap file [%s]\n",
                        swapname);
            );
            /* Munmap file */
            if (munmap(ptr, qfits_memory_p_size[pos])!=0) {
                qfits_mem_debug( perror("munmap"); );
            }
            /* Close swap file */
            if (close(qfits_memory_p_swapfd[pos])==-1) {
                qfits_mem_debug( perror("close"); );
            }
            /* Remove swap file */
            if (remove(swapname)!=0) {
                qfits_mem_debug( perror("remove"); );
            }
            qfits_memory_table.alloc_swap -= qfits_memory_p_size[pos] ;
            qfits_memory_table.nswapfiles -- ;
            break ;
        case MEMTYPE_MMAP:
            /* --- MEMORY-MAPPED pointer */
            /* Decrease reference count */
            qfits_memory_p_mm_refcount[pos] -- ;
            /* Decrease total number of mappings */
            qfits_memory_table.n_mm_mappings -- ;
            /* Non-null ref count means the file stays mapped */
            if (qfits_memory_p_mm_refcount[pos]>0) {
                qfits_mem_debug(
                    fprintf(stderr, "qfits_mem: decref on %s (%d mappings)\n",
                            qfits_memory_p_mm_filename[pos],
                            qfits_memory_p_mm_refcount[pos]);
                );
                return ;
            }
            /* Ref count reached zero: unmap the file */
            qfits_mem_debug(
                    fprintf(stderr,
                        "qfits_mem: unmapping file %s\n",
                        qfits_memory_p_mm_filename[pos]);
            );
            munmap((char*)qfits_memory_p_val[pos],
                    qfits_memory_p_size[pos]);
            /* Decrease total number of mapped files */
            qfits_memory_table.n_mm_files -- ;
            break ;
        default:
            qfits_mem_debug(
                    fprintf(stderr, "qfits_mem: unknown memory cell type???");
            );
            break ;
    }

    if (qfits_memory_p_memtype[pos]!=MEMTYPE_MMAP) {
        /* Adjust allocated totals */
        qfits_memory_table.alloc_total -= qfits_memory_p_size[pos] ;

        /* Print out message in debug mode */
        qfits_mem_debug(
            fprintf(stderr, "qfits_mem: free(%p) %ld bytes in %s (%d)\n",
                    ptr,
                    (long)qfits_memory_p_size[pos],
                    filename,
                    lineno);
        );
    }
    /* Remove cell from main table */
    qfits_memory_remcell(pos) ;
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Free memory.
  @param    ptr         Pointer to free.
  @param    filename    Name of the file where the dealloc took place.
  @param    lineno      Line number in the file.
  @return   void

  Free the memory associated to a given pointer. Prints out a warning on stderr
  if the requested pointer is NULL or cannot be found in the extended memory
  table.
 */
/*----------------------------------------------------------------------------*/
void qfits_memory_free(
        void        *   ptr,
        const char  *   filename,
        int             lineno)
{
    int     i ;
    int     pos ;
    char *  swapname ;
    int     nptrs ;
    int     ii;

    /* If QFITS_MEMORY_MODE is 0 or 1, do not use the qfits_memory model  */
    if ((QFITS_MEMORY_MODE == 0) || (QFITS_MEMORY_MODE == 1)) {
        free(ptr);
        return ;
    }

    /* Do nothing for a NULL pointer */
    if (ptr==NULL) {
        /* Output a warning */
        fprintf(stderr, "qfits_mem: free requested on NULL ptr -- %s (%d)\n",
                filename, lineno);
        return ;
    }

    /* Locate pointer in main table */
    nptrs = 0 ;
    pos = -1 ;
    i = PTR_HASH(ptr);
    for (ii=0 ; ii<QFITS_MEMORY_MAXPTRS ; ii++) {
        if (++i == QFITS_MEMORY_MAXPTRS) i = 0;
        if (qfits_memory_p_val[i] == NULL) continue ;
        nptrs++ ;
        if (qfits_memory_p_val[i] == ptr) {
            pos=i ;
            break ;
        }
        if (qfits_memory_p_memtype[i]==MEMTYPE_MMAP) {
            if (((char*)qfits_memory_p_val[i]<=(char*)ptr) &&
                (((char*)qfits_memory_p_val[i] +
                  qfits_memory_p_size[i]) >= (char*)ptr)) {
                pos = i ;
                break ;
            }
        }
        if (nptrs>=qfits_memory_table.ncells) break ;
    }
    if (pos==-1) {
        fprintf(stderr,
                "qfits_mem: %s (%d) free requested on unallocated ptr (%p)\n",
                filename, lineno, ptr);
        /* Pointer sent to system's free() function, maybe it should not? */
        free(ptr);
        return ;
    }

    /* Deallocate pointer */
    switch (qfits_memory_p_memtype[pos]) {
        case MEMTYPE_RAM:
            /* --- RAM pointer */
            /* Free normal memory pointer */
            free(ptr);
            qfits_memory_table.alloc_ram -= qfits_memory_p_size[pos] ;
            break ;
        case MEMTYPE_SWAP:
            /* --- SWAP pointer */
            swapname = qfits_memory_tmpfilename(qfits_memory_p_swapfileid[pos]);
            qfits_mem_debug(
                    fprintf(stderr, "qfits_mem: deallocating swap file [%s]\n",
                        swapname);
            );
            /* Munmap file */
            if (munmap(ptr, qfits_memory_p_size[pos])!=0) {
                qfits_mem_debug( perror("munmap"); );
            }
            /* Close swap file */
            if (close(qfits_memory_p_swapfd[pos])==-1) {
                qfits_mem_debug( perror("close"); );
            }
            /* Remove swap file */
            if (remove(swapname)!=0) {
                qfits_mem_debug( perror("remove"); );
            }
            qfits_memory_table.alloc_swap -= qfits_memory_p_size[pos] ;
            qfits_memory_table.nswapfiles -- ;
            break ;
        case MEMTYPE_MMAP:
            /* --- MEMORY-MAPPED pointer */
            /* Decrease reference count */
            qfits_memory_p_mm_refcount[pos] -- ;
            /* Decrease total number of mappings */
            qfits_memory_table.n_mm_mappings -- ;
            /* Non-null ref count means the file stays mapped */
            if (qfits_memory_p_mm_refcount[pos]>0) {
                qfits_mem_debug(
                    fprintf(stderr, "qfits_mem: decref on %s (%d mappings)\n",
                            qfits_memory_p_mm_filename[pos],
                            qfits_memory_p_mm_refcount[pos]);
                );
                return ;
            }
            /* Ref count reached zero: unmap the file */
            qfits_mem_debug(
                    fprintf(stderr,
                        "qfits_mem: unmapping file %s\n",
                        qfits_memory_p_mm_filename[pos]);
            );
            munmap((char*)qfits_memory_p_val[pos],
                    qfits_memory_p_size[pos]);
            /* Decrease total number of mapped files */
            qfits_memory_table.n_mm_files -- ;
            break ;
        default:
            qfits_mem_debug(
                    fprintf(stderr, "qfits_mem: unknown memory cell type???");
            );
            break ;
    }

    if (qfits_memory_p_memtype[pos]!=MEMTYPE_MMAP) {
        /* Adjust allocated totals */
        qfits_memory_table.alloc_total -= qfits_memory_p_size[pos] ;

        /* Print out message in debug mode */
        qfits_mem_debug(
            fprintf(stderr, "qfits_mem: free(%p) %ld bytes in %s (%d)\n",
                    ptr,
                    (long)qfits_memory_p_size[pos],
                    filename,
                    lineno);
        );
    }
    /* Remove cell from main table */
    qfits_memory_remcell(pos) ;
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Re-Allocate memory.
  @param    ptr         Pointer to free.
  @param    size        Size (in bytes) to allocate.
  @param    filename    Name of the file where the alloc took place.
  @param    lineno      Line number in the file.
  @return   1 newly allocated pointer.

  This function is a replacement call for realloc. It should never be called
  directly but through a macro instead, as:

  @code
  qfits_memory_realloc(nmemb, size, __FILE__, __LINE__)
  @endcode
 */
/*----------------------------------------------------------------------------*/
void * qfits_memory_realloc(
        void        *   ptr,
        size_t          size,
        const char  *   filename,
        int             lineno)
{
    void    *   ptr2 ;
    size_t      small_sz ;
    size_t      ptr_sz ;
    int         pos = -1 ;
    int         i ;

    /* If QFITS_MEMORY_MODE is 0 or 1, do not use the qfits_memory model  */
    if (QFITS_MEMORY_MODE == 0) return realloc(ptr, size) ;
    else if (QFITS_MEMORY_MODE == 1) {
        ptr2 = realloc(ptr, size) ;
        if (ptr2 == NULL) exit(1) ;
        else return ptr2 ;
    }

    if (ptr == NULL) return qfits_memory_malloc(size, filename, lineno) ;

    /* Get the pointer size */
    for (i=0 ; i<QFITS_MEMORY_MAXPTRS ; i++) {
        if (qfits_memory_p_val[i] == NULL) continue ;
        if (qfits_memory_p_val[i] == ptr) {
            pos = i ;
            break ;
        }
    }
    if (pos==-1) {
        fprintf(stderr,
            "qfits_mem: %s (%d) realloc requested on unallocated ptr (%p)\n",
            filename, lineno, ptr);
        /* Pointer sent to system's realloc() function, maybe it should not? */
        return realloc(ptr, size) ;
    }
    ptr_sz = qfits_memory_p_size[pos] ;

    /* Compute the smaller size */
    small_sz = size < ptr_sz ? size : ptr_sz ;

    /* Allocate the new pointer */
    ptr2 = qfits_memory_malloc(size, filename, lineno) ;

    /* Copy the common data */
    memcpy(ptr2, ptr, small_sz) ;

    /* Free the passed ptr */
    qfits_memory_free(ptr, filename, lineno) ;

    /* Return  */
    return ptr2 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Duplicate a string using calloc.
  @param    s       String to duplicate.
  @param    filename    Name of the file where the call took place.
  @param    lineno      Line number in the file.
  @return   1 newly allocated character string.

  This function calls in turn calloc to perform the allocation. It should
  never be called directly but only through a macro, like:

  @code
  qfits_memory_strdup(s, __FILE__, __LINE__)
  @endcode

  This function calls qfits_memory_malloc() to do the allocation.
 */
/*----------------------------------------------------------------------------*/
char * qfits_memory_strdup(
        const char  *   s,
        const char  *   filename,
        int             lineno)
{
    char    *   t ;

    /* If QFITS_MEMORY_MODE is 0 or 1, do not use the qfits_memory model  */
    if (QFITS_MEMORY_MODE == 0) return strdup_(s) ;
    else if (QFITS_MEMORY_MODE == 1) {
        t = strdup_(s) ;
        if (t == NULL) exit(1) ;
        else return t ;
    }

    if (s==NULL) return NULL ;
    t = qfits_memory_malloc(1+strlen(s), filename, lineno);
    return strcpy(t, s);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Display memory status information.
  @return   void

  This function is meant for debugging purposes, but it is recommended to
  call it at the end of every executable making use of the extended memory
  features.
 */
/*----------------------------------------------------------------------------*/
void qfits_memory_status(void)
{
    int     i ;

    /* If QFITS_MEMORY_MODE is 0 or 1, do not use the qfits_memory model  */
    if ((QFITS_MEMORY_MODE == 0) || (QFITS_MEMORY_MODE == 1)) return ;

#if (QFITS_MEMORY_DEBUG>=1)
    fprintf(stderr, "#----- memory diagnostics -----\n") ;

    fprintf(stderr,
            "#- Peak memory usage\n"
            "ALL_maxalloc_kb     %ld\n"
            "ALL_maxpointers     %d\n",
            (long)(qfits_memory_table.alloc_max/1024),
            qfits_memory_table.max_cells);
    fprintf(stderr,
            "#- Local implementation\n"
            "TAB_ptrs            %d\n"
            "TAB_size            %u bytes\n",
            QFITS_MEMORY_MAXPTRS,
            (unsigned)sizeof(qfits_memory_table));
#ifdef __linux__
    fprintf(stderr,
            "#- Linux specific\n"
            "LINUX_pagesize      %d bytes\n"
            "LINUX_RLIMIT_DATA   %d kb\n",
            qfits_memory_table.pagesize,
            qfits_memory_table.rlimit_data);
#endif
#endif

    if (qfits_memory_table.ncells<1) return ;
    fprintf(stderr, "#----- memory diagnostics -----\n") ;

    fprintf(stderr,
            "#- ALL status\n"
            "ALL_npointers       %d\n"
            "ALL_size            %ld\n"
            "ALL_maxalloc_kb     %ld\n"
            "ALL_maxpointers     %d\n",
            qfits_memory_table.ncells,
            (long)qfits_memory_table.alloc_total,
            (long)(qfits_memory_table.alloc_max/1024),
            qfits_memory_table.max_cells);

    if (qfits_memory_table.alloc_ram > 0) {
        fprintf(stderr,
                "#- RAM status\n"
                "RAM_alloc           %ld\n",
                (long)qfits_memory_table.alloc_ram);
    }
    if (qfits_memory_table.alloc_swap > 0) {
        fprintf(stderr,
                "#- SWP status\n"
                "SWP_alloc           %ld\n"
                "SWP_files           %d\n",
                (long)qfits_memory_table.alloc_swap,
                qfits_memory_table.nswapfiles);
    }

    if (qfits_memory_table.n_mm_files>0) {
        fprintf(stderr,
                "#- MAP status\n"
                "MAP_files           %d\n"
                "MAP_mappings        %d\n",
                qfits_memory_table.n_mm_files,
                qfits_memory_table.n_mm_mappings);
    }

    fprintf(stderr, "#- pointer details\n");
    for (i=0 ; i<QFITS_MEMORY_MAXPTRS; i++) {
        qfits_memory_dumpcell(i, stderr);
    }
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Tell if there is still some memory allocated
  @return   1 if the memory table is tempty, 0 if no,
            -1 if the memory model is off
 */
/*----------------------------------------------------------------------------*/
int qfits_memory_is_empty(void)
{
    if ((QFITS_MEMORY_MODE == 0) || (QFITS_MEMORY_MODE == 1)) return -1 ;
    if (qfits_memory_table.ncells<1) return 1 ;
    else return 0 ;
}

/**@}*/

/*
 * posted to comp.sys.next.programmer:
 *
 *
 * From: moser@ifor.math.ethz.ch (Dominik Moser,CLV A4,2 40 19,720 49 89)
 * Subject: Re: Compile problems (pgp 2.6.3i)
 * Date: 10 Jul 1996 06:50:42 GMT
 * Organization: Swiss Federal Institute of Technology (ETHZ)
 * References: <4rrhvj$6fr@bagan.srce.hr>
 * Message-ID: <4rvjs2$6oh@elna.ethz.ch>
 *
 * Most systems don't have this (yet)

 */
static char * strdup_(const char * str)
{
    char    *   p ;

    if ((p = malloc(strlen(str)+1)) == NULL)
    return((char *) NULL);

    (void) strcpy(p, str);

    return(p);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Hash a string to an unsigned value.
  @param    key     String to hash
  @return   1 unsigned value as a hash for the given string.

  This hash function has been taken from an Article in Dr Dobbs Journal. This
  is normally a collision-free function, distributing keys evenly. The key is
  stored anyway in the struct so that collision can be avoided by comparing the
  key itself in last resort.
 */
/*----------------------------------------------------------------------------*/
static unsigned qfits_memory_hash(char * key)
{
    int         len ;
    unsigned    hash ;
    int         i ;

    len = strlen(key);
    for (hash=0, i=0 ; i<len ; i++) {
        hash += (unsigned)key[i] ;
        hash += (hash<<10);
        hash ^= (hash>>6) ;
    }
    hash += (hash <<3);
    hash ^= (hash >>11);
    hash += (hash <<15);
    return hash ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Initialize extended memory features.
  @return   void

  This function is implicitly called by the first malloc() or calloc() or
  strdup_() execution. It allocates a minimal number of memory cells into
  the global extended memory table. It also install atexit routines the first
  time it is called, and increases the number of possible descriptors to the
  maximum.
 */
/*----------------------------------------------------------------------------*/
static void qfits_memory_init(void)
{
    struct rlimit rlim ;

    qfits_mem_debug(
        fprintf(stderr,
                "qfits_mem: initializing main table size=%d ptrs (%ld bytes)\n",
                QFITS_MEMORY_MAXPTRS,
                (long)sizeof(qfits_memory_table));
    );
    /* Initialize memory table */
    memset(&qfits_memory_table, 0, sizeof(qfits_memory_table));

    /* Install cleanup routine at exit */
    atexit(qfits_memory_cleanup);

    /* Increase number of descriptors to maximum */
    getrlimit(RLIMIT_NOFILE, &rlim) ;
    qfits_mem_debug(
        fprintf(stderr, "qfits_mem: increasing from %ld to %ld file handles\n",
                (long)rlim.rlim_cur,
                (long)rlim.rlim_max);
    );
    rlim.rlim_cur = rlim.rlim_max ;
    setrlimit(RLIMIT_NOFILE, &rlim) ;

#ifdef __linux__
    /* Get RLIMIT_DATA on Linux */
    getrlimit(RLIMIT_DATA, &rlim);
    qfits_memory_table.rlimit_data = rlim.rlim_cur ;
    qfits_mem_debug(
        fprintf(stderr, "qfits_mem: got RLIMIT_DATA=%d\n",
                qfits_memory_table.rlimit_data);
    );
    /* Get page size on Linux */
    qfits_memory_table.pagesize = sysconf(_SC_PAGESIZE);

#endif
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Removes all swap files.
  @return   void

  This routine will delete all swap files from the temporary area.
 */
/*----------------------------------------------------------------------------*/
static void qfits_memory_cleanup(void)
{
    int     reg ;

    if (qfits_memory_table.file_reg>0) {
        qfits_mem_debug(
            fprintf(stderr, "qfits_mem: cleaning up swap files... ");
        );
        /*
         * Call remove() on all possible VM files. If the file exists, it
         * is effectively removed. It it does not, ignore the error.
         * This is not the cleanest way of doing it, but this function is
         * meant to be called also in cases of emergency (e.g. segfault),
         * so it should not rely on a correct memory table.
         */
        for (reg=0 ; reg<qfits_memory_table.file_reg ; reg++) {
            remove(qfits_memory_tmpfilename(reg+1));
        }
        qfits_mem_debug(
            fprintf(stderr, "qfits_mem: done cleaning swap files\n");
        );
    }
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Add allocation cell to qfits_memory_table.
  @param    pointer         Pointer value.
  @param    size            Pointer size.
  @param    filename        Name of the C source file where alloc was done.
  @param    lineno          Line # where the allocation took place.
  @param    memtype         Memory type: RAM or SWAP.
  @param    swapfileid      Associated swap file ID (if any).
  @param    swapfd          Associated swap file descriptor (if any).
  @param    mm_filename     Mapped file name (if any)
  @return   the index in the qfits_memory_table of the added cell

  Add a memory cell in the xtended memory table to register that a new
  allocation took place.
  This call is not protected against illegal parameter values, so make sure
  the passed values are correct!
 */
/*----------------------------------------------------------------------------*/
static int qfits_memory_addcell(
        void        *   pointer,
        size_t          size,
        const char  *   filename,
        int             lineno,
        char            memtype,
        int             swapfileid,
        int             swapfd,
        char        *   mm_filename)
{
    int pos, ii ;

    /* Check there is still some space left */
    if (qfits_memory_table.ncells >= QFITS_MEMORY_MAXPTRS) {
        fprintf(stderr, "fatal qfits_memory error: reached max pointers (%d)\n",
                QFITS_MEMORY_MAXPTRS);
        exit(-1);
    }
    /* Find an available slot */
    pos = PTR_HASH(pointer);
    for (ii = 0 ; ii<QFITS_MEMORY_MAXPTRS ; ii++) {
        if (++pos == QFITS_MEMORY_MAXPTRS) pos = 0;
        if (qfits_memory_p_val[pos] == NULL) break ;
    }
    qfits_mem_debug(
            fprintf(stderr, "qfits_mem: freecell found at pos %d\n", pos);
            );

    /* Store information */
    qfits_memory_p_val[pos] = pointer ;
    qfits_memory_p_size[pos] = size ;

    /* Filename and line number */
#if (QFITS_MEMORY_DEBUG>=1)
    qfits_memory_p_filename[pos] = filename ;
    qfits_memory_p_lineno[pos] = lineno ;
#endif

    qfits_memory_p_memtype[pos] = memtype ;
    qfits_memory_p_swapfileid[pos] = swapfileid ;
    qfits_memory_p_swapfd[pos] = swapfd ;

    if (mm_filename!=NULL) {
        strncpy(qfits_memory_p_mm_filename[pos], mm_filename,MAPFILENAMESZ);
        qfits_memory_p_mm_hash[pos] = qfits_memory_hash(mm_filename);
        qfits_memory_p_mm_refcount[pos] = 1 ;
    } else {
        qfits_memory_p_mm_filename[pos][0] = 0 ;
        qfits_memory_p_mm_hash[pos] = 0 ;
        qfits_memory_p_mm_refcount[pos] = 0 ;
    }
    qfits_memory_table.ncells ++ ;
    if (qfits_memory_table.ncells > qfits_memory_table.max_cells)
        qfits_memory_table.max_cells = qfits_memory_table.ncells ;
    return pos ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Remove a memory cell from the xtended memory table.
  @param    pos     Position of the pointer in the table.
  @return   int 0 if Ok, -1 if error occurred.

  Remove the specified cell in qfits_memory_table.
  This call is not protected against illegal parameter values, so make sure
  the passed values are correct!
 */
/*----------------------------------------------------------------------------*/
static int qfits_memory_remcell(int pos)
{
    qfits_mem_debug(
        fprintf(stderr, "qfits_mem: removing cell from pos %d (cached)\n", pos);
    );
    /* Set pointer to NULL */
    qfits_memory_p_val[pos] = NULL ;
    /* Decrement number of allocated pointers */
    qfits_memory_table.ncells -- ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Dump a memory cell to an open file pointer.
  @param    cell    Cell to dump.
  @param    out     Open file pointer to dump to.
  @return   void

  This function is meant for debugging purposes only. It takes in input a
  pointer to a memory cell and dumps it to the requested file pointer (it
  is Ok to provide stdout or stderr as file pointers). If the passed
  position is invalid or the table pointer is NULL, this function returns
  immediately.
 */
/*----------------------------------------------------------------------------*/
static void qfits_memory_dumpcell(
        int         pos,
        FILE    *   out)
{
    if (pos<0 || pos>=QFITS_MEMORY_MAXPTRS) return ;
    if (qfits_memory_p_val[pos]==NULL) return ;

    if (qfits_memory_p_memtype[pos] == MEMTYPE_MMAP) {
#if (QFITS_MEMORY_DEBUG>=1)
        fprintf(out,
            "M(%p) - %s (%d) maps [%s] for %ld bytes",
            qfits_memory_p_val[pos],
            qfits_memory_p_filename[pos],
            qfits_memory_p_lineno[pos],
            qfits_memory_p_mm_filename[pos],
            (long)qfits_memory_p_size[pos]);
#else
        fprintf(out,
            "M(%p) maps [%s] for %ld bytes",
            qfits_memory_p_val[pos],
            qfits_memory_p_mm_filename[pos],
            (long)qfits_memory_p_size[pos]);
#endif
    } else {
#if (QFITS_MEMORY_DEBUG>=1)
        fprintf(out, "%c(%p) - %s (%d) for %ld bytes",
            qfits_memory_p_memtype[pos],
            qfits_memory_p_val[pos],
            qfits_memory_p_filename[pos],
            qfits_memory_p_lineno[pos],
            (long)qfits_memory_p_size[pos]);
#else
        fprintf(out, "%c(%p) for %ld bytes",
            qfits_memory_p_memtype[pos],
            qfits_memory_p_val[pos],
            (long)qfits_memory_p_size[pos]);
#endif
    }
    if (qfits_memory_p_memtype[pos]==MEMTYPE_SWAP) {
        fprintf(out, " swf[%s][%d]",
                qfits_memory_tmpfilename(qfits_memory_p_swapfileid[pos]),
                qfits_memory_p_swapfd[pos]);
    }
    fprintf(out, "\n");
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute filename associated to a temporary file ID.
  @param    reg     Registration number of temporary file name.
  @return   pointer to statically allocated char string.

  This function computes the valid file name associated to a temporary file
  ID. It computes the result, stores it in an internal static string and
  returns a pointer to it.
 */
/*----------------------------------------------------------------------------*/
static char * qfits_memory_tmpfilename(int reg)
{
    static char qfits_mem_tmpfilename[TMPFILENAMESZ] ;
    /* Create file name using tmp directory as a base */
    sprintf(qfits_mem_tmpfilename, "%s/vmswap_%05ld_%05x",
            qfits_memory_tmpdirname, (long)getpid(), reg) ;
    return qfits_mem_tmpfilename ;
}
#endif

#ifndef QFITS_BYTESWAP_H
#define QFITS_BYTESWAP_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/





/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_byteswap  Low-level byte-swapping routines
 *
 *  This module offers access to byte-swapping routines.
 *  Generic routines are offered that should work everywhere.
 *  Assembler is also included for x86 architectures, and dedicated
 *  assembler calls for processors > 386.
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Swap a 16-bit number
  @param    w A 16-bit (short) number to byte-swap.
  @return   The swapped version of w, w is untouched.

  This function swaps a 16-bit number, returned the swapped value without
  modifying the passed argument. Assembler included for x86 architectures.
 */
/*----------------------------------------------------------------------------*/
unsigned short qfits_swap_bytes_16(unsigned short w)
{
#ifdef CPU_X86
    __asm("xchgb %b0,%h0" :
            "=q" (w) :
            "0" (w));
    return w ;
#else
    return (((w) & 0x00ff) << 8 | ((w) & 0xff00) >> 8);
#endif
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Swap a 32-bit number
  @param    dw A 32-bit (long) number to byte-swap.
  @return   The swapped version of dw, dw is untouched.

  This function swaps a 32-bit number, returned the swapped value without
  modifying the passed argument. Assembler included for x86 architectures
  and optimized for processors above 386.
 */
/*----------------------------------------------------------------------------*/
unsigned int qfits_swap_bytes_32(unsigned int dw)
{
#ifdef CPU_X86
#if CPU_X86 > 386
    __asm("bswap   %0":
            "=r" (dw)   :
#else
    __asm("xchgb   %b0,%h0\n"
        " rorl    $16,%0\n"
        " xchgb   %b0,%h0":
        "=q" (dw)      :
#endif
        "0" (dw));
    return dw ;
#else
    return ((((dw) & 0xff000000) >> 24) | (((dw) & 0x00ff0000) >>  8) |
            (((dw) & 0x0000ff00) <<  8) | (((dw) & 0x000000ff) << 24));
#endif
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Swaps bytes in a variable of given size
  @param    p pointer to void (generic pointer)
  @param    s size of the element to swap, pointed to by p
  @return    void

  This byte-swapper is portable and works for any even variable size.
  It is not truly the most efficient ever, but does its job fine
  everywhere this file compiles.
 */
/*----------------------------------------------------------------------------*/
void qfits_swap_bytes(void * p, int s)
{
    unsigned char tmp, *a, *b ;

    a = (unsigned char*)p ;
    b = a + s ;

    while (a<b) {
        tmp = *a ;
        *a++ = *--b ;
        *b = tmp ;
    }
}
#endif
/**@}*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifndef QFITS_CACHE_H
#define QFITS_CACHE_H








/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/** Query the number of extensions */
#define QFITS_QUERY_N_EXT        (1<<30)
/** Query the offset to header start */
#define QFITS_QUERY_HDR_START    (1<<29)
/** Query the offset to data start */
#define QFITS_QUERY_DAT_START    (1<<28)
/** Query header size in bytes */
#define QFITS_QUERY_HDR_SIZE    (1<<27)
/** Query data size in bytes */
#define QFITS_QUERY_DAT_SIZE    (1<<26)

void qfits_cache_purge(void) ;
int qfits_query(const char *, int) ;

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/* Define this symbol to get debug symbols -- not recommended! */
#define QFITS_CACHE_DEBUG    0
#if QFITS_CACHE_DEBUG
#define qdebug( code ) { code }
#else
#define qdebug( code )
#endif

/*
 * Cache size:
 * Maximum number of FITS file informations stored in the cache.
 */
#define QFITS_CACHESZ        128

/*
 * This static definition declares the maximum possible number of
 * extensions in a FITS file. It only has effects in the qfits_cache_add
 * function where a table is statically allocated for efficiency reasons.
 * If the number of extensions grows over this limit, change the value of
 * this constant. If the number of extensions is a priori unknown but can
 * grow much larger than a predictable value, the best solution is to
 * implement a dynamic memory allocation in qfits_cache_add.
 */
#define QFITS_MAX_EXTS            10192

/*-----------------------------------------------------------------------------
                                   New types
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*
  This structure stores all informations about a given FITS file.
  It is strictly internal to this module.
 */
/*----------------------------------------------------------------------------*/
typedef struct _qfits_cache_cell_ {
    char    *    name ;    /* File name     */
    ino_t       inode ; /* Inode */
    time_t        mtime;  /* Last modification date */
    int            filesize; /* File size in bytes */
    time_t        ctime;  /* Last modification date */

    int            exts ;    /* # of extensions in file */

    int        *    ohdr ;    /* Offsets to headers */
    int        *    shdr ;    /* Header sizes */
    int        *    data ;    /* Offsets to data */
    int        *    dsiz ;    /* Data sizes */

    int         fsize ; /* File size in blocks (2880 bytes) */
} qfits_cache_cell ;

static qfits_cache_cell qfits_cache[QFITS_CACHESZ] ;
static int qfits_cache_last = -1 ;
static int qfits_cache_entries = 0 ;
static int qfits_cache_init = 0 ;



/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_cache     FITS caching capabilities
 *
 * This modules implements a cache for FITS access routines.
 * The first time a FITS file is seen by the library, all corresponding
 * pointers are cached here. This speeds up multiple accesses to large
 * files by magnitudes.
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Purge the qfits cache.
  @return    void

  This function is useful for programs running for a long period,
  to clean up the cache. Ideally in a daemon, it should be called
  by a timer at regular intervals. Notice that since the cache is
  fairly small, you should not need to care too much about this.
 */
/*----------------------------------------------------------------------------*/
void qfits_cache_purge(void)
{
    int    i ;

    qdebug(
        printf("qfits: purging cache...\n");
    );

    for (i=0 ; i<QFITS_CACHESZ; i++) {
        if (qfits_cache[i].name!=NULL) {
            qfits_free(qfits_cache[i].name);
            qfits_cache[i].name = NULL ;
            qfits_free(qfits_cache[i].ohdr);
            qfits_free(qfits_cache[i].data);
            qfits_free(qfits_cache[i].shdr);
            qfits_free(qfits_cache[i].dsiz);
            qfits_cache_entries -- ;
        }
    }
    if (qfits_cache_entries!=0) {
        qdebug(
            printf("qfits: internal error in cache consistency\n");
        );
        exit(-1);
    }
    qfits_cache_last = -1 ;
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Query a FITS file offset from the cache.
  @param    filename    Name of the file to examine.
  @param    what        What should be queried (see below).
  @return    an integer offset, or -1 if an error occurred.

  This function queries the cache for FITS offset information. If the
  requested file name has never been seen before, it is completely parsed
  to extract all offset informations, which are then stored in the cache.
  The next query will get the informations from the cache, avoiding
  a complete re-parsing of the file. This is especially useful for large
  FITS files with lots of extensions, because querying the extensions
  is an expensive operation.

  This operation has side-effects: the cache is an automatically
  allocated structure in memory, that can only grow. Every request
  on a new FITS file will make it grow. The structure is pretty
  light-weight in memory, but nonetheless this is an issue for daemon-type
  programs which must run over long periods. The solution is to clean
  the cache using qfits_cache_purge() at regular intervals. This is left
  to the user of this library.

  To request information about a FITS file, you must pass an integer
  built from the following symbols:

  - @c QFITS_QUERY_N_EXT
  - @c QFITS_QUERY_HDR_START
  - @c QFITS_QUERY_DAT_START
  - @c QFITS_QUERY_HDR_SIZE
  - @c QFITS_QUERY_DAT_SIZE

  Querying the number of extensions present in a file is done
  simply with:

  @code
  next = qfits_query(filename, QFITS_QUERY_N_EXT);
  @endcode

  Querying the offset to the i-th extension header is done with:

  @code
  off = qfits_query(filename, QFITS_QUERY_HDR_START | i);
  @endcode

  i.e. you must OR (|) the extension number with the
  @c QFITS_QUERY_HDR_START symbol. Requesting offsets to extension data is
  done in the same way:

  @code
  off = qfits_query(filename, QFITS_QUERY_DAT_START | i);
  @endcode

  Notice that extension 0 is the main header and main data part
  of the FITS file.
 */
/*----------------------------------------------------------------------------*/
int qfits_query(const char * filename, int what)
{
    int    rank ;
    int    which ;
    int    answer ;

    qdebug(
        printf("qfits: cache req %s\n", filename);
    );
    if ((rank=qfits_is_cached(filename))==-1) {
        rank = qfits_cache_add(filename);
    }
    if (rank==-1) {
        qdebug(
            printf("qfits: error adding %s to cache\n", filename);
        );
        return -1 ;
    }

    /* See what was requested */
    answer=-1 ;
    if (what & QFITS_QUERY_N_EXT) {
        answer = qfits_cache[rank].exts ;
        qdebug(
            printf("qfits: query n_exts\n");
            printf("qfits: -> %d\n", answer);
        );
    } else if (what & QFITS_QUERY_HDR_START) {
        which = what & (~QFITS_QUERY_HDR_START);
        if (which>=0 && which<=qfits_cache[rank].exts) {
            answer = qfits_cache[rank].ohdr[which] * FITS_BLOCK_SIZE ;
        }
        qdebug(
            printf("qfits: query offset to header %d\n", which);
            printf("qfits: -> %d (%d bytes)\n", answer/2880, answer);
        );
    } else if (what & QFITS_QUERY_DAT_START) {
        which = what & (~QFITS_QUERY_DAT_START);
        if (which>=0 && which<=qfits_cache[rank].exts) {
            answer = qfits_cache[rank].data[which] * FITS_BLOCK_SIZE ;
        }
        qdebug(
            printf("qfits: query offset to data %d\n", which);
            printf("qfits: -> %d (%d bytes)\n", answer/2880, answer);
        );
    } else if (what & QFITS_QUERY_HDR_SIZE) {
        which = what & (~QFITS_QUERY_HDR_SIZE);
        if (which>=0 && which<=qfits_cache[rank].exts) {
            answer = qfits_cache[rank].shdr[which] * FITS_BLOCK_SIZE ;
        }
        qdebug(
            printf("qfits: query sizeof header %d\n", which);
            printf("qfits: -> %d (%d bytes)\n", answer/2880, answer);
        );
    } else if (what & QFITS_QUERY_DAT_SIZE) {
        which = what & (~QFITS_QUERY_DAT_SIZE);
        if (which>=0 && which<=qfits_cache[rank].exts) {
            answer = qfits_cache[rank].dsiz[which] * FITS_BLOCK_SIZE ;
        }
        qdebug(
            printf("qfits: query sizeof data %d\n", which);
            printf("qfits: -> %d (%d bytes)\n", answer/2880, answer);
        );
    }
    return answer ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Add pointer information about a file into the qfits cache.
  @param    filename    Name of the file to examine.
  @return    index to the file information in the cache, or -1 if failure.

  This function picks a file name, and examines the corresponding FITS file
  to deduce all relevant pointers in the file (byte offsets). These byte
  offsets are later used to speed up header lookups. Example: requesting
  some keyword information in the header of the n-th extension will first
  fseek the file to the header start, then search from this position
  onwards. This means that the input FITS file is only parsed for extension
  positions once.

  What this function does is:

  - Open the file, read the first FITS block (@c FITS_BLOCK_SIZE bytes)
  - Check the file is FITS (must have SIMPLE  = at top)
  - Register start of first header at offset 0.
  - Look for END keyword, register start of first data section
    if NAXIS>0.
  - If the EXTEND=T line was found, continue looking for extensions.
  - For each consecutive extension, register extension header start
    and extension data start.
 */
/*----------------------------------------------------------------------------*/
static int qfits_cache_add(const char * filename)
{
    FILE    *    in ;
    int            off_hdr[QFITS_MAX_EXTS];
    int            off_dat[QFITS_MAX_EXTS];
    char        buf[FITS_BLOCK_SIZE] ;
    char    *    buf_c ;
    int            n_blocks ;
    int            found_it ;
    int            xtend ;
    int            naxis ;
    char    *    read_val ;
    int            last ;
    int            end_of_file ;
    int            data_bytes ;
    int            skip_blocks ;
    struct stat sta ;
    int         seeked ;
    int            i ;

    qfits_cache_cell * qc ;

    /* Initialize cache if not done yet (done only once) */
    if (qfits_cache_init==0) {
        qfits_cache_init++ ;
        qfits_cache_activate();
    }

    /* Stat file to get its size */
    if (stat(filename, &sta)!=0) {
        qdebug(
            printf("qfits: cannot stat file %s\n", filename);
        );
        return -1 ;
    }

    /* Open input file */
    if ((in=fopen(filename, "r"))==NULL) {
        qdebug(
            printf("qfits: cannot open file %s\n", filename);
        );
        return -1 ;
    }

    /* Read first block in */
    if (fread(buf, 1, FITS_BLOCK_SIZE, in)!=FITS_BLOCK_SIZE) {
        qdebug(
            printf("qfits: error reading first block from %s\n", filename);
        );
        fclose(in);
        return -1 ;
    }
    /* Identify FITS magic number */
    if (buf[0]!='S' ||
        buf[1]!='I' ||
        buf[2]!='M' ||
        buf[3]!='P' ||
        buf[4]!='L' ||
        buf[5]!='E' ||
        buf[6]!=' ' ||
        buf[7]!=' ' ||
        buf[8]!='=') {
        qdebug(
            printf("qfits: file %s is not FITS\n", filename);
        );
        fclose(in);
        return -1 ;
    }

    /*
     * Browse through file to identify primary HDU size and see if there
     * might be some extensions. The size of the primary data zone will
     * also be estimated from the gathering of the NAXIS?? values and
     * BITPIX.
     */

    /* Rewind input file, END card might be in first block */
    rewind(in);

    /* Initialize all counters */
    n_blocks = 0 ;
    found_it = 0 ;
    xtend = 0 ;
    naxis = 0 ;
    data_bytes = 1 ;

    /* Start looking for END card */
    while (found_it==0) {
        /* Read one FITS block */
        if (fread(buf, 1, FITS_BLOCK_SIZE, in)!=FITS_BLOCK_SIZE) {
            qdebug(
                printf("qfits: error reading file %s\n", filename);
            );
            fclose(in);
            return -1 ;
        }
        n_blocks ++ ;
        /* Browse through current block */
        buf_c = buf ;
        for (i=0 ; i<FITS_NCARDS ; i++) {

            /* Look for BITPIX keyword */
            if (buf_c[0]=='B' &&
                buf_c[1]=='I' &&
                buf_c[2]=='T' &&
                buf_c[3]=='P' &&
                buf_c[4]=='I' &&
                buf_c[5]=='X' &&
                buf_c[6]==' ') {
                read_val = qfits_getvalue(buf_c);
                data_bytes *= (int)atoi(read_val) / 8 ;
                if (data_bytes<0) data_bytes *= -1 ;
            } else
            /* Look for NAXIS keyword */
            if (buf_c[0]=='N' &&
                buf_c[1]=='A' &&
                buf_c[2]=='X' &&
                buf_c[3]=='I' &&
                buf_c[4]=='S') {

                if (buf_c[5]==' ') {
                    /* NAXIS keyword */
                    read_val = qfits_getvalue(buf_c);
                    naxis = (int)atoi(read_val);
                } else {
                    /* NAXIS?? keyword (axis size) */
                    read_val = qfits_getvalue(buf_c);
                    data_bytes *= (int)atoi(read_val);
                }
            } else
            /* Look for EXTEND keyword */
            if (buf_c[0]=='E' &&
                buf_c[1]=='X' &&
                buf_c[2]=='T' &&
                buf_c[3]=='E' &&
                buf_c[4]=='N' &&
                buf_c[5]=='D' &&
                buf_c[6]==' ') {
                /* The EXTEND keyword is present: might be some extensions */
                read_val = qfits_getvalue(buf_c);
                if (read_val[0]=='T' || read_val[0]=='1') {
                    xtend=1 ;
                }
            } else
            /* Look for END keyword */
            if (buf_c[0] == 'E' &&
                buf_c[1] == 'N' &&
                buf_c[2] == 'D' &&
                buf_c[3] == ' ') {
                found_it = 1 ;
            }
            buf_c += FITS_LINESZ ;
        }
    }

    /*
     * Prepare qfits cache for addition of a new entry
     */
    qfits_cache_last++ ;
    /* Rotate buffer if needed */
    if (qfits_cache_last >= QFITS_CACHESZ) {
        qfits_cache_last = 0 ;
    }
    /* Alias to current pointer in cache for easier reading */
    qc = &(qfits_cache[qfits_cache_last]);

    /* Clean cache cell if needed */
    if (qc->name!=NULL) {
        qfits_free(qc->name) ;
        qc->name = NULL ;
        qfits_free(qc->ohdr);
        qfits_free(qc->data);
        qfits_free(qc->shdr);
        qfits_free(qc->dsiz);
        qfits_cache_entries -- ;
    }

    /* Initialize cache cell */
    qc->exts=0 ;
    qc->name = qfits_strdup(filename);
    qc->inode= sta.st_ino ;

    /* Set first HDU offsets */
    off_hdr[0] = 0 ;
    off_dat[0] = n_blocks ;

    /* Last is the pointer to the last added extension, plus one. */
    last = 1 ;

    if (xtend) {
        /* Look for extensions */
        qdebug(
            printf("qfits: searching for extensions in %s\n", filename);
        );

        /*
         * Register all extension offsets
         */
        end_of_file = 0 ;
        while (end_of_file==0) {
            /*
             * Skip the previous data section if pixels were declared
             */
            if (naxis>0) {
                /* Skip as many blocks as there are declared pixels */
                skip_blocks = data_bytes/FITS_BLOCK_SIZE ;
                if ((data_bytes % FITS_BLOCK_SIZE)!=0) {
                    skip_blocks ++ ;
                }
                seeked = fseek(in, skip_blocks*FITS_BLOCK_SIZE, SEEK_CUR);
                if (seeked<0) {
                    qdebug(
                        printf("qfits: error seeking file %s\n", filename);
                    );
                    qfits_free(qc->name);
                    fclose(in);
                    return -1 ;
                }
                /* Increase counter of current seen blocks. */
                n_blocks += skip_blocks ;
            }

            /* Look for extension start */
            found_it=0 ;
            while ((found_it==0) && (end_of_file==0)) {
                if (fread(buf,1,FITS_BLOCK_SIZE,in)!=FITS_BLOCK_SIZE) {
                    /* Reached end of file */
                    end_of_file=1 ;
                    break ;
                }
                n_blocks ++ ;
                /* Search for XTENSION at block top */
                if (buf[0]=='X' &&
                    buf[1]=='T' &&
                    buf[2]=='E' &&
                    buf[3]=='N' &&
                    buf[4]=='S' &&
                    buf[5]=='I' &&
                    buf[6]=='O' &&
                    buf[7]=='N' &&
                    buf[8]=='=') {
                    /* Got an extension */
                    found_it=1 ;
                    off_hdr[last] = n_blocks-1 ;
                }
            }
            if (end_of_file) break ;

            /*
             * Look for extension END
             * Rewind one block backwards, END might be in same section as
             * XTENSION start.
             */
            if (fseek(in, -FITS_BLOCK_SIZE, SEEK_CUR)==-1) {
                qdebug(
                    printf("qfits: error fseeking file backwards\n");
                ) ;
                qfits_free(qc->name);
                fclose(in);
                return -1 ;
            }
            n_blocks -- ;
            found_it=0 ;
            data_bytes = 1 ;
            naxis = 0 ;
            while ((found_it==0) && (end_of_file==0)) {
                if (fread(buf,1,FITS_BLOCK_SIZE,in)!=FITS_BLOCK_SIZE) {
                    qdebug(
                    printf("qfits: XTENSION without END in %s\n", filename);
                    );
                    end_of_file=1;
                    break ;
                }
                n_blocks++ ;

                /* Browse current block */
                buf_c = buf ;
                for (i=0 ; i<FITS_NCARDS ; i++) {
                    /* Look for BITPIX keyword */
                    if (buf_c[0]=='B' &&
                        buf_c[1]=='I' &&
                        buf_c[2]=='T' &&
                        buf_c[3]=='P' &&
                        buf_c[4]=='I' &&
                        buf_c[5]=='X' &&
                        buf_c[6]==' ') {
                        read_val = qfits_getvalue(buf_c);
                        data_bytes *= (int)atoi(read_val) / 8 ;
                        if (data_bytes<0) data_bytes *= -1 ;
                    } else
                    /* Look for NAXIS keyword */
                    if (buf_c[0]=='N' &&
                        buf_c[1]=='A' &&
                        buf_c[2]=='X' &&
                        buf_c[3]=='I' &&
                        buf_c[4]=='S') {

                        if (buf_c[5]==' ') {
                            /* NAXIS keyword */
                            read_val = qfits_getvalue(buf_c);
                            naxis = (int)atoi(read_val);
                        } else {
                            /* NAXIS?? keyword (axis size) */
                            read_val = qfits_getvalue(buf_c);
                            data_bytes *= (int)atoi(read_val);
                        }
                    } else
                    /* Look for END keyword */
                    if (buf_c[0]=='E' &&
                        buf_c[1]=='N' &&
                        buf_c[2]=='D' &&
                        buf_c[3]==' ') {
                        /* Got the END card */
                        found_it=1 ;
                        /* Update registered extension list */
                        off_dat[last] = n_blocks ;
                        last ++ ;
                        qc->exts ++ ;
                        break ;
                    }
                    buf_c+=FITS_LINESZ ;
                }
            }
        }
    }

    /* Close file */
    fclose(in);

    /* Check last */
    if (last >= QFITS_MAX_EXTS) {
        return -1 ;
    }

    /* Allocate buffers in cache */
    qc->ohdr = qfits_malloc(last * sizeof(int));
    qc->data = qfits_malloc(last * sizeof(int));
    qc->shdr = qfits_malloc(last * sizeof(int));
    qc->dsiz = qfits_malloc(last * sizeof(int));
    /* Store retrieved pointers in the cache */
    for (i=0 ; i<last ; i++) {
        /* Offsets to start */
        qc->ohdr[i] = off_hdr[i];
        qc->data[i] = off_dat[i];

        /* Sizes */
        qc->shdr[i] = off_dat[i] - off_hdr[i]  ;
        if (i==last-1) {
            qc->dsiz[i] = (sta.st_size/FITS_BLOCK_SIZE) - off_dat[i] ;
        } else {
            qc->dsiz[i] = off_hdr[i+1] - off_dat[i] ;
        }
    }
    qc->fsize = sta.st_size / FITS_BLOCK_SIZE ;
    /* Add last modification date */
    qc->mtime = sta.st_mtime ;
    qc->filesize  = sta.st_size ;
    qc->ctime = sta.st_ctime ;
    qfits_cache_entries ++ ;

    qdebug(
        qfits_cache_dump();
    );
    /* Return index of the added file in the cache */
    return qfits_cache_last ;
}

static void qfits_cache_dump(void)
{
    int i, j ;

    printf("qfits: dumping cache...\n");

    printf("cache contains %d entries\n", qfits_cache_entries);
    for (i=0 ; i<QFITS_CACHESZ ; i++) {
        if (qfits_cache[i].name!=NULL) {
            printf("qfits: -----> entry: %d\n", i);
            printf("qfits: name  %s\n", qfits_cache[i].name);
            printf("qfits: exts  %d\n", qfits_cache[i].exts);
            printf("qfits: size  %d\n", qfits_cache[i].fsize);
            printf("qfits: ohdr  %d\n"
                   "qfits: shdr  %d\n"
                   "qfits: data  %d\n"
                   "qfits: dsiz  %d\n",
                   qfits_cache[i].ohdr[0],
                   qfits_cache[i].shdr[0],
                   qfits_cache[i].data[0],
                   qfits_cache[i].dsiz[0]);
            if (qfits_cache[i].exts>0) {
                for (j=1 ; j<=qfits_cache[i].exts ; j++) {
                    printf("qfits: %s [%d]\n", qfits_cache[i].name, j);
                    printf("qfits: ohdr  %d\n"
                           "qfits: shdr  %d\n"
                           "qfits: data  %d\n"
                           "qfits: dsiz  %d\n",
                           qfits_cache[i].ohdr[j],
                           qfits_cache[i].shdr[j],
                           qfits_cache[i].data[j],
                           qfits_cache[i].dsiz[j]);
                }
            }
        }
    }
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    initialize cache buffer with minimum size
 */
/*----------------------------------------------------------------------------*/
static void qfits_cache_activate(void)
{
    int i ;
    qdebug(
        printf("qfits: activating cache...\n");
    );
    /* Set all slots to NULL */
    for (i=0 ; i<QFITS_CACHESZ ; i++) {
        qfits_cache[i].name = NULL ;
    }
    /* Register purge function with atexit */
    atexit(qfits_cache_purge);
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find out if a file is in the cache already
  @param    filename    file name
  @return   int 1 if in the cache, 0 if not
 */
/*----------------------------------------------------------------------------*/
static int qfits_is_cached(const char * filename)
{
    int            i, n ;
    struct stat sta ;

    /* Stat input file */
    if (stat(filename, &sta)!=0) {
        return -1 ;
    }
    n=0 ;
    /* Loop over all cache entries */
    for (i=0 ; i<QFITS_CACHESZ ; i++) {
        /* If entry is valid (name is not NULL) */
        if (qfits_cache[i].name!=NULL) {
            /* One more entry found */
            n++ ;
            /* If inode is the same */
            if ((qfits_cache[i].inode == sta.st_ino) &&
                (qfits_cache[i].mtime == sta.st_mtime) &&
                (qfits_cache[i].filesize  == sta.st_size) &&
                (qfits_cache[i].ctime == sta.st_ctime)) {
                /* This is the requested file */
                return i ;
            }
        }
        /* Early exit: all entries have been browsed */
        if (n>=qfits_cache_entries) {
            return -1 ;
        }
    }
    return -1 ;
}

#endif
#ifndef QFITS_TOOLS_H
#define QFITS_TOOLS_H

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/* Unknown type for FITS value */
#define QFITS_UNKNOWN       0
/* Boolean type for FITS value */
#define QFITS_BOOLEAN       1
/* Int type for FITS value */
#define    QFITS_INT        2
/* Float type for FITS value */
#define QFITS_FLOAT         3
/* Complex type for FITS value */
#define QFITS_COMPLEX       4
/* String type for FITS value */
#define QFITS_STRING        5


/*-----------------------------------------------------------------------------
                            Global variables
 -----------------------------------------------------------------------------*/

/*
 * The following global variables are only used for regular expression
 * matching of integers and floats. These definitions are private to
 * this module.
 */
/** A regular expression matching a floating-point number */
static char regex_float[] =
    "^[+-]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eEdD][+-]?[0-9]+)?$";

/** A regular expression matching an integer */
static char regex_int[] = "^[+-]?[0-9]+$";

/** A regular expression matching a complex number (int or float) */
static char regex_cmp[] =
"^[+-]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eEdD][+-]?[0-9]+)?[ ]+[+-]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eEdD][+-]?[0-9]+)?$";

/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_tools Simple FITS access routines
 *
 *  This module offers a number of very basic low-level FITS access
 *  routines.
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                            Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Retrieve the value of a key in a FITS header
  @param    filename    Name of the FITS file to browse
  @param    keyword     Name of the keyword to find
  @return   pointer to statically allocated character string

  Provide the name of a FITS file and a keyword to look for. The input
  file is memory-mapped and the first keyword matching the requested one is
  located. The value corresponding to this keyword is copied to a
  statically allocated area, so do not modify it or free it.

  The input keyword is first converted to upper case and expanded to
  the HIERARCH scheme if given in the shortFITS notation.

  This function is pretty fast due to the mmapping. Due to buffering
  on most Unixes, it is possible to call many times this function in a
  row on the same file and do not suffer too much from performance
  problems. If the file contents are already in the cache, the file
  will not be re-opened every time.

  It is possible, though, to modify this function to perform several
  searches in a row. See the source code.

  Returns NULL in case the requested keyword cannot be found.
 */
/*----------------------------------------------------------------------------*/
char * qfits_query_hdr(const char * filename, const char * keyword)
{
    return qfits_query_ext(filename, keyword, 0);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Retrieve the value of a keyin a FITS extension header.
  @param    filename    name of the FITS file to browse.
  @param    keyword     name of the FITS key to look for.
  @param    xtnum       xtension number
  @return   pointer to statically allocated character string

  Same as qfits_query_hdr but for extensions. xtnum starts from 1 to
  the number of extensions. If xtnum is zero, this function is
  strictly identical to qfits_query_hdr().
 */
/*----------------------------------------------------------------------------*/
char * qfits_query_ext(const char * filename, const char * keyword, int xtnum)
{
    char    *   exp_key ;
    char    *   where ;
    char    *   start ;
    char    *   value ;
    char        test1, test2 ;
    int         i ;
    int         len ;
    int         different ;
    int         seg_start ;
    int         seg_size ;
    long        bufcount ;
    size_t      size ;

    /* Bulletproof entries */
    if (filename==NULL || keyword==NULL || xtnum<0) return NULL ;

    /* Expand keyword */
    exp_key = qfits_expand_keyword(keyword);

    /*
     * Find out offsets to the required extension
     * Record the xtension start and stop offsets
     */
    if (qfits_get_hdrinfo(filename, xtnum, &seg_start, &seg_size)==-1) {
        return NULL ;
    }

    /*
     * Get a hand on requested buffer
     */

    start = qfits_falloc((char *)filename, seg_start, &size);
    if (start==NULL) return NULL ;

    /*
     * Look for keyword in header
     */

    bufcount=0 ;
    where = start ;
    len = (int)strlen(exp_key);
    while (1) {
        different=0 ;
        for (i=0 ; i<len ; i++) {
            if (where[i]!=exp_key[i]) {
                different++ ;
                break ;
            }
        }
        if (!different) {
            /* Get 2 chars after keyword */
            test1=where[len];
            test2=where[len+1];
            /* If first subsequent character is the equal sign, bingo. */
            if (test1=='=') break ;
            /* If subsequent char is equal sign, bingo */
            if (test1==' ' && (test2=='=' || test2==' '))
                break ;
        }
        /* Watch out for header end */
        if ((where[0]=='E') &&
            (where[1]=='N') &&
            (where[2]=='D') &&
            (where[3]==' ')) {
            /* Detected header end */
            qfits_fdealloc(start, seg_start, size) ;
            return NULL ;
        }
        /* Forward one line */
        where += 80 ;
        bufcount += 80 ;
        if (bufcount>seg_size) {
            /* File is damaged or not FITS: bailout */
            qfits_fdealloc(start, seg_start, size) ;
            return NULL ;
        }
    }

    /* Found the keyword, now get its value */
    value = qfits_getvalue(where);
    qfits_fdealloc(start, seg_start, size) ;
    return value;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Counts the number of extensions in a FITS file
  @param    filename    Name of the FITS file to browse.
  @return   int
  Counts how many extensions are in the file. Returns 0 if no
  extension is found, and -1 if an error occurred.
 */
/*----------------------------------------------------------------------------*/
int qfits_query_n_ext(const char * filename)
{
    return qfits_query(filename, QFITS_QUERY_N_EXT);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Counts the number of planes in a FITS extension.
  @param    filename    Name of the FITS file to browse.
  @param    extnum        Extensin number
  @return   int
  Counts how many planes are in the extension. Returns 0 if no plane is found,
  and -1 if an error occurred.
 */
/*----------------------------------------------------------------------------*/
int qfits_query_nplanes(const char * filename, int extnum)
{
    char    *    sval ;
    int            next ;
    int            naxes ;
    int            nplanes ;

    /* Check file existence */
    if (filename == NULL) return -1 ;
    /* Check validity of extnum */
    next = qfits_query_n_ext(filename) ;
    if (extnum>next) {
        qfits_error("invalid extension specified") ;
        return -1 ;
    }

    /* Find the number of axes  */
    naxes = 0 ;
    if ((sval = qfits_query_ext(filename, "NAXIS", extnum)) == NULL) {
        qfits_error("missing key in header: NAXIS");
        return -1 ;
    }
    naxes = atoi(sval) ;

    /* Check validity of naxes */
    if ((naxes < 2) || (naxes > 3)) return -1 ;

    /* Two dimensions cube */
    if (naxes == 2) nplanes = 1 ;
    else {
        /* For 3D cubes, get the third dimension size   */
        if ((sval = qfits_query_ext(filename, "NAXIS3", extnum))==NULL) {
            qfits_error("missing key in header: NAXIS3");
            return -1 ;
        }
        nplanes = atoi(sval);
        if (nplanes < 1) nplanes = 0 ;
    }
    return nplanes ;
}

#define PRETTY_STRING_STATICBUFS    8
/*----------------------------------------------------------------------------*/
/**
  @brief    Clean out a FITS string value.
  @param    s pointer to allocated FITS value string.
  @return   pointer to statically allocated character string

  From a string FITS value like 'marvin o''hara', remove head and tail
  quotes, replace double '' with simple ', trim blanks on each side,
  and return the result in a statically allocated area.

  Examples:

  - ['o''hara'] becomes [o'hara]
  - ['  H    '] becomes [H]
  - ['1.0    '] becomes [1.0]

 */
/*----------------------------------------------------------------------------*/
char * qfits_pretty_string(const char * s)
{
    static char     pretty_buf[PRETTY_STRING_STATICBUFS][81] ;
    static int      flip=0 ;
    char        *   pretty ;
    int             i,j ;

    /* bulletproof */
    if (s==NULL) return NULL ;

    /* Switch between static buffers */
    pretty = pretty_buf[flip];
    flip++ ;
    if (flip==PRETTY_STRING_STATICBUFS)
        flip=0 ;

    pretty[0] = (char)0 ;
    if (s[0]!='\'') return (char *)s ;

    /* skip first quote */
    i=1 ;
    j=0 ;
    /* trim left-side blanks */
    while (s[i]==' ') {
        if (i==(int)strlen(s)) break ;
        i++ ;
    }
    if (i>=(int)(strlen(s)-1)) return pretty ;
    /* copy string, changing double quotes to single ones */
    while (i<(int)strlen(s)) {
        if (s[i]=='\'') {
            i++ ;
        }
        pretty[j]=s[i];
        i++ ;
        j++ ;
    }
    /* NULL-terminate the pretty string */
    pretty[j+1]=(char)0;
    /* trim right-side blanks */
    j = (int)strlen(pretty)-1;
    while (pretty[j]==' ') j-- ;
    pretty[j+1]=(char)0;
    return pretty;
}
#undef PRETTY_STRING_STATICBUFS

/*----------------------------------------------------------------------------*/
/**
  @brief    Identify if a FITS value is boolean
  @param    s FITS value as a string
  @return   int 0 or 1

  Identifies if a FITS value is boolean.
 */
/*----------------------------------------------------------------------------*/
int qfits_is_boolean(const char * s)
{
    if (s==NULL) return 0 ;
    if (s[0]==0) return 0 ;
    if ((int)strlen(s)>1) return 0 ;
    if (s[0]=='T' || s[0]=='F') return 1 ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Identify if a FITS value is an int.
  @param    s FITS value as a string
  @return   int 0 or 1

  Identifies if a FITS value is an integer.
 */
/*----------------------------------------------------------------------------*/
int qfits_is_int(const char * s)
{
    regex_t re_int ;
    int     status ;

    if (s==NULL) return 0 ;
    if (s[0]==0) return 0 ;
    if (regcomp(&re_int, &regex_int[0], REG_EXTENDED|REG_NOSUB)!=0) {
        qfits_error("internal error: compiling int rule");
        exit(-1);
    }
    status = regexec(&re_int, s, 0, NULL, 0) ;
    regfree(&re_int) ;
    return (status) ? 0 : 1 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Identify if a FITS value is float.
  @param    s FITS value as a string
  @return   int 0 or 1

  Identifies if a FITS value is float.
 */
/*----------------------------------------------------------------------------*/
int qfits_is_float(const char * s)
{
    regex_t re_float;
    int     status ;

    if (s==NULL) return 0 ;
    if (s[0]==0) return 0 ;
    if (regcomp(&re_float, &regex_float[0], REG_EXTENDED|REG_NOSUB)!=0) {
        qfits_error("internal error: compiling float rule");
        exit(-1);
    }
    status = regexec(&re_float, s, 0, NULL, 0) ;
    regfree(&re_float) ;
    return (status) ? 0 : 1 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Identify if a FITS value is complex.
  @param    s FITS value as a string
  @return   int 0 or 1

  Identifies if a FITS value is complex.
 */
/*----------------------------------------------------------------------------*/
int qfits_is_complex(const char * s)
{
    regex_t re_cmp ;
    int     status ;

    if (s==NULL) return 0 ;
    if (s[0]==0) return 0 ;
    if (regcomp(&re_cmp, &regex_cmp[0], REG_EXTENDED|REG_NOSUB)!=0) {
        qfits_error("internal error: compiling complex rule");
        exit(-1);
    }
    status = regexec(&re_cmp, s, 0, NULL, 0) ;
    regfree(&re_cmp) ;
    return (status) ? 0 : 1 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Identify if a FITS value is string.
  @param    s FITS value as a string
  @return   int 0 or 1

  Identifies if a FITS value is a string.
 */
/*----------------------------------------------------------------------------*/
int qfits_is_string(const char * s)
{
    if (s==NULL) return 0 ;
    if (s[0]=='\'') return 1 ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Identify the type of a FITS value given as a string.
  @param    s FITS value as a string
  @return   integer naming the FITS type

  Returns the following value:

  - QFITS_UNKNOWN (0) for an unknown type.
  - QFITS_BOOLEAN (1) for a boolean type.
  - QFITS_INT (2) for an integer type.
  - QFITS_FLOAT (3) for a floating-point type.
  - QFITS_COMPLEX (4) for a complex number.
  - QFITS_STRING (5) for a FITS string.
 */
/*----------------------------------------------------------------------------*/
int qfits_get_type(const char * s)
{
    if (s==NULL) return QFITS_UNKNOWN ;
    if (qfits_is_boolean(s)) return QFITS_BOOLEAN ;
    if (qfits_is_int(s)) return QFITS_INT ;
    if (qfits_is_float(s)) return QFITS_FLOAT ;
    if (qfits_is_complex(s)) return QFITS_COMPLEX ;
    return QFITS_STRING ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Query a card in a FITS (main) header by a given key
  @param    filename    Name of the FITS file to check.
  @param    keyword     Where to read a card in the header.
  @return   Allocated string containing the card or NULL
 */
/*----------------------------------------------------------------------------*/
char * qfits_query_card(
        const char  *   filename,
        const char  *   keyword)
{
    char    *   exp_key ;
    int         fd ;
    char    *   buf ;
    char    *   buf2 ;
    char    *   where ;
    int         hs ;
    char    *   card ;

    /* Bulletproof entries */
    if (filename==NULL || keyword==NULL) return NULL ;

    /* Expand keyword */
    exp_key = qfits_expand_keyword(keyword) ;

    /* Memory-map the FITS header of the input file  */
    qfits_get_hdrinfo(filename, 0, NULL, &hs) ;
    if (hs < 1) {
        qfits_error("error getting FITS header size for %s", filename);
        return NULL ;
    }
    fd = open(filename, O_RDWR) ;
    if (fd == -1) return NULL ;
    buf = (char*)mmap(0,
                      hs,
                      PROT_READ | PROT_WRITE,
                      MAP_SHARED,
                      fd,
                      0) ;
    if (buf == (char*)-1) {
        perror("mmap") ;
        close(fd) ;
        return NULL ;
    }

    /* Apply search for the input keyword */
    buf2 = qfits_malloc(hs+1) ;
    memcpy(buf2, buf, hs) ;
    buf2[hs] = (char)0 ;
    where = buf2 ;
    do {
        where = strstr(where, exp_key);
        if (where == NULL) {
            close(fd);
            munmap(buf,hs);
            qfits_free(buf2) ;
            return NULL ;
        }
        if ((where-buf2)%80) where++ ;
    } while ((where-buf2)%80) ;

    where = buf + (int)(where - buf2) ;

    /* Create the card */
    card = qfits_malloc(81*sizeof(char)) ;
    strncpy(card, where, 80) ;
    card[80] = (char)0 ;

    /* Free and return */
    close(fd) ;
    munmap(buf, hs) ;
    qfits_free(buf2) ;
    return card ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Replace a card in a FITS (main) header by a given card
  @param    filename    Name of the FITS file to modify.
  @param    keyword     Where to substitute a card in the header.
  @param    substitute  What to replace the line with.
  @return   int 0 if Ok, -1 otherwise

  Replaces a whole card (80 chars) in a FITS header by a given FITS
  line (80 chars). The replacing line is assumed correctly formatted
  and containing at least 80 characters. The file is modified: it must
  be accessible in read/write mode.

  The input keyword is first converted to upper case and expanded to
  the HIERARCH scheme if given in the shortFITS notation.

  Returns 0 if everything worked Ok, -1 otherwise.
 */
/*----------------------------------------------------------------------------*/
int qfits_replace_card(
        const char  *   filename,
        const char  *   keyword,
        const char  *   substitute)
{
    char    *   exp_key ;
    int         fd ;
    char    *   buf ;
    char    *   buf2 ;
    char    *   where ;
    int         hs ;


    /* Bulletproof entries */
    if (filename==NULL || keyword==NULL || substitute==NULL) return -1 ;

    /* Expand keyword */
    exp_key = qfits_expand_keyword(keyword);
    /*
     * Memory-map the FITS header of the input file
     */

    qfits_get_hdrinfo(filename, 0, NULL, &hs) ;
    if (hs < 1) {
        qfits_error("error getting FITS header size for %s", filename);
        return -1 ;
    }
    fd = open(filename, O_RDWR) ;
    if (fd == -1) {
        return -1 ;
    }
    buf = (char*)mmap(0,
                      hs,
                      PROT_READ | PROT_WRITE,
                      MAP_SHARED,
                      fd,
                      0) ;
    if (buf == (char*)-1) {
        perror("mmap") ;
        close(fd) ;
        return -1 ;
    }

    /* Apply search and replace for the input keyword lists */
    buf2 = qfits_malloc(hs+1) ;
    memcpy(buf2, buf, hs) ;
    buf2[hs] = (char)0 ;
    where = buf2 ;
    do {
        where = strstr(where, exp_key);
        if (where == NULL) {
            close(fd);
            munmap(buf,hs);
            qfits_free(buf2) ;
            return -1 ;
        }
        if ((where-buf2)%80) where++ ;
    } while ((where-buf2)%80) ;

    where = buf + (int)(where - buf2) ;

    /* Replace current placeholder by blanks */
    memset(where, ' ', 80);
    /* Copy substitute into placeholder */
    memcpy(where, substitute, strlen(substitute));

    close(fd) ;
    munmap(buf, hs) ;
    qfits_free(buf2) ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Return the current QFITS version
  @return   the QFITS version
 */
/*----------------------------------------------------------------------------*/
const char * qfits_version(void)
{
    return (const char *)PACKAGE_VERSION ;
}
#endif
#ifndef QFITS_CARD_H
#define QFITS_CARD_H

/*-----------------------------------------------------------------------------
                              Static functions
 -----------------------------------------------------------------------------*/

static char * expkey_strupc(const char *) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_card   Card handling functions
 *
 * This module contains various routines to help parsing a single FITS
 * card into its components: key, value, comment.
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Write out a card to a string on 80 chars.
  @param    line    Allocated output character buffer.
  @param    key     Key to write.
  @param    val     Value to write.
  @param    com     Comment to write.
  @return   void

  Write out a key, value and comment into an allocated character buffer.
  The buffer must be at least 80 chars to receive the information.
  Formatting is done according to FITS standard.
 */
/*----------------------------------------------------------------------------*/
void qfits_card_build(
        char        *   line,
        const char  *   key,
        const char  *   val,
        const char  *   com)
{
    int     len ;
    int     hierarch = 0 ;
    char    cval[81];
    char    cval2[81];
    char    cval_q[81];
    char    ccom[81];
    char    safe_line[512];
    int     i, j ;

    if (line==NULL || key==NULL) return ;

    /* Set the line with zeroes */
    memset(line, ' ', 80);
    if (key==NULL) return ;

    /* END keyword*/
    if (!strcmp(key, "END")) {
        /* Write key and return */
        sprintf(line, "END") ;
        return ;
    }
    /* HISTORY, COMMENT and blank keywords */
    if (!strcmp(key, "HISTORY") ||
        !strcmp(key, "COMMENT") ||
        !strncmp(key, "        ", 8)) {
        /* Write key */
        sprintf(line, "%s ", key);
        if (val==NULL) return ;

        /* There is a value to write, copy it correctly */
        len = strlen(val);
        /* 72 is 80 (FITS line size) - 8 (sizeof COMMENT or HISTORY) */
        if (len>72) len=72 ;
        strncpy(line+8, val, len);
        return ;
    }

    /* Check for NULL values */
    if (val==NULL) cval[0]=(char)0;
    else if (strlen(val)<1) cval[0]=(char)0;
    else strcpy(cval, val);

    /* Check for NULL comments */
    if (com==NULL) strcpy(ccom, "no comment");
    else strcpy(ccom, com);

    /* Set hierarch flag */
    if (!strncmp(key, "HIERARCH", 8)) hierarch ++ ;

    /* Boolean, int, float or complex */
    if (qfits_is_int(cval) ||
            qfits_is_float(cval) ||
            qfits_is_boolean(cval) ||
            qfits_is_complex(cval)) {
        if (hierarch) sprintf(safe_line, "%-29s= %s / %s", key, cval, ccom);
        else sprintf(safe_line, "%-8.8s= %20s / %-48s", key, cval, ccom);
        strncpy(line, safe_line, 80);
        line[80]=(char)0;
        return ;
    }

    /* Blank or NULL values */
    if (cval[0]==0) {
        if (hierarch) {
            sprintf(safe_line, "%-29s=                    / %s", key, ccom);
        } else {
        sprintf(safe_line, "%-8.8s=                      / %-48s", key, ccom);
        }
        strncpy(line, safe_line, 80);
        line[80]=(char)0;
        return ;
    }

    /* Can only be a string - Make simple quotes ['] as double [''] */
    memset(cval_q, 0, 81);
    strcpy(cval2, qfits_pretty_string(cval));
    j=0 ;
    i=0 ;
    while (cval2[i] != (char)0) {
        if (cval2[i]=='\'') {
            cval_q[j]='\'';
            j++ ;
            cval_q[j]='\'';
        } else {
            cval_q[j] = cval2[i];
        }
        i++ ;
        j++ ;
    }

    if (hierarch) {
        sprintf(safe_line, "%-29s= '%s' / %s", key, cval_q, ccom);
        if (strlen(key) + strlen(cval_q) + 3 >= 80)
            safe_line[79] = '\'';
    } else {
        sprintf(safe_line, "%-8.8s= '%-8s' / %s", key, cval_q, ccom);
    }
    strncpy(line, safe_line, 80);

    /* Null-terminate in any case */
    line[80]=(char)0;
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find the keyword in a key card (80 chars)
  @param    line allocated 80-char line from a FITS header
  @return    statically allocated char *

  Find out the part of a FITS line corresponding to the keyword.
  Returns NULL in case of error. The returned pointer is statically
  allocated in this function, so do not modify or try to free it.
 */
/*----------------------------------------------------------------------------*/
char * qfits_getkey(const char * line)
{
    static char     key[81];
    int                i ;

    if (line==NULL) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getkey: NULL input line\n");
#endif
        return NULL ;
    }

    /* Special case: blank keyword */
    if (!strncmp(line, "        ", 8)) {
        strcpy(key, "        ");
        return key ;
    }
    /* Sort out special cases: HISTORY, COMMENT, END do not have = in line */
    if (!strncmp(line, "HISTORY ", 8)) {
        strcpy(key, "HISTORY");
        return key ;
    }
    if (!strncmp(line, "COMMENT ", 8)) {
        strcpy(key, "COMMENT");
        return key ;
    }
    if (!strncmp(line, "END ", 4)) {
        strcpy(key, "END");
        return key ;
    }

    memset(key, 0, 81);
    /* General case: look for the first equal sign */
    i=0 ;
    while (line[i]!='=' && i<80) i++ ;
    if (i>=80) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getkey: cannot find equal sign\n");
#endif
        return NULL ;
    }
    i-- ;
    /* Equal sign found, now backtrack on blanks */
    while (line[i]==' ' && i>=0) i-- ;
    if (i<=0) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getkey: error backtracking on blanks\n");
#endif
        return NULL ;
    }
    i++ ;

    /* Copy relevant characters into output buffer */
    strncpy(key, line, i) ;
    /* Null-terminate the string */
    key[i+1] = (char)0;
    return key ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find the value in a key card (80 chars)
  @param    line allocated 80-char line from a FITS header
  @return    statically allocated char *

  Find out the part of a FITS line corresponding to the value.
  Returns NULL in case of error, or if no value can be found. The
  returned pointer is statically allocated in this function, so do not
  modify or try to free it.
 */
/*----------------------------------------------------------------------------*/
char * qfits_getvalue(const char * line)
{
    static char value[81] ;
    int     i ;
    int     from, to ;
    int     inq ;

    if (line==NULL) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getvalue: NULL input line\n");
#endif
        return NULL ;
    }

    /* Special cases */

    /* END has no associated value */
    if (!strncmp(line, "END ", 4)) {
        return NULL ;
    }
    /*
     * HISTORY has for value everything else on the line, stripping
     * blanks before and after. Blank HISTORY is also accepted.
     */
    memset(value, 0, 81);

    if (!strncmp(line, "HISTORY ", 8) || !strncmp(line, "        ", 8)) {
        i=7 ;
        /* Strip blanks from the left side */
        while (line[i]==' ' && i<80) i++ ;
        if (i>=80) return NULL ; /* Blank HISTORY */
        from=i ;

        /* Strip blanks from the right side */
        to=79 ;
        while (line[to]==' ') to-- ;
        /* Copy relevant characters into output buffer */
        strncpy(value, line+from, to-from+1);
        /* Null-terminate the string */
        value[to-from+1] = (char)0;
        return value ;
    } else if (!strncmp(line, "COMMENT ", 8)) {
        /* COMMENT is like HISTORY */
        /* Strip blanks from the left side */
        i=7 ;
        while (line[i]==' ' && i<80) i++ ;
        if (i>=80) return NULL ;
        from=i ;

        /* Strip blanks from the right side */
        to=79 ;
        while (line[to]==' ') to-- ;

        if (to<from) {
#ifdef DEBUG_FITSHEADER
            printf("qfits_getvalue: inconsistent value search in COMMENT\n");
#endif
            return NULL ;
        }
        /* Copy relevant characters into output buffer */
        strncpy(value, line+from, to-from+1);
        /* Null-terminate the string */
        value[to-from+1] = (char)0;
        return value ;
    }
    /* General case - Get past the keyword */
    i=0 ;
    while (line[i]!='=' && i<80) i++ ;
    if (i>80) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getvalue: no equal sign found on line\n");
#endif
        return NULL ;
    }
    i++ ;
    while (line[i]==' ' && i<80) i++ ;
    if (i>80) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getvalue: no value past the equal sign\n");
#endif
        return NULL ;
    }
    from=i;

    /* Now value section: Look for the first slash '/' outside a string */
    inq = 0 ;
    while (i<80) {
        if (line[i]=='\'')
            inq=!inq ;
        if (line[i]=='/')
            if (!inq)
                break ;
        i++;
    }
    i-- ;

    /* Backtrack on blanks */
    while (line[i]==' ' && i>=0) i-- ;
    if (i<0) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getvalue: error backtracking on blanks\n");
#endif
        return NULL ;
    }
    to=i ;

    if (to<from) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getvalue: from>to?\n");
        printf("line=[%s]\n", line);
#endif
        return NULL ;
    }
    /* Copy relevant characters into output buffer */
    strncpy(value, line+from, to-from+1);
    /* Null-terminate the string */
    value[to-from+1] = (char)0;
    return value ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find the comment in a key card (80 chars)
  @param    line allocated 80-char line from a FITS header
  @return    statically allocated char *

  Find out the part of a FITS line corresponding to the comment.
  Returns NULL in case of error, or if no comment can be found. The
  returned pointer is statically allocated in this function, so do not
  modify or try to free it.
 */
/*----------------------------------------------------------------------------*/
char * qfits_getcomment(const char * line)
{
    static char comment[81];
    int    i ;
    int    from, to ;
    int    inq ;

    if (line==NULL) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getcomment: null line in input\n");
#endif
        return NULL ;
    }

    /* Special cases: END, HISTORY, COMMENT and blank have no comment */
    if (!strncmp(line, "END ", 4)) return NULL ;
    if (!strncmp(line, "HISTORY ", 8)) return NULL ;
    if (!strncmp(line, "COMMENT ", 8)) return NULL ;
    if (!strncmp(line, "        ", 8)) return NULL ;

    memset(comment, 0, 81);
    /* Get past the keyword */
    i=0 ;
    while (line[i]!='=' && i<80) i++ ;
    if (i>=80) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getcomment: no equal sign on line\n");
#endif
        return NULL ;
    }
    i++ ;

    /* Get past the value until the slash */
    inq = 0 ;
    while (i<80) {
        if (line[i]=='\'')
            inq = !inq ;
        if (line[i]=='/')
            if (!inq)
                break ;
        i++ ;
    }
    if (i>=80) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getcomment: no slash found on line\n");
#endif
        return NULL ;
    }
    i++ ;
    /* Get past the first blanks */
    while (line[i]==' ') i++ ;
    from=i ;

    /* Now backtrack from the end of the line to the first non-blank char */
    to=79 ;
    while (line[to]==' ') to-- ;

    if (to<from) {
#ifdef DEBUG_FITSHEADER
        printf("qfits_getcomment: from>to?\n");
#endif
        return NULL ;
    }
    /* Copy relevant characters into output buffer */
    strncpy(comment, line+from, to-from+1);
    /* Null-terminate the string */
    comment[to-from+1] = (char)0;
    return comment ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Expand a keyword from shortFITS to HIERARCH notation.
  @param    keyword        Keyword to expand.
  @return    1 pointer to statically allocated string.

  This function expands a given keyword from shortFITS to HIERARCH
  notation, bringing it to uppercase at the same time.

  Examples:

  @verbatim
  det.dit          expands to     HIERARCH ESO DET DIT
  ins.filt1.id     expands to     HIERARCH ESO INS FILT1 ID
  @endverbatim

  If the input keyword is a regular FITS keyword (i.e. it contains
  not dots '.') the result is identical to the input.
 */
/*----------------------------------------------------------------------------*/
char * qfits_expand_keyword(const char * keyword)
{
    static char expanded[81];
    char        ws[81];
    char    *    token ;

    /* Bulletproof entries */
    if (keyword==NULL) return NULL ;
    /* If regular keyword, copy the uppercased input and return */
    if (strstr(keyword, ".")==NULL) {
        strcpy(expanded, expkey_strupc(keyword));
        return expanded ;
    }
    /* Regular shortFITS keyword */
    sprintf(expanded, "HIERARCH ESO");
    strcpy(ws, expkey_strupc(keyword));
    token = strtok(ws, ".");
    while (token!=NULL) {
        strcat(expanded, " ");
        strcat(expanded, token);
        token = strtok(NULL, ".");
    }
    return expanded ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Uppercase a string
  @param    s   string
  @return   string
 */
/*----------------------------------------------------------------------------*/
static char * expkey_strupc(const char * s)
{
    static char l[1024+1];
    int i ;

    if (s==NULL) return NULL ;
    memset(l, 0, 1024+1);
    i=0 ;
    while (s[i] && i<1024) {
        l[i] = (char)toupper((int)s[i]);
        i++ ;
    }
    l[1024]=(char)0;
    return l ;
}
#endif

#ifndef QFITS_ERROR_H
#define QFITS_ERROR_H

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/* Max number of error handling functions registered */
#define QFITS_ERR_MAXERRDISP        8
/* Max size of an error message */
#define QFITS_ERR_MSGSIZE            1024

/*-----------------------------------------------------------------------------
                               Static declarations
 -----------------------------------------------------------------------------*/

/* Type of a display function only defined for legibility here */
typedef void (*qfits_err_dispfunc)(char *) ;

/* Default display function prints out msg to stderr */
static void qfits_err_display_stderr(char * s)
{ fprintf(stderr, "qfits: %s\n", s); }

/* Static control structure, completely private */
static struct {
    qfits_err_dispfunc     disp[QFITS_ERR_MAXERRDISP] ;
    int                 n ;
    int                    active ;
} qfits_err_control = {{qfits_err_display_stderr}, 1, 0} ;

static int qfits_err_statget(void) ;
static int qfits_err_statset(int) ;
static int qfits_err_register( void (*dispfn)(char*) ) ;
static void qfits_err_main_display(char *) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_error     Messaging functionalities
 *
 *   This module is responsible for error message display. It allows
 *  to re-direct all messages to a given set of functions for display.
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

/* Public warning/error functions */
void qfits_warning(const char *fmt, ...)
{
    char msg[QFITS_ERR_MSGSIZE] ;
    char all[QFITS_ERR_MSGSIZE] ;
    va_list ap ;

    /* Check if display is activated */
    if (qfits_err_control.active==0) {
        return ;
    }
    va_start(ap, fmt) ;
    vsprintf(msg, fmt, ap) ;
    va_end(ap);

    sprintf(all, "*** %s", msg);
    qfits_err_main_display(all);
    return ;
}
void qfits_error(const char *fmt, ...)
{
    char msg[QFITS_ERR_MSGSIZE] ;
    char all[QFITS_ERR_MSGSIZE] ;
    va_list ap ;

    /* Check if display is activated */
    if (qfits_err_control.active==0) {
        return ;
    }
    va_start(ap, fmt) ;
    vsprintf(msg, fmt, ap) ;
    va_end(ap);

    sprintf(all, "error: %s", msg);
    qfits_err_main_display(all);
    return ;
}
/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    used for message display
  @param    msg     message
  @return   nothing
  It calls registered display functions one after another.
 */
/*----------------------------------------------------------------------------*/
static void qfits_err_main_display(char * msg)
{
    int    i ;

    /* Check if there is a message in input */
    if (msg==NULL)
        return ;

    /* Loop on all registered functions and call them */
    for (i=0 ; i<qfits_err_control.n ; i++) {
        if (qfits_err_control.disp[i]) {
            qfits_err_control.disp[i](msg);
        }
    }
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the current status of error display.
  @return    int 1 if error display is active, 0 if not.

  This function returns the current error display status. If it returns 1,
  it means that all calls to qfits_error/qfits_warning will display
  messages using the registered functions, otherwise they do nothing.
 */
/*----------------------------------------------------------------------------*/
static int qfits_err_statget(void)
{
    return qfits_err_control.active ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Set the current status of error display.
  @param    sta        New status to be set.
  @return    int giving the previous display status.

  This function sets the current error display status to the required
  value, and returns the previous value. It is useful to store the
  previous value, in view of restoring it afterwards, e.g. to make a
  function silent on all its errors. Example:

  @code
  int prev_stat = qfits_err_statset(0) ;
  function_call() ;
  qfits_err_statset(prev_stat);
  @endcode
 */
/*----------------------------------------------------------------------------*/
static int qfits_err_statset(int sta)
{
    int prev ;
    prev = qfits_err_control.active ;
    qfits_err_control.active=sta ;
    return prev ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Register a function to display error/warning messages.
  @param    dispfn    Display function (see doc below).
  @return    int 0 if function was registered, -1 if not.

  @code
  void display_function(char * msg);
  @endcode

  They are simple functions that expect a ready-made error message
  and return void. They can do whatever they want with the message
  (log it to a file, send it to a GUI, to the syslog, ...). The
  message is built using a printf-like statement in qfits_error and
  qfits_warning, then passed to all registered display functions.

  A maximum of QFITS_ERR_MAXERRDISP can be registered (see source code).
  If the limit has been reached, this function will signal it by
  returning -1.
 */
/*----------------------------------------------------------------------------*/
static int qfits_err_register(qfits_err_dispfunc dispfn)
{
    if (qfits_err_control.n==QFITS_ERR_MAXERRDISP) {
        /* Cannot register any more function */
        return -1 ;
    }
    qfits_err_control.disp[qfits_err_control.n] = dispfn ;
    qfits_err_control.n ++ ;
    return 0 ;
}
#endif
#ifndef QFITS_FILENAME_H
#define QFITS_FILENAME_H

/*-----------------------------------------------------------------------------
                                  Define
 -----------------------------------------------------------------------------*/

/* Maximum size of a filename buffer */
#define MAXNAMESZ       4096

/*----------------------------------------------------------------------------*/
/**
 * @defgroup   qfits_filename   Get various names (filenames, dir names,...)
 * The following functions are useful to cut out a filename into its components.
 * All functions work with statically allocated memory, i.e. the pointers they
 * return are not permanent but likely to be overwritten at each function call.
 * If you need a returned value later on, you should store it into a local
 * variable.
 *
 * Example:
 *
 * @code
 * char * s ;
 * s = qfits_get_dir_name("/mnt/cdrom/data/image.fits")
 * @endcode
 *
 * s contains now "/mnt/cdrom/data" but will loose these contents at the next
 * function call. To retain its value, you can either do:
 *
 * @code
 * char s[1024];
 * strcpy(s, qfits_get_dir_name("/mnt/cdrom/data/image.fits"));
 * @endcode
 *
 * or:
 *
 * @code
 * char * s;
 * s = strdup(qfits_get_dir_name("/mnt/cdrom/data/image.fits"));
 * ...
 * free(s);
 * @endcode
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Find the directory name in the given string.
  @param    filename    Full path name to scan.
  @return   Pointer to statically allocated string.

  Provide a full path name and you get in return a pointer to a statically
  allocated string containing the name of the directory only, without trailing
  slash. If the input string does not contain a slash (i.e. it is not a full
  path), the returned string is '.', corresponding to the current working
  directory. Since the returned string is statically allocated, do not try to
  free it or modify it.

  This function does not check for the existence of the path or the file.

  Examples:
  @verbatim
    qfits_get_dir_name("/cdrom/data/image.fits") returns "/cdrom/data"
    qfits_get_dir_name("filename.fits") returns "."
  @endverbatim
 */
/*----------------------------------------------------------------------------*/
char * qfits_get_dir_name(const char * filename)
{
    static char path[MAXNAMESZ];
    char *last_slash;

    if (strlen(filename)>MAXNAMESZ) return NULL ;
    strcpy(path, filename);
    /* Find last '/'.  */
    last_slash = path != NULL ? strrchr (path, '/') : NULL;

    if (last_slash == path)
    /* The last slash is the first character in the string.  We have to
    return "/".  */
        ++last_slash;
    else if (last_slash != NULL && last_slash[1] == '\0')
        /* The '/' is the last character, we have to look further.  */
        last_slash = memchr (path, last_slash - path, '/');

    if (last_slash != NULL)
        /* Terminate the path.  */
        last_slash[0] = '\0';
    else
        /* This assignment is ill-designed but the XPG specs require to
        return a string containing "." in any case no directory part is
        found and so a static and constant string is required.  */
        strcpy(path, ".");
    return path;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find out the base name of a file (i.e. without prefix path)
  @param    filename    Full path name to scan.
  @return   Pointer to char within the input string.

  Provide a full path name and you get in return a pointer to a statically
  allocated string containing the name of the file only, without prefixing
  directory names. If the input string does not contain a slash (i.e. it is
  not a full path), the returned string is a copy of the input string.

  This function does not check for the existence of the path or the file.

  Examples:
  @verbatim
    qfits_get_base_name("/cdrom/data/image.fits") returns "image.fits"
    qfits_get_base_name("filename.fits") returns "filename.fits"
  @endverbatim
 */
/*----------------------------------------------------------------------------*/
char * qfits_get_base_name(const char *filename)
{
    char *p ;
    p = strrchr (filename, '/');
    return p ? p + 1 : (char *) filename;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find out the root part of a basename (name without extension).
  @param    filename    File name to scan.
  @return   Pointer to statically allocated string.

  Find out the root part of a file name, i.e. the file name without extension.
  Since in Unix a file name can have several dots, only a number of extensions
  are supported. This includes:

  - .fits and .FITS
  - .tfits and .TFITS
  - .paf and .PAF
  - .ascii and .ASCII
  - .dat and .DAT
  - .txt and .TXT

  This function does not check for the existence of the path or the file.

  Examples:
  @verbatim
    qfits_get_root_name("/cdrom/filename.fits") returns "/cdrom/filename"
    qfits_get_root_name("filename.paf") returns "filename"
    qfits_get_root_name("filename") returns "filename"
    qfits_get_root_name("filename.ext") returns "filename.ext"
  @endverbatim

  Since the returned string is statically allocated in this module, do not try
  to free it or modify its contents.
 */
/*----------------------------------------------------------------------------*/
char * qfits_get_root_name(const char * filename)
{
    static char path[MAXNAMESZ+1];
    char * lastdot ;

    if (strlen(filename)>MAXNAMESZ) return NULL ;
    memset(path, MAXNAMESZ, 0);
    strcpy(path, filename);
    lastdot = strrchr(path, '.');
    if (lastdot == NULL) return path ;
    if ((!strcmp(lastdot, ".fits")) || (!strcmp(lastdot, ".FITS")) ||
        (!strcmp(lastdot, ".paf")) || (!strcmp(lastdot, ".PAF")) ||
        (!strcmp(lastdot, ".dat")) || (!strcmp(lastdot, ".DAT")) ||
        (!strcmp(lastdot, ".txt")) || (!strcmp(lastdot, ".TXT")) ||
        (!strcmp(lastdot, ".tfits")) || (!strcmp(lastdot, ".TFITS")) ||
        (!strcmp(lastdot, ".ascii")) || (!strcmp(lastdot, ".ASCII")))
    {
        lastdot[0] = (char)0;
    }
    return path ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find out the extension of a file name
  @param    filename    File name without path prefix.
  @return   Pointer to char within the input string.

  Find out the extension of a given file name. Notice that the input character
  string must not contain a path prefix (typically, you feed in the output of
  @c qfits_get_base_name).

  Works with all kinds of extensions: returns the part of the string after the
  last dot. If no dot is found in the input string, NULL is returned.

  Examples:
  @verbatim
    qfits_get_ext_name("filename.fits") returns "fits"
    qfits_get_ext_name("hello.c") returns "c"
    qfits_get_ext_name("readme") returns NULL
  @endverbatim
 */
/*----------------------------------------------------------------------------*/
char * qfits_get_ext_name(const char * filename)
{
    char * p;
    p = strrchr(filename, '.');
    return p ? p+1 : NULL ;
}
#endif
#ifndef QFITS_FLOAT_H
#define QFITS_FLOAT_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------------
                                   Macros
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Test a floating-point variable for NaN value.
  @param    n   Number to test (float or double)
  @return   1 if n is NaN, 0 else.

  This macro is needed to support both float and double variables
  as input parameter. It checks on the size of the input variable
  to branch to the float or double version.

  Portability is an issue for this function which is present on
  most Unixes but not all, under various libraries (C lib on BSD,
  Math lib on Linux, sunmath on Solaris, ...). Integrating the
  code for this function makes qfits independent from any math
  library.
 */
/*----------------------------------------------------------------------------*/
#define qfits_isnan(n) ((sizeof(n)==sizeof(float)) ? _qfits_isnanf(n) : \
                        (sizeof(n)==sizeof(double)) ? _qfits_isnand(n) : -1)

/*----------------------------------------------------------------------------*/
/**
  @brief    Test a floating-point variable for Inf value.
  @param    n   Number to test (float or double)
  @return   1 if n is Inf or -Inf, 0 else.

  This macro is needed to support both float and double variables
  as input parameter. It checks on the size of the input variable
  to branch to the float or double version.

  Portability is an issue for this function which is missing on most
  Unixes. Most of the time, another function called finite() is
  offered to perform the opposite task, but it is not consistent
  among platforms and found in various libraries. Integrating the
  code for this function makes qfits independent from any math
  library.
 */
/*----------------------------------------------------------------------------*/
#define qfits_isinf(n) ((sizeof(n)==sizeof(float)) ? _qfits_isinff(n) : \
                        (sizeof(n)==sizeof(double)) ? _qfits_isinfd(n) : -1)

/*-----------------------------------------------------------------------------
                                   New types
 -----------------------------------------------------------------------------*/

#ifndef WORDS_BIGENDIAN
/* Little endian ordering */
typedef union _ieee_double_pattern_ {
    double d ;
    struct {
        unsigned int lsw ;
        unsigned int msw ;
    } p ;
} ieee_double_pattern ;
#else
/* Big endian ordering */
typedef union _ieee_double_pattern_ {
    double d ;
    struct {
        unsigned int msw ;
        unsigned int lsw ;
    } p ;
} ieee_double_pattern ;
#endif

typedef union _ieee_float_pattern_ {
    float f ;
    int   i ;
} ieee_float_pattern ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_float This module implements the qfits_isnan()
 *                          and qfits_isinf() macros
 *
 *  The isnan() and isinf() macros are unfortunately not yet part of
 *  the standard C math library everywhere. They can usually be found
 *  in different places, if they are offered at all, and require the
 *  application to link against the math library. To avoid portability
 *  problems and linking against -lm, this module implements a fast
 *  and portable way of finding out whether a floating-point value
 *  (float or double) is a NaN or an Inf.
 *
 *  Instead of calling isnan() and isinf(), the programmer including
 *  this file should call qfits_isnan() and qfits_isinf().
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

int _qfits_isnanf(float f)
{
    ieee_float_pattern ip ;
    int ix ;

    ip.f = f ;
    ix = ip.i ;
    ix &= 0x7fffffff ;
    ix = 0x7f800000 - ix ;
    return (int)(((unsigned int)(ix))>>31);
}

int _qfits_isinff(float f)
{
    ieee_float_pattern ip ;
    int ix, t ;

    ip.f = f ;
    ix = ip.i ;
    t = ix & 0x7fffffff;
    t ^= 0x7f800000;
    t |= -t;
    return ~(t >> 31) & (ix >> 30);
}

int _qfits_isnand(double d)
{
    ieee_double_pattern id ;
    int hx, lx ;

    id.d = d ;
    lx = id.p.lsw ;
    hx = id.p.msw ;

    hx &= 0x7fffffff;
    hx |= (unsigned int)(lx|(-lx))>>31;
    hx = 0x7ff00000 - hx;
    return (int)(((unsigned int)hx)>>31);
}

int _qfits_isinfd(double d)
{
    ieee_double_pattern id ;
    int hx, lx ;

    id.d = d ;
    lx = id.p.lsw ;
    hx = id.p.msw ;

    lx |= (hx & 0x7fffffff) ^ 0x7ff00000;
    lx |= -lx;
    return ~(lx >> 31) & (hx >> 30);
}
#endif
#ifndef QFITS_HEADER_H
#define QFITS_HEADER_H

/*----------------------------------------------------------------------------*/
/**
  @brief    FITS header object

  This structure represents a FITS header in memory. It is actually no
  more than a thin layer on top of the keytuple object. No field in this
  structure should be directly modifiable by the user, only through
  accessor functions.
 */
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/*
  @brief    keytuple object (internal)

  This structure represents a FITS card (key, val, comment) in memory.
  A FITS header is a list of such structs. Although the struct is here
  made public for convenience, it is not supposed to be directly used.
  High-level FITS routines should do the job just fine, without need
  to know about this struct.
 */
/*----------------------------------------------------------------------------*/
typedef struct _keytuple_ {

    char    *   key ;   /** Key: unique string in a list */
    char    *   val ;   /** Value, always as a string */
    char    *   com ;   /** Comment associated to key */
    char    *   lin ;   /** Initial line in FITS header if applicable */
    int         typ ;   /** Key type */

    /** Implemented as a doubly-linked list */
    struct _keytuple_ * next ;
    struct _keytuple_ * prev ;
} keytuple ;

/*----------------------------------------------------------------------------*/
/**
  @enum        keytype
  @brief    Possible key types

  This enum stores all possible types for a FITS keyword. These determine
  the order of appearance in a header, they are a crucial point for
  DICB (ESO) compatibility. This classification is internal to this
  module.
 */
/*----------------------------------------------------------------------------*/
typedef enum _keytype_ {
    keytype_undef            =0,

    keytype_top                =1,

    /* Mandatory keywords */
    /* All FITS files */
    keytype_bitpix            =2,
    keytype_naxis            =3,

    keytype_naxis1            =11,
    keytype_naxis2            =12,
    keytype_naxis3            =13,
    keytype_naxis4            =14,
    keytype_naxisi            =20,
    /* Random groups only */
    keytype_group            =30,
    /* Extensions */
    keytype_pcount            =31,
    keytype_gcount            =32,
    /* Main header */
    keytype_extend            =33,
    /* Images */
    keytype_bscale            =34,
    keytype_bzero            =35,
    /* Tables */
    keytype_tfields            =36,
    keytype_tbcoli            =40,
    keytype_tformi            =41,

    /* Other primary keywords */
    keytype_primary            =100,

    /* HIERARCH ESO keywords ordered according to DICB */
    keytype_hierarch_dpr    =200,
    keytype_hierarch_obs    =201,
    keytype_hierarch_tpl    =202,
    keytype_hierarch_gen    =203,
    keytype_hierarch_tel    =204,
    keytype_hierarch_ins    =205,
    keytype_hierarch_det    =206,
    keytype_hierarch_log    =207,
    keytype_hierarch_pro    =208,
    /* Other HIERARCH keywords */
    keytype_hierarch        =300,

    /* HISTORY and COMMENT */
    keytype_history            =400,
    keytype_comment            =500,
    /* END */
    keytype_end                =1000
} keytype ;

/*-----------------------------------------------------------------------------
                        Private to this module
 -----------------------------------------------------------------------------*/

static keytuple * keytuple_new(const char *, const char *, const char *,
        const char *);
static void keytuple_del(keytuple *);
static void keytuple_dmp(const keytuple *);
static keytype keytuple_type(const char *);
static int qfits_header_makeline(char *, const keytuple *, int) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_header    FITS header handling
 *
 * This file contains definition and related methods for the FITS header
 * structure. This structure is meant to remain opaque to the user, who
 * only accesses it through the dedicated functions.
 *
 * The 'keytuple' type is strictly internal to this module.
 * It describes FITS cards as tuples (key,value,comment,line), where key
 * is always a non-NULL character string, value and comment are
 * allowed to be NULL. 'line' is a string containing the line as it
 * has been read from the input FITS file (raw). It is set to NULL if
 * the card is modified later. This allows in output two options:
 * either reconstruct the FITS lines by printing key = value / comment
 * in a FITS-compliant way, or output the lines as they were found in
 * input, except for the modified ones.
 *
 * The following functions are associated methods
 * to this data structure:
 *
 * - keytuple_new()      constructor
 * - keytuple_del()      destructor
 * - keytuple_dmp()      dumps a keytuple to stdout
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/
/*----------------------------------------------------------------------------*/
/**
  @brief    FITS header constructor
  @return    1 newly allocated (empty) qfits_header object.

  This is the main constructor for a qfits_header object. It returns
  an allocated linked-list handler with an empty card list.
 */
/*----------------------------------------------------------------------------*/
qfits_header * qfits_header_new(void)
{
    qfits_header    *    h ;
    h = qfits_malloc(sizeof(qfits_header));
    h->first = NULL ;
    h->last  = NULL ;
    h->n = 0 ;

    h->current = NULL ;
    h->current_idx = -1 ;

    return h;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    FITS header default constructor.
  @return    1 newly allocated qfits_header object.

  This is a secondary constructor for a qfits_header object. It returns
  an allocated linked-list handler containing two cards: the first one
  (SIMPLE=T) and the last one (END).

 */
/*----------------------------------------------------------------------------*/
qfits_header * qfits_header_default(void)
{
    qfits_header    *    h ;
    h = qfits_header_new() ;
    qfits_header_append(h, "SIMPLE", "T", "Fits format", NULL);
    qfits_header_append(h, "END", NULL, NULL, NULL);
    return h;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Add a new card to a FITS header
  @param    hdr qfits_header object to modify
  @param    key FITS key
  @param    val FITS value
  @param    com FITS comment
  @param    lin FITS original line if exists
  @return    void

  This function adds a new card into a header, at the one-before-last
  position, i.e. the entry just before the END entry if it is there.
  The key must always be a non-NULL string, all other input parameters
  are allowed to get NULL values.
 */
/*----------------------------------------------------------------------------*/
void qfits_header_add(
        qfits_header    *   hdr,
        const char      *   key,
        const char      *   val,
        const char      *   com,
        const char      *   lin)
{
    keytuple    *    k ;
    keytuple    *    kbf ;
    keytuple    *    first ;
    keytuple    *    last ;

    if (hdr==NULL || key==NULL) return ;
    if (hdr->n<2) return ;

    first = (keytuple*)hdr->first ;
    last  = (keytuple*)hdr->last ;

    if (((keytype)first->typ != keytype_top) ||
        ((keytype)last->typ != keytype_end)) return ;

    /* Create new key tuple */
    k = keytuple_new(key, val, com, lin);

    /* Find the last keytuple with same key type */
    kbf = first ;
    while (kbf!=NULL) {
        if ((k->typ>=kbf->typ) && (k->typ<kbf->next->typ)) break ;
        kbf = kbf->next ;
    }
    if (kbf==NULL) kbf = last->prev ;

    /* Hook it into list */
    k->next = kbf->next ;
    (kbf->next)->prev = k ;
    kbf->next = k ;
    k->prev = kbf ;

    hdr->n ++ ;
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    add a new card to a FITS header
  @param    hdr     qfits_header object to modify
  @param    after    Key to specify insertion place
  @param    key     FITS key
  @param    val     FITS value
  @param    com     FITS comment
  @param    lin     FITS original line if exists
  @return    void

  Adds a new card to a FITS header, after the specified key. Nothing
  happens if the specified key is not found in the header. All fields
  can be NULL, except after and key.
 */
/*----------------------------------------------------------------------------*/
void qfits_header_add_after(
        qfits_header    *   hdr,
        const char      *   after,
        const char      *   key,
        const char      *   val,
        const char      *   com,
        const char      *   lin)
{
    keytuple    *   kreq;
    keytuple    *   k;
    char        *   exp_after ;

    if (hdr==NULL || after==NULL || key==NULL) return ;

    exp_after = qfits_expand_keyword(after);
    /* Locate where the entry is requested */
    kreq = (keytuple*)(hdr->first) ;
    while (kreq!=NULL) {
        if (!strcmp(kreq->key, exp_after)) break ;
        kreq = kreq->next ;
    }
    if (kreq==NULL) return ;
    k = keytuple_new(key, val, com, lin);

    k->next = kreq->next ;
    kreq->next->prev = k ;
    kreq->next = k ;
    k->prev = kreq ;
    hdr->n ++ ;
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Append a new card to a FITS header.
  @param    hdr qfits_header object to modify
  @param    key FITS key
  @param    val FITS value
  @param    com FITS comment
  @param    lin FITS original line if exists
  @return    void

  Adds a new card in a FITS header as the last one. All fields can be
  NULL except key.
 */
/*----------------------------------------------------------------------------*/
void qfits_header_append(
        qfits_header    *   hdr,
        const char      *   key,
        const char      *   val,
        const char      *   com,
        const char      *   lin)
{
    keytuple    *    k;
    keytuple    *    last ;

    if (hdr==NULL || key==NULL) return ;

    k = keytuple_new(key, val, com, lin);
    if (hdr->n==0) {
        hdr->first = hdr->last = k ;
        hdr->n = 1 ;
        return ;
    }
    last  = (keytuple*)hdr->last ;
    last->next = k ;
    k->prev = last ;
    hdr->last = k ;
    hdr->n++ ;
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Delete a card in a FITS header.
  @param    hdr qfits_header to modify
  @param    key specifies which card to remove
  @return    void

  Removes a card from a FITS header. The first found card that matches
  the key is removed.
 */
/*----------------------------------------------------------------------------*/
void qfits_header_del(qfits_header * hdr, const char * key)
{
    keytuple    *   k ;
    char        *   xkey ;

    if (hdr==NULL || key==NULL) return ;

    xkey = qfits_expand_keyword(key);
    k = (keytuple*)hdr->first ;
    while (k!=NULL) {
        if (!strcmp(k->key, xkey)) break ;
        k = k->next ;
    }
    if (k==NULL)
        return ;
    if(k == hdr->first) {
        hdr->first = k->next ;
    } else {
        k->prev->next = k->next ;
        k->next->prev = k->prev ;
    }
    keytuple_del(k);
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Modifies a FITS card.
  @param    hdr qfits_header to modify
  @param    key FITS key
  @param    val FITS value
  @param    com FITS comment
  @return    void

  Finds the first card in the header matching 'key', and replaces its
  value and comment fields by the provided values. The initial FITS
  line is set to NULL in the card.
 */
/*----------------------------------------------------------------------------*/
void qfits_header_mod(
        qfits_header    *   hdr,
        const char      *   key,
        const char      *   val,
        const char      *   com)
{
    keytuple    *   k ;
    char        *   xkey ;

    if (hdr==NULL || key==NULL) return ;

    xkey = qfits_expand_keyword(key);
    k = (keytuple*)hdr->first ;
    while (k!=NULL) {
        if (!strcmp(k->key, xkey)) break ;
        k=k->next ;
    }
    if (k==NULL) return ;

    if (k->val) qfits_free(k->val);
    if (k->com) qfits_free(k->com);
    if (k->lin) qfits_free(k->lin);
    k->val = NULL ;
    k->com = NULL ;
    k->lin = NULL ;
    if (val) {
        if (strlen(val)>0) k->val = qfits_strdup(val);
    }
    if (com) {
        if (strlen(com)>0) k->com = qfits_strdup(com);
    }
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Sort a FITS header
  @param    hdr     Header to sort (modified)
  @return   -1 in error case, 0 otherwise
 */
/*----------------------------------------------------------------------------*/
int qfits_header_sort(qfits_header ** hdr)
{
    qfits_header    *   sorted ;
    keytuple        *   k ;
    keytuple        *   kbf ;
    keytuple        *   next ;
    keytuple        *   last ;

    /* Test entries */
    if (hdr == NULL) return -1 ;
    if (*hdr == NULL) return -1 ;
    if ((*hdr)->n < 2) return 0 ;

    /* Create the new FITS header */
    sorted = qfits_header_new() ;

    /* Move the first keytuple to the sorted empty header */
    k = (keytuple*)(*hdr)->first ;
    next = k->next ;
    sorted->first = sorted->last = k ;
    k->next = k->prev = NULL ;
    sorted->n = 1 ;

    /* Loop over the other tuples */
    while (next != NULL) {
        k = next ;
        next = k->next ;

        /* Find k's place in sorted */
        kbf = (keytuple*)sorted->first ;
        while (kbf!=NULL) {
            if (k->typ < kbf->typ) break ;
            kbf = kbf->next ;
        }

        /* Hook k into sorted list */
        if (kbf == NULL) {
            /* k is last in sorted */
            last = sorted->last ;
            sorted->last = k ;
            k->next = NULL ;
            k->prev = last ;
            last->next = k ;
        } else {
            /* k goes just before kbf */
            k->next = kbf ;
            k->prev = kbf->prev ;
            if (kbf->prev != NULL) (kbf->prev)->next = k ;
            else sorted->first = k ;
            kbf->prev = k ;
        }
        (sorted->n) ++ ;
    }

    /* Replace the input header by the sorted one */
    (*hdr)->first = (*hdr)->last = NULL ;
    qfits_header_destroy(*hdr) ;
    *hdr = sorted ;

    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Copy a FITS header
  @param    src    Header to replicate
  @return    Pointer to newly allocated qfits_header object.

  Makes a strict copy of all information contained in the source
  header. The returned header must be freed using qfits_header_destroy.
 */
/*----------------------------------------------------------------------------*/
qfits_header * qfits_header_copy(const qfits_header * src)
{
    qfits_header    *   fh_copy ;
    keytuple        *   k ;

    if (src==NULL) return NULL ;

    fh_copy = qfits_header_new();
    k = (keytuple*)src->first ;
    while (k!=NULL) {
        qfits_header_append(fh_copy, k->key, k->val, k->com, k->lin) ;
        k = k->next ;
    }
    return fh_copy ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    qfits_header destructor
  @param    hdr qfits_header to deallocate
  @return    void

  Frees all memory associated to a given qfits_header object.
 */
/*----------------------------------------------------------------------------*/
void qfits_header_destroy(qfits_header * hdr)
{
    keytuple * k ;
    keytuple * kn ;

    if (hdr==NULL) return ;

    k = (keytuple*)hdr->first ;
    while (k!=NULL) {
        kn = k->next ;
        keytuple_del(k);
        k = kn ;
    }
    qfits_free(hdr);
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Return the value associated to a key, as a string
  @param    hdr qfits_header to parse
  @param    key key to find
  @return    pointer to statically allocated string

  Finds the value associated to the given key and return it as a
  string. The returned pointer is statically allocated, so do not
  modify its contents or try to free it.

  Returns NULL if no matching key is found or no value is attached.
 */
/*----------------------------------------------------------------------------*/
char * qfits_header_getstr(const qfits_header * hdr, const char * key)
{
    keytuple    *   k ;
    char        *   xkey ;

    if (hdr==NULL || key==NULL) return NULL ;

    xkey = qfits_expand_keyword(key);
    k = (keytuple*)hdr->first ;
    while (k!=NULL) {
        if (!strcmp(k->key, xkey)) break ;
        k=k->next ;
    }
    if (k==NULL) return NULL ;
    return k->val ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Return the i-th key/val/com/line tuple in a header.
  @param    hdr        Header to consider
  @param    idx        Index of the requested card
  @param    key        Output key
  @param    val        Output value
  @param    com        Output comment
  @param    lin        Output initial line
  @return    int 0 if Ok, -1 if error occurred.

  This function is useful to browse a FITS header object card by card.
  By iterating on the number of cards (available in the 'n' field of
  the qfits_header struct), you can retrieve the FITS lines and their
  components one by one. Indexes run from 0 to n-1. You can pass NULL
  values for key, val, com or lin if you are not interested in a
  given field.

  @code
  int i ;
  char key[FITS_LINESZ+1] ;
  char val[FITS_LINESZ+1] ;
  char com[FITS_LINESZ+1] ;
  char lin[FITS_LINESZ+1] ;

  for (i=0 ; i<hdr->n ; i++) {
      qfits_header_getitem(hdr, i, key, val, com, lin);
    printf("card[%d] key[%s] val[%s] com[%s]\n", i, key, val, com);
  }
  @endcode

  This function has primarily been written to interface a qfits_header
  object to other languages (C++/Python). If you are working within a
  C program, you should use the other header manipulation routines
  available in this module.
 */
/*----------------------------------------------------------------------------*/
int qfits_header_getitem(
        const qfits_header  *   hdr,
        int                     idx,
        char                *   key,
        char                *   val,
        char                *   com,
        char                *   lin)
{
    keytuple    *   k ;
    int             count ;

    if (hdr==NULL) return -1 ;
    if (key==NULL && val==NULL && com==NULL && lin==NULL) return 0 ;
    if (idx<0 || idx>hdr->n) return -1 ;

    /* Get pointer to keytuple */
    if (idx == 0) {
	    ((qfits_header *)hdr)->current_idx = 0 ;
	    ((qfits_header *)hdr)->current = hdr->first ;
	    k = hdr->current ;
	} else if (idx == hdr->current_idx + 1) {
	    ((qfits_header *)hdr)->current = ((keytuple*) (hdr->current))->next ;
	    ((qfits_header *)hdr)->current_idx++ ;
	    k = hdr->current ;
	} else {
	    count=0 ;
	    k = (keytuple*)hdr->first ;
	    while (count<idx) {
            k = k->next ;
            count++ ;
        }
	}

    /* Fill return values */
    if (key!=NULL) strcpy(key, k->key);
    if (val!=NULL) {
        if (k->val!=NULL) strcpy(val, k->val);
        else val[0]=0 ;
    }
    if (com!=NULL) {
        if (k->com!=NULL) strcpy(com, k->com);
        else com[0]=0 ;
    }
    if (lin!=NULL) {
        if (k->lin!=NULL) strcpy(lin, k->lin);
        else lin[0]=0 ;
    }
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Return the comment associated to a key, as a string
  @param    hdr qfits_header to parse
  @param    key key to find
  @return    pointer to statically allocated string

  Finds the comment associated to the given key and return it as a
  string. The returned pointer is statically allocated, so do not
  modify its contents or try to free it.

  Returns NULL if no matching key is found or no comment is attached.
 */
/*----------------------------------------------------------------------------*/
char * qfits_header_getcom(const qfits_header * hdr, const char * key)
{
    keytuple    *   k ;
    char        *   xkey ;

    if (hdr==NULL || key==NULL) return NULL ;

    xkey = qfits_expand_keyword(key);
    k = (keytuple*)hdr->first ;
    while (k!=NULL) {
        if (!strcmp(k->key, xkey)) break ;
        k=k->next ;
    }
    if (k==NULL) return NULL ;
    return k->com ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Return the value associated to a key, as an int
  @param    hdr qfits_header to parse
  @param    key key to find
  @param    errval default value to return if nothing is found
  @return    int

  Finds the value associated to the given key and return it as an
  int. Returns errval if no matching key is found or no value is
  attached.
 */
/*----------------------------------------------------------------------------*/
int qfits_header_getint(
        const qfits_header  *   hdr,
        const char          *   key,
        int                     errval)
{
    char    *   c ;
    int         d ;

    if (hdr==NULL || key==NULL) return errval ;

    c = qfits_header_getstr(hdr, key);
    if (c==NULL) return errval ;
    if (sscanf(c, "%d", &d)!=1) return errval ;
    return d ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Return the value associated to a key, as a double
  @param    hdr qfits_header to parse
  @param    key key to find
  @param    errval default value to return if nothing is found
  @return    double

  Finds the value associated to the given key and return it as a
  double. Returns errval if no matching key is found or no value is
  attached.
 */
/*----------------------------------------------------------------------------*/
double qfits_header_getdouble(
        const qfits_header  *   hdr,
        const char          *   key,
        double                  errval)
{
    char    *    c ;

    if (hdr==NULL || key==NULL) return errval ;

    c = qfits_header_getstr(hdr, key);
    if (c==NULL) return errval ;
    return atof(c);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Return the value associated to a key, as a boolean (int).
  @param    hdr qfits_header to parse
  @param    key key to find
  @param    errval default value to return if nothing is found
  @return    int

  Finds the value associated to the given key and return it as a
  boolean. Returns errval if no matching key is found or no value is
  attached. A boolean is here understood as an int taking the value 0
  or 1. errval can be set to any other integer value to reflect that
  nothing was found.

  errval is returned if no matching key is found or no value is
  attached.

  A true value is any character string beginning with a 'y' (yes), a
  't' (true) or the digit '1'. A false value is any character string
  beginning with a 'n' (no), a 'f' (false) or the digit '0'.
 */
/*----------------------------------------------------------------------------*/
int qfits_header_getboolean(
        const qfits_header  *   hdr,
        const char          *   key,
        int                     errval)
{
    char    *    c ;
    int            ret ;

    if (hdr==NULL || key==NULL) return errval ;

    c = qfits_header_getstr(hdr, key);
    if (c==NULL) return errval ;
    if (strlen(c)<1) return errval ;

    if (c[0]=='y' || c[0]=='Y' || c[0]=='1' || c[0]=='t' || c[0]=='T') {
        ret = 1 ;
    } else if (c[0]=='n' || c[0]=='N' || c[0]=='0' || c[0]=='f' || c[0]=='F') {
        ret = 0 ;
    } else {
        ret = errval ;
    }
    return ret;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Dump a FITS header to an opened file.
  @param    hdr     FITS header to dump
  @param    out     Opened file pointer
  @return   int 0 if Ok, -1 otherwise
  Dumps a FITS header to an opened file pointer.
 */
/*----------------------------------------------------------------------------*/
int qfits_header_dump(
        const qfits_header  *   hdr,
        FILE                *   out)
{
    keytuple    *   k ;
    char            line[81];
    int             n_out ;

    if (hdr==NULL) return -1 ;
    if (out==NULL) out=stdout ;

    k = (keytuple*)hdr->first ;
    n_out = 0 ;
    while (k!=NULL) {
        /* Make line from information in the node */
        qfits_header_makeline(line, k, 1);
        if ((fwrite(line, 1, 80, out))!=80) {
            fprintf(stderr, "error dumping FITS header");
            return -1 ;
        }
        n_out ++;
        k=k->next;
    }
    /* If printing out to a regular file, blank pad */
    if (out!=stdout && out!=stderr) {
        /* Blank-pad the output */
        memset(line, ' ', 80);
        while (n_out % 36) {
            fwrite(line, 1, 80, out);
            n_out++ ;
        }
    }
    return 0 ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    keytuple constructor
  @param    key        Key associated to key tuple (cannot be NULL).
  @param    val        Value associated to key tuple.
  @param    com        Comment associated to key tuple.
  @param    lin        Initial line read from FITS header (if applicable).
  @return    1 pointer to newly allocated keytuple.

  This function is a keytuple creator. NULL values and zero-length strings
  are valid parameters for all but the key field. The returned object must
  be deallocated using keytuple_del().

 */
/*----------------------------------------------------------------------------*/
static keytuple * keytuple_new(
        const char * key,
        const char * val,
        const char * com,
        const char * lin)
{
    keytuple    *    k ;

    if (key==NULL) return NULL ;

    /* Allocate space for new structure */
    k = qfits_malloc(sizeof(keytuple));
    /* Hook a copy of the new key */
    k->key = qfits_strdup(qfits_expand_keyword(key)) ;
    /* Hook a copy of the value if defined */
    k->val = NULL ;
    if (val!=NULL) {
        if (strlen(val)>0) k->val = qfits_strdup(val);
    }
    /* Hook a copy of the comment if defined */
    k->com = NULL ;
    if (com!=NULL) {
        if (strlen(com)>0) k->com = qfits_strdup(com) ;
    }
    /* Hook a copy of the initial line if defined */
    k->lin = NULL ;
    if (lin!=NULL) {
        if (strlen(lin)>0) k->lin = qfits_strdup(lin);
    }
    k->next = NULL ;
    k->prev = NULL ;
    k->typ = keytuple_type(key);

    return k;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    keytuple type identification routine
  @param    key        String representing a FITS keyword.
  @return    A key type (see keytype enum)

  This function identifies the type of a FITS keyword when given the
  keyword as a string. Keywords are expected literally as they are
  found in a FITS header on disk.

 */
/*----------------------------------------------------------------------------*/
static keytype keytuple_type(const char * key)
{
    keytype kt ;

    kt = keytype_undef ;
    /* Assign type to key tuple */
    if (!strcmp(key, "SIMPLE") || !strcmp(key, "XTENSION")) kt = keytype_top ;
    else if (!strcmp(key, "END"))                   kt = keytype_end ;
    else if (!strcmp(key, "BITPIX"))                kt = keytype_bitpix ;
    else if (!strcmp(key, "NAXIS"))                 kt = keytype_naxis ;
    else if (!strcmp(key, "NAXIS1"))                kt = keytype_naxis1 ;
    else if (!strcmp(key, "NAXIS2"))                kt = keytype_naxis2 ;
    else if (!strcmp(key, "NAXIS3"))                kt = keytype_naxis3 ;
    else if (!strcmp(key, "NAXIS4"))                kt = keytype_naxis4 ;
    else if (!strncmp(key, "NAXIS", 5))             kt = keytype_naxisi ;
    else if (!strcmp(key, "GROUP"))                 kt = keytype_group ;
    else if (!strcmp(key, "PCOUNT"))                kt = keytype_pcount ;
    else if (!strcmp(key, "GCOUNT"))                kt = keytype_gcount ;
    else if (!strcmp(key, "EXTEND"))                kt = keytype_extend ;
    else if (!strcmp(key, "BSCALE"))                kt = keytype_bscale ;
    else if (!strcmp(key, "BZERO"))                 kt = keytype_bzero ;
    else if (!strcmp(key, "TFIELDS"))               kt = keytype_tfields ;
    else if (!strncmp(key, "TBCOL", 5))             kt = keytype_tbcoli ;
    else if (!strncmp(key, "TFORM", 5))             kt = keytype_tformi ;
    else if (!strncmp(key, "HIERARCH ESO DPR", 16)) kt = keytype_hierarch_dpr ;
    else if (!strncmp(key, "HIERARCH ESO OBS", 16)) kt = keytype_hierarch_obs ;
    else if (!strncmp(key, "HIERARCH ESO TPL", 16)) kt = keytype_hierarch_tpl ;
    else if (!strncmp(key, "HIERARCH ESO GEN", 16)) kt = keytype_hierarch_gen ;
    else if (!strncmp(key, "HIERARCH ESO TEL", 16)) kt = keytype_hierarch_tel ;
    else if (!strncmp(key, "HIERARCH ESO INS", 16)) kt = keytype_hierarch_ins ;
    else if (!strncmp(key, "HIERARCH ESO LOG", 16)) kt = keytype_hierarch_log ;
    else if (!strncmp(key, "HIERARCH ESO PRO", 16)) kt = keytype_hierarch_pro ;
    else if (!strncmp(key, "HIERARCH", 8))          kt = keytype_hierarch ;
    else if (!strcmp(key, "HISTORY"))               kt = keytype_history ;
    else if (!strcmp(key, "COMMENT"))               kt = keytype_comment ;
    else if ((int)strlen(key)<9)                    kt = keytype_primary ;
    return kt ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Keytuple destructor.
  @param    k    Keytuple to deallocate.
  @return    void
  Keytuple destructor.
 */
/*----------------------------------------------------------------------------*/
static void keytuple_del(keytuple * k)
{
    if (k==NULL) return ;
    if (k->key) qfits_free(k->key);
    if (k->val) qfits_free(k->val);
    if (k->com) qfits_free(k->com);
    if (k->lin) qfits_free(k->lin);
    qfits_free(k);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Keytuple dumper.
  @param    k    Keytuple to dump
  @return    void

  This function dumps a key tuple to stdout. It is meant for debugging
  purposes only.
 */
/*----------------------------------------------------------------------------*/
static void keytuple_dmp(const keytuple * k)
{
    if (!k) return ;
    printf("[%s]=[", k->key);
    if (k->val) printf("%s", k->val);
    printf("]");
    if (k->com) printf("/[%s]", k->com);
    printf("\n");
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Build a FITS line from the information contained in a card.
  @param    line pointer to allocated string to be filled
  @param    node pointer to card node in qfits_header linked-list
  @param    conservative flag to indicate conservative behaviour
  @return    int 0 if Ok, anything else otherwise

  Build a FITS line (80 chars) from the information contained in a
  node of a qfits_header linked-list. If the mode is set to
  conservative, the original FITS line will be used wherever present.
  If conservative is set to 0, a new line will be formatted.
 */
/*----------------------------------------------------------------------------*/
static int qfits_header_makeline(
        char            *   line,
        const keytuple  *   k,
        int                 conservative)
{
    char blankline[81];
    int     i ;

    if (line==NULL || k==NULL) return -1 ;

    /* If a previous line information is there, use it as is */
    if (conservative) {
        if (k->lin != NULL) {
            memcpy(line, k->lin, 80);
            line[80]=(char)0;
            return 0 ;
        }
    }
    /* Got to build keyword from scratch */
    memset(blankline, 0, 81);
    qfits_card_build(blankline, k->key, k->val, k->com);
    memset(line, ' ', 80);
    i=0 ;
    while (blankline[i] != (char)0) {
        line[i] = blankline[i] ;
        i++ ;
    }
    line[80]=(char)0;
    return 0;
}
#endif
#ifndef QFITS_IMAGE_H
#define QFITS_IMAGE_H
/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/










/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/** Symbol to set returned pixel type to float */
#define PTYPE_FLOAT        0
/** Symbol to set returned pixel type to int */
#define PTYPE_INT        1
/** Symbol to set returned pixel type to double */
#define PTYPE_DOUBLE    2

/* FITS pixel depths */
/* FITS BITPIX=8 */
#define BPP_8_UNSIGNED        (8)
/* FITS BITPIX=16 */
#define BPP_16_SIGNED        (16)
/* FITS BITPIX=32 */
#define BPP_32_SIGNED        (32)
/* FITS BITPIX=-32 */
#define BPP_IEEE_FLOAT      (-32)
/* FITS BITPIX=-64 */
#define BPP_IEEE_DOUBLE     (-64)
/* Default BITPIX for output */
#define BPP_DEFAULT         BPP_IEEE_FLOAT

/*-----------------------------------------------------------------------------
                                   New types
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Alias for unsigned char

/*----------------------------------------------------------------------------*/
/**
  @brief    qfits loader control object

  This structure serves two purposes: input and output for the qfits
  pixel loading facility. To request pixels from a FITS file, you
  need to allocate (statically or dynamically) such a structure and
  fill up the input fields (filename, xtension number, etc.) to specify
  the pixels you want from the file.

  Before performing the actual load, you must pass the initialized
  structure to qfitsloader_init() which will check whether the operation
  is feasible or not (check its returned value).

  If the operation was deemed feasible, you can proceed to load the pixels,
  passing the same structure to qfits_loadpix() which will fill up the
  output fields of the struct. Notice that a pixel buffer will have been
  allocated (through malloc or mmap) and placed into the structure. You
  need to call free() on this pointer when you are done with it,
  typically in the image or cube destructor.

  The qfitsloader_init() function is also useful to probe a FITS file
  for useful informations, like getting the size of images in the file,
  the pixel depth, or data offset.

  Example of a code that prints out various informations about
  a plane to load, without actually loading it:

  @code
int main(int argc, char * argv[])
{
    qfitsloader    ql ;

    ql.filename = argv[1] ;
    ql.xtnum    = 0 ;
    ql.pnum     = 0 ;

    if (qfitsloader_init(&ql)!=0) {
        printf("cannot read info about %s\n", argv[1]);
        return -1 ;
    }

    printf(    "file         : %s\n"
            "xtnum        : %d\n"
            "pnum         : %d\n"
            "# xtensions  : %d\n"
            "size X       : %d\n"
            "size Y       : %d\n"
            "planes       : %d\n"
            "bitpix       : %d\n"
            "datastart    : %d\n"
            "datasize     : %d\n"
            "bscale       : %g\n"
            "bzero        : %g\n",
            ql.filename,
            ql.xtnum,
            ql.pnum,
            ql.exts,
            ql.lx,
            ql.ly,
            ql.np,
            ql.bitpix,
            ql.seg_start,
            ql.seg_size,
            ql.bscale,
            ql.bzero);
    return 0 ;
}
  @endcode

/**
  @brief    qfits dumper control object

  This structure offers various control parameters to dump a pixel
  buffer to a FITS file. The buffer will be dumped as requested
  to the requested file in append mode. Of course, the requested file
  must be writeable for the operation to succeed.

  The following example demonstrates how to save a linear ramp sized
  100x100 to a FITS file with BITPIX=16. Notice that this code only
  dumps the pixel buffer, no header information is provided in this
  case.

  @code
    int   i, j ;
    int * ibuf ;
    qfitsdumper    qd ;

    // Fill a buffer with 100x100 int pixels
    ibuf = malloc(100 * 100 * sizeof(int));
    for (j=0 ; j<100 ; j++) {
        for (i=0 ; i<100 ; i++) {
            ibuf[i+j*100] = i+j ;
        }
    }

    qd.filename  = "out.fits" ;     // Output file name
    qd.npix      = 100 * 100 ;      // Number of pixels
    qd.ptype     = PTYPE_INT ;      // Input buffer type
    qd.ibuf      = ibuf ;           // Set buffer pointer
    qd.out_ptype = BPP_16_SIGNED ;  // Save with BITPIX=16

    // Dump buffer to file (error checking omitted for clarity)
    qfits_pixdump(&qd);

    free(ibuf);
  @endcode

  If the provided output file name is "STDOUT" (all capitals), the
  function will dump the pixels to the stdout steam (usually the console,
  could have been re-directed).
 */
/*----------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------------
                                Defines
 -----------------------------------------------------------------------------*/

#define QFITSLOADERINIT_MAGIC   0xcafe
/* Compute the number of bytes per pixel for a given BITPIX value */
#define BYTESPERPIXEL(x)    (   ((x) == BPP_8_UNSIGNED) ?     1 : \
                                ((x) == BPP_16_SIGNED)  ?     2 : \
                                ((x) == BPP_32_SIGNED)  ?     4 : \
                                ((x) == BPP_IEEE_FLOAT) ?     4 : \
                                ((x) == BPP_IEEE_DOUBLE) ?    8 : 0 )

/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_image Pixel loader for FITS images.
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                            Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Initialize a qfitsloader control object.
  @param    ql  qfitsloader object to initialize.
  @return   int 0 if Ok, -1 if error occurred.

  This function expects a qfitsloader object with a number of input
  fields correctly filled in. The minimum fields to set are:

  - filename: Name of the file to examine.
  - xtnum: Extension number to examine (0 for main section).
  - pnum: Plane number in the requested extension.
  - map : loading mode - flag to know if the file has to be mapped

  You can go ahead with these fields only if you only want to get
  file information for this plane in this extension. If you want
  to later load the plane, you must additionally fill the 'ptype'
  field to a correct value (PTYPE_INT, PTYPE_FLOAT, PTYPE_DOUBLE)
  before calling qfits_loadpix() so that it knows which conversion
  to perform.

  This function is basically a probe sent on a FITS file to ask
  qfits if loading these data would be Ok or not. The actual loading
  is performed by qfits_loadpix() afterwards.
 */
/*----------------------------------------------------------------------------*/
int qfitsloader_init(qfitsloader * ql)
{
    qfits_header    *   fh ;

    int     n_ext ;
    int     seg_start ;
    int     seg_size ;
    int     bitpix, naxis, naxis1, naxis2, naxis3 ;
    char *  xt_type ;
    char *  sval ;
    struct stat sta ;

    /* Check passed object is allocated */
    if (ql==NULL) {
        qfits_error("pixio: NULL loader");
        return -1 ;
    }
    /* Object must contain a filename */
    if (ql->filename == NULL) {
        qfits_error("pixio: NULL filename in loader");
        return -1 ;
    }
    /* Check requested file exists and contains data */
    if (stat(ql->filename, &sta)!=0) {
        qfits_error("no such file: %s", ql->filename);
        return -1 ;
    }
    if (sta.st_size<1) {
        qfits_error("empty file: %s", ql->filename);
        return -1 ;
    }

    /* Requested extension number must be positive */
    if (ql->xtnum<0) {
        qfits_error("pixio: negative xtnum in loader");
        return -1 ;
    }
    /* Requested returned pixel type must be legal */
    if ((ql->ptype!=PTYPE_INT) &&
        (ql->ptype!=PTYPE_FLOAT) &&
        (ql->ptype!=PTYPE_DOUBLE)) {
        qfits_error("pixio: invalid ptype in loader");
        return -1 ;
    }
    /* Check requested file is FITS */
    if (qfits_is_fits(ql->filename)!=1) {
        qfits_error("pixio: not a FITS file: %s", ql->filename);
        return -1 ;
    }
    /* Get number of extensions for this file */
    n_ext = qfits_query_n_ext(ql->filename);
    if (n_ext==-1) {
        qfits_error("pixio: cannot get number of extensions in %s",
                    ql->filename);
        return -1 ;
    }
    /* Check requested extension falls within range */
    if (ql->xtnum > n_ext) {
        qfits_error("pixio: requested extension %d but file %s has %d\n",
                    ql->xtnum,
                    ql->filename,
                    n_ext);
        return -1 ;
    }
    ql->exts = n_ext ;
    /* Get segment offset and size for the requested buffer */
    if (qfits_get_datinfo(ql->filename,
                          ql->xtnum,
                          &seg_start,
                          &seg_size)!=0) {
        qfits_error("pixio: cannot get seginfo for %s extension %d",
                    ql->filename,
                    ql->xtnum);
        return -1 ;
    }
    /* Check segment size is consistent with file size */
    if (sta.st_size < (seg_start+seg_size)) {
        return -1 ;
    }
    ql->seg_start = seg_start ;
    ql->seg_size  = seg_size ;

    /* Get file header */
    fh = qfits_header_readext(ql->filename, ql->xtnum);
    if (fh==NULL) {
        qfits_error("pixio: cannot read header from ext %d in %s",
                    ql->xtnum,
                    ql->filename);
        return -1 ;
    }
    /* If the requested image is within an extension */
    if (ql->xtnum>0) {
        /* Check extension is an image */
        xt_type = qfits_header_getstr(fh, "XTENSION");
        if (xt_type==NULL) {
            qfits_error("pixio: cannot read extension type for ext %d in %s",
                        ql->xtnum,
                        ql->filename);
            qfits_header_destroy(fh);
            return -1 ;
        }
        xt_type = qfits_pretty_string(xt_type);
        if (strcmp(xt_type, "IMAGE")) {
            qfits_error(
                "pixio: not an image -- extension %d in %s has type [%s]",
                ql->xtnum,
                ql->filename,
                xt_type);
            qfits_header_destroy(fh);
            return -1 ;
        }
    }

    /* Get file root informations */
    bitpix = qfits_header_getint(fh, "BITPIX", -1);
    naxis  = qfits_header_getint(fh, "NAXIS",  -1);
    naxis1 = qfits_header_getint(fh, "NAXIS1", -1);
    naxis2 = qfits_header_getint(fh, "NAXIS2", -1);
    naxis3 = qfits_header_getint(fh, "NAXIS3", -1);

    /* Get BSCALE and BZERO if available */
    sval = qfits_header_getstr(fh, "BSCALE");
    if (sval==NULL) {
        ql->bscale = 1.0 ;
    } else {
        ql->bscale = atof(sval);
    }
    sval = qfits_header_getstr(fh, "BZERO");
    if (sval==NULL) {
        ql->bzero = 0.0 ;
    } else {
        ql->bzero = atof(sval);
    }

    /* Destroy header */
    qfits_header_destroy(fh);

    /* Check BITPIX is present */
    if (bitpix==-1) {
        qfits_error("pixio: missing BITPIX in file %s", ql->filename);
        return -1 ;
    }
    /* Check BITPIX is valid */
    if ((bitpix!=   8) &&
        (bitpix!=  16) &&
        (bitpix!=  32) &&
        (bitpix!= -32) &&
        (bitpix!= -64)) {
        qfits_error("pixio: invalid BITPIX (%d) in file %s",
                    bitpix,
                    ql->filename);
        return -1 ;
    }
    ql->bitpix = bitpix ;

    /* Check NAXIS is present and valid */
    if (naxis<0) {
        qfits_error("pixio: no NAXIS in file %s", ql->filename);
        return -1 ;
    }
    if (naxis==0) {
        qfits_error("pixio: no pixel in ext %d of file %s",
                    ql->xtnum,
                    ql->filename);
        return -1 ;
    }
    if (naxis>3) {
        qfits_error("pixio: NAXIS>3 (%d) unsupported", naxis);
        return -1 ;
    }
    /* NAXIS1 must always be present */
    if (naxis1<0) {
        qfits_error("pixio: no NAXIS1 in file %s", ql->filename);
        return -1 ;
    }
    /* Update dimension fields in loader struct */
    ql->lx = 1 ;
    ql->ly = 1 ;
    ql->np = 1 ;

    switch (naxis) {
        case 1:
        ql->lx = naxis1 ;
        break ;

        case 3:
        if (naxis3<0) {
            qfits_error("pixio: no NAXIS3 in file %s", ql->filename);
            return -1 ;
        }
        ql->np = naxis3 ;
        /* No break statement: intended to go through next case */

        case 2:
        if (naxis2<0) {
            qfits_error("pixio: no NAXIS2 in file %s", ql->filename);
            return -1 ;
        }
        ql->ly = naxis2 ;
        ql->lx = naxis1 ;
        break ;
    }
    /* Check that requested plane number falls within range */
    if (ql->pnum >= ql->np) {
        qfits_error("pixio: requested plane %d but NAXIS3=%d",
                    ql->pnum,
                    ql->np);
        return -1 ;
    }

    /* Everything Ok, fields have been filled along. */
    /* Mark the structure as initialized */
    ql->_init = QFITSLOADERINIT_MAGIC ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a pixel buffer for one complete image.
  @param    ql  Allocated and initialized qfitsloader control object.
  @return   int 0 if Ok, -1 if error occurred.
  @see      qfits_loadpix_window
*/
/*----------------------------------------------------------------------------*/
int qfits_loadpix(qfitsloader * ql)
{
    if (ql==NULL) return -1 ;
    return qfits_loadpix_window(ql, 1, 1, ql->lx, ql->ly) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a pixel buffer for one image window
  @param    ql  Allocated and initialized qfitsloader control object.
  @param    llx
  @param    lly     Position of the window (start with (1,1))
  @param    urx
  @param    ury
  @return   int 0 if Ok, -1 if error occurred.

  This function performs a load of a pixel buffer into memory. It
  expects an allocated and initialized qfitsloader object in input.
  See qfitsloader_init() about initializing the object.

  This function will fill up the ibuf/fbuf/dbuf field, depending
  on the requested pixel type (resp. int, float or double).

  If llx lly urx and ury do not specify the whole image, ql->map must be 0,
  we do not want to mmap a file an load only a part of it.
 */
/*----------------------------------------------------------------------------*/
int qfits_loadpix_window(
        qfitsloader     *   ql,
        int                 llx,
        int                 lly,
        int                 urx,
        int                 ury)
{
    byte    *   fptr ;
    size_t      fsize ;
    int         datastart ;
    int         imagesize, window_size, linesize ;
    FILE    *   lfile ;
    int         dataread ;
    int         nx, ny ;
    int         i ;

    /* Check inputs */
    if (ql==NULL) return -1 ;
    if (ql->_init != QFITSLOADERINIT_MAGIC) {
        qfits_error("pixio: called with unitialized obj");
        return -1 ;
    }
    if (llx>urx || lly>ury || llx<1 || lly<1 || urx>ql->lx || ury>ql->ly) {
        qfits_error("pixio: invalid window specification");
        return -1 ;
    }

    /* No map if only a zone is specified */
    if (llx != 1 || lly != 1 || urx != ql->lx || ury != ql->ly) {
        if (ql->map == 1) {
            qfits_error("pixio: cannot mmap for a part of the image");
            return -1 ;
        }
    }

    /* Initialise */
    nx = urx-llx+1 ;
    ny = ury-lly+1 ;
    imagesize = ql->lx * ql->ly * BYTESPERPIXEL(ql->bitpix);
    window_size = nx * ny * BYTESPERPIXEL(ql->bitpix);
    linesize = nx * BYTESPERPIXEL(ql->bitpix);
    datastart = ql->seg_start + ql->pnum * imagesize ;

    /* Check loading mode */
    if (ql->map) {
        /* Map the file in */
        fptr = (byte*)qfits_falloc(ql->filename, datastart, &fsize);
        if (fptr==NULL) {
            qfits_error("pixio: cannot falloc(%s)", ql->filename);
            return -1 ;
        }
    } else {
        /* Open the file */
        if ((lfile=fopen(ql->filename, "r"))==NULL) {
            qfits_error("pixio: cannot open %s", ql->filename);
            return -1 ;
        }
        /* Go to the start of the image */
        if (fseek(lfile, datastart, SEEK_SET)!=0) {
            qfits_error("pixio: cannot seek %s", ql->filename);
            fclose(lfile);
            return -1 ;
        }
        /* Go to the start of the zone */
        if (fseek(lfile, (llx-1+(lly-1)*ql->lx)*BYTESPERPIXEL(ql->bitpix),
                    SEEK_CUR)!=0) {
            qfits_error("pixio: cannot seek %s", ql->filename);
            fclose(lfile);
            return -1 ;
        }

        fptr = (byte*)qfits_malloc(window_size) ;

        /* Only a window is specified */
        if (llx != 1 || lly != 1 || urx != ql->lx || ury != ql->ly) {
            /* Loop on the lines */
            for (i=0 ; i<ny ; i++) {
                /* Read the file */
                dataread=fread(fptr+i*linesize, sizeof(byte), linesize, lfile);
                if (dataread!=linesize) {
                    qfits_free(fptr) ;
                    fclose(lfile);
                    qfits_error("pixio: cannot read from %s", ql->filename);
                    return -1 ;
                }
                /* Go to the next line */
                if (fseek(lfile,ql->lx*BYTESPERPIXEL(ql->bitpix)-linesize,
                           SEEK_CUR)!=0){
                    qfits_error("pixio: cannot seek %s", ql->filename);
                    fclose(lfile);
                    return -1 ;
                }
            }
            fclose(lfile);
        } else {
        /* The whole image is specified */
            dataread = fread(fptr, sizeof(byte), window_size, lfile);
            fclose(lfile);
            if (dataread!=window_size) {
                qfits_free(fptr) ;
                qfits_error("pixio: cannot read from %s", ql->filename);
                return -1 ;
            }
        }
    }

    /* Initialize buffer pointers */
    ql->ibuf = NULL ;
    ql->fbuf = NULL ;
    ql->dbuf = NULL ;

    /*
     * Special cases: mapped file is identical to requested format.
     * This is only possible on big-endian machines, since FITS is
     * big-endian only.
     */
#ifdef WORDS_BIGENDIAN
    if (ql->ptype==PTYPE_FLOAT && ql->bitpix==-32) {
        ql->fbuf = (float*)fptr ;
        return 0 ;
    }
    if (ql->ptype==PTYPE_DOUBLE && ql->bitpix==-64) {
        ql->dbuf = (double*)fptr ;
        return 0 ;
    }
    if (ql->ptype==PTYPE_INT && ql->bitpix==32) {
        ql->ibuf = (int*)fptr ;
        return 0 ;
    }
#endif

    /* General case: fallback to dedicated conversion function */
    switch (ql->ptype) {
        case PTYPE_FLOAT:
        ql->fbuf = qfits_pixin_float(   fptr,
                                        nx * ny,
                                        ql->bitpix,
                                        ql->bscale,
                                        ql->bzero);
        break ;

        case PTYPE_INT:
        ql->ibuf = qfits_pixin_int( fptr,
                                    nx * ny,
                                    ql->bitpix,
                                    ql->bscale,
                                    ql->bzero);
        break ;

        case PTYPE_DOUBLE:
        ql->dbuf = qfits_pixin_double(  fptr,
                                        nx * ny,
                                        ql->bitpix,
                                        ql->bscale,
                                        ql->bzero);
        break ;
    }

    if (ql->map) {
        qfits_fdealloc((char*)fptr, datastart, fsize) ;
    } else {
        qfits_free(fptr);
    }

    if (ql->ibuf==NULL && ql->fbuf==NULL && ql->dbuf==NULL) {
        qfits_error("pixio: error during conversion");
        return -1 ;
    }
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Dump a pixel buffer to an output FITS file in append mode.
  @param    qd  qfitsdumper control object.
  @return   int 0 if Ok, -1 otherwise.

  This function takes in input a qfitsdumper control object. This object
  must be allocated beforehand and contain valid references to the data
  to save, and how to save it.

  The minimum fields to fill are:

  - filename: Name of the FITS file to dump to.
  - npix: Number of pixels in the buffer to be dumped.
  - ptype: Type of the passed buffer (PTYPE_FLOAT, PTYPE_INT, PTYPE_DOUBLE)
  - out_ptype: Requested FITS BITPIX for the output.

  One of the following fields must point to the corresponding pixel
  buffer:

  - ibuf for an int pixel buffer (ptype=PTYPE_INT)
  - fbuf for a float pixel buffer (ptype=PTYPE_FLOAT)
  - dbuf for a double pixel buffer (ptype=PTYPE_DOUBLE)

  This is a fairly low-level function, in the sense that it does not
  check that the output file already contains a proper header or even
  that the file it is appending to is indeed a FITS file. It will
  convert the pixel buffer to the requested BITPIX type and append
  the data to the file, without padding with zeros. See qfits_zeropad()
  about padding.

  If the given output file name is "STDOUT" (all caps), the dump
  will be performed to stdout.
 */
/*----------------------------------------------------------------------------*/
int qfits_pixdump(qfitsdumper * qd)
{
    FILE    *   f_out ;
    byte    *   buf_out ;
    int         buf_free ;
    int         buf_sz ;

    /* Check inputs */
    if (qd==NULL) return -1 ;
    if (qd->filename==NULL) return -1 ;
    switch (qd->ptype) {
        case PTYPE_FLOAT:
        if (qd->fbuf==NULL) return -1 ;
        break ;

        case PTYPE_DOUBLE:
        if (qd->dbuf==NULL) return -1 ;
        break ;

        case PTYPE_INT:
        if (qd->ibuf==NULL) return -1 ;
        break ;

        default:
        return -1 ;
    }
    if (qd->npix <= 0) {
        qfits_error("Negative or NULL number of pixels specified");
        return -1 ;
    }

    /*
     * Special cases: input buffer is identical to requested format.
     * This is only possible on big-endian machines, since FITS is
     * big-endian only.
     */
    buf_out = NULL ;
    buf_free = 1 ;
#ifdef WORDS_BIGENDIAN
    if (qd->ptype==PTYPE_FLOAT && qd->out_ptype==-32) {
        buf_out = (byte*)qd->fbuf ;
        buf_free=0 ;
    } else if (qd->ptype==PTYPE_DOUBLE && qd->out_ptype==-64) {
        buf_out = (byte*)qd->dbuf ;
        buf_free=0 ;
    } else if (qd->ptype==PTYPE_INT && qd->out_ptype==32) {
        buf_out = (byte*)qd->ibuf ;
        buf_free=0 ;
    }
#endif
    buf_sz = qd->npix * BYTESPERPIXEL(qd->out_ptype);

    /* General case */
    if (buf_out==NULL) {
        switch (qd->ptype) {
            /* Convert buffer */
            case PTYPE_FLOAT:
            buf_out = qfits_pixdump_float(  qd->fbuf,
                                            qd->npix,
                                            qd->out_ptype);
            break ;

            /* Convert buffer */
            case PTYPE_INT:
            buf_out = qfits_pixdump_int(    qd->ibuf,
                                            qd->npix,
                                            qd->out_ptype);
            break ;

            /* Convert buffer */
            case PTYPE_DOUBLE:
            buf_out = qfits_pixdump_double( qd->dbuf,
                                             qd->npix,
                                             qd->out_ptype);
            break ;
        }
    }
    if (buf_out==NULL) {
        qfits_error("cannot dump pixel buffer");
        return -1 ;
    }

    /* Dump buffer */
    if (!strncmp(qd->filename, "STDOUT", 6)) {
        f_out = stdout ;
    } else {
        f_out = fopen(qd->filename, "a");
    }
    if (f_out==NULL) {
        qfits_free(buf_out);
        return -1 ;
    }
    fwrite(buf_out, buf_sz, 1, f_out);
    if (buf_free) {
        qfits_free(buf_out);
    }
    if (f_out!=stdout) {
        fclose(f_out);
    }
    return 0 ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert a float pixel buffer to a byte buffer.
  @param    buf     Input float buffer.
  @param    npix    Number of pixels in the input buffer.
  @param    ptype   Requested output BITPIX type.
  @return   1 pointer to a newly allocated byte buffer.

  This function converts the given float buffer to a buffer of bytes
  suitable for dumping to a FITS file (i.e. big-endian, in the
  requested pixel type). The returned pointer must be deallocated
  using the qfits_free() function.
 */
/*----------------------------------------------------------------------------*/
static byte * qfits_pixdump_float(float * buf, int npix, int ptype)
{
    byte    *   buf_out ;
    register byte * op ;
    int         i ;
    int         lpix ;
    short       spix ;
    double      dpix ;

    buf_out = qfits_malloc(npix * BYTESPERPIXEL(ptype));
    op = buf_out ;
    switch (ptype) {
        case 8:
        /* Convert from float to 8 bits */
        for (i=0 ; i<npix ; i++) {
            if (buf[i]>255.0) {
                *op++ = (byte)0xff ;
            } else if (buf[i]<0.0) {
                *op++ = (byte)0x00 ;
            } else {
                *op++ = (byte)buf[i];
            }
        }
        break ;

        case 16:
        /* Convert from float to 16 bits */
        for (i=0 ; i<npix ; i++) {
            if (buf[i]>32767.0) {
                *op++ = (byte)0x7f ;
                *op++ = (byte)0xff ;
            } else if (buf[i]<-32768.0) {
                *op++ = (byte)0x80 ;
                *op++ = 0x00 ;
            } else {
                spix = (short)buf[i];
                *op++ = (spix >> 8) ;
                *op++ = (spix & (byte)0xff) ;
            }
        }
        break ;

        case 32:
        /* Convert from float to 32 bits */
        for (i=0 ; i<npix ; i++) {
            if (buf[i] > 2147483647.0) {
                *op++ = (byte)0x7f ;
                *op++ = (byte)0xff ;
                *op++ = (byte)0xff ;
                *op++ = (byte)0xff ;
            } else if (buf[i]<-2147483648.0) {
                *op++ = (byte)0x80 ;
                *op++ = (byte)0x00 ;
                *op++ = (byte)0x00 ;
                *op++ = (byte)0x00 ;
            } else {
                lpix = (int)buf[i] ;
                *op++ = (byte)(lpix >> 24) ;
                *op++ = (byte)(lpix >> 16) & 0xff ;
                *op++ = (byte)(lpix >> 8 ) & 0xff ;
                *op++ = (byte)(lpix) & 0xff ;
            }
        }
        break ;

        case -32:
        /* Convert from float to float */
        memcpy(op, buf, npix * sizeof(float));
#ifndef WORDS_BIGENDIAN
        for (i=0 ; i<npix ; i++) {
            qfits_swap_bytes(op, 4);
            op++ ;
            op++ ;
            op++ ;
            op++ ;
        }
#endif
        break ;

        case -64:
        /* Convert from float to double */
        for (i=0 ; i<npix ; i++) {
            dpix = (double)buf[i] ;
#ifndef WORDS_BIGENDIAN
            qfits_swap_bytes(&dpix, 8) ;
#endif
            memcpy(op, &dpix, 8);
            op += 8 ;
        }
        break ;

        default:
        qfits_error("not supported yet");
        buf_out = NULL ;
        break ;
    }
    return buf_out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert an int pixel buffer to a byte buffer.
  @param    buf     Input int buffer.
  @param    npix    Number of pixels in the input buffer.
  @param    ptype   Requested output BITPIX type.
  @return   1 pointer to a newly allocated byte buffer.

  This function converts the given int buffer to a buffer of bytes
  suitable for dumping to a FITS file (i.e. big-endian, in the
  requested pixel type). The returned pointer must be deallocated
  using the qfits_free() function.
 */
/*----------------------------------------------------------------------------*/
static byte * qfits_pixdump_int(int * buf, int npix, int ptype)
{
    byte    *   buf_out ;
    register byte * op ;
    int i ;

    short   spix ;
    float   fpix ;
    double  dpix ;

    buf_out = qfits_malloc(npix * BYTESPERPIXEL(ptype));
    op = buf_out ;
    switch (ptype) {
        case 8:
        /* Convert from int32 to 8 bits */
        for (i=0 ; i<npix ; i++) {
            if (buf[i]>255) {
                *op++ = (byte)0xff ;
            } else if (buf[i]<0) {
                *op++ = (byte)0x00 ;
            } else {
                *op++ = (byte)buf[i] ;
            }
        }
        break ;

        case 16:
        /* Convert from int32 to 16 bits */
        for (i=0 ; i<npix ; i++) {
            if (buf[i]>32767) {
                spix = 32767 ;
            } else if (buf[i]<-32768) {
                spix = -32768 ;
            } else {
                spix = (short)buf[i] ;
            }
#ifndef WORDS_BIGENDIAN
            qfits_swap_bytes(&spix, 2);
#endif
            memcpy(op, &spix, 2);
            op += 2 ;
        }
        break ;

        case 32:
        /* Convert from int32 to 32 bits */
        memcpy(op, buf, npix * sizeof(int));
#ifndef WORDS_BIGENDIAN
        for (i=0 ; i<npix ; i++) {
            qfits_swap_bytes(op, 4);
            op+=4 ;
        }
#endif
        break ;

        case -32:
        /* Convert from int32 to float */
        for (i=0 ; i<npix ; i++) {
            fpix = (float)buf[i] ;
#ifndef WORDS_BIGENDIAN
            qfits_swap_bytes(&fpix, 4);
#endif
            memcpy(op, &fpix, 4) ;
            op += 4 ;
        }
        break ;

        case -64:
        /* Convert from int32 to double */
        for (i=0 ; i<npix ; i++) {
            dpix = (double)buf[i] ;
#ifndef WORDS_BIGENDIAN
            qfits_swap_bytes(&dpix, 8) ;
#endif
            memcpy(op, &dpix, 8);
            op += 8 ;
        }
        break ;

        default:
        qfits_error("not supported yet");
        buf_out = NULL ;
        break ;
    }
    return buf_out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert a double pixel buffer to a byte buffer.
  @param    buf     Input double buffer.
  @param    npix    Number of pixels in the input buffer.
  @param    ptype   Requested output BITPIX type.
  @return   1 pointer to a newly allocated byte buffer.

  This function converts the given double buffer to a buffer of bytes
  suitable for dumping to a FITS file (i.e. big-endian, in the
  requested pixel type). The returned pointer must be deallocated
  using the qfits_free() function.
 */
/*----------------------------------------------------------------------------*/
static byte * qfits_pixdump_double(double * buf, int npix, int ptype)
{
    byte    *   buf_out ;
    register byte  * op ;
    int i ;

    short   spix ;
    float   fpix ;
    int     lpix ;

    buf_out = qfits_malloc(npix * BYTESPERPIXEL(ptype));
    op = buf_out ;
    switch (ptype) {
        case 8:
        /* Convert from double to 8 bits */
        for (i=0 ; i<npix ; i++) {
            if (buf[i]>255.0) {
                *op++ = (byte)0xff ;
            } else if (buf[i]<0.0) {
                *op++ = (byte)0x00 ;
            } else {
                *op++ = (byte)buf[i] ;
            }
        }
        break ;

        case 16:
        /* Convert from double to 16 bits */
        for (i=0 ; i<npix ; i++) {
            if (buf[i]>32767.0) {
                spix = 32767 ;
            } else if (buf[i]<-32768.0) {
                spix = -32768 ;
            } else {
                spix = (short)buf[i] ;
            }
#ifndef WORDS_BIGENDIAN
            qfits_swap_bytes(&spix, 2);
#endif
            memcpy(op, &spix, 2);
            op += 2 ;
        }
        break ;

        case 32:
        /* Convert from double to 32 bits */
        for (i=0 ; i<npix ; i++) {
            if (buf[i] > 2147483647.0) {
                lpix = 2147483647 ;
            } else if (buf[i] < -2147483648.0) {
                lpix = -2147483647 ;
            } else {
                lpix = (int)buf[i];
            }
#ifndef WORDS_BIGENDIAN
            qfits_swap_bytes(&lpix, 4);
#endif
            memcpy(op, &lpix, 4);
            op += 4 ;
        }
        break ;

        case -32:
        /* Convert from double to float */
        for (i=0 ; i<npix ; i++) {
            fpix = (float)buf[i] ;
#ifndef WORDS_BIGENDIAN
            qfits_swap_bytes(&fpix, 4);
#endif
            memcpy(op, &fpix, 4) ;
            op += 4 ;
        }
        break ;

        case -64:
        /* Convert from double to double */
        memcpy(op, buf, npix * 8) ;
#ifndef WORDS_BIGENDIAN
        for (i=0 ; i<npix ; i++) {
            qfits_swap_bytes(op, 8);
            op += 8 ;
        }
#endif
        break ;

        default:
        qfits_error("not supported yet");
        buf_out = NULL ;
        break ;
    }
    return buf_out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a pixel buffer as floats.
  @param    p_source    Pointer to source buffer (as bytes).
  @param    npix        Number of pixels to load.
  @param    bitpix      FITS BITPIX in original file.
  @param    bscale      FITS BSCALE in original file.
  @param    bzero       FITS BZERO in original file.
  @return   1 pointer to a newly allocated buffer of npix floats.

  This function takes in input a pointer to a byte buffer as given
  in the original FITS file (big-endian format). It converts the
  buffer to an array of float (whatever representation is used for
  floats by this platform is used) and returns the newly allocated
  buffer, or NULL if an error occurred.

  The returned buffer must be deallocated using qfits_free().
 */
/*----------------------------------------------------------------------------*/
static float * qfits_pixin_float(
        byte    *   p_source,
        int         npix,
        int         bitpix,
        double      bscale,
        double      bzero)
{
    int         i ;
    float   *   baseptr ;
    float   *   p_dest ;
    double      dpix ;
    short       spix ;
    int         lpix ;
    float       fpix ;
    byte        XLpix[8] ;


    baseptr = p_dest = qfits_malloc(npix * sizeof(float));
    switch (bitpix) {

        case 8:
        /* No swapping needed */
        for (i=0 ; i<npix ; i++) {
            p_dest[i] = (float)((double)p_source[i] * bscale + bzero) ;
        }
        break ;

        case 16:
        for (i=0 ; i<npix ; i++) {
            memcpy(&spix, p_source, 2);
            p_source += 2 ;
#ifndef WORDS_BIGENDIAN
            spix = qfits_swap_bytes_16(spix);
#endif
            *p_dest++ = (float)(bscale * (double)spix + bzero) ;
        }
        break ;

        case 32:
        for (i=0 ; i<npix ; i++) {
            memcpy(&lpix, p_source, 4);
            p_source += 4 ;
#ifndef WORDS_BIGENDIAN
            lpix = qfits_swap_bytes_32(lpix);
#endif
            *p_dest++ = (float)(bscale * (double)lpix + bzero);
        }
        break ;

        case -32:
        for (i=0 ; i<npix ; i++) {
            memcpy(&lpix, p_source, 4);
            p_source += 4 ;
#ifndef WORDS_BIGENDIAN
            lpix = qfits_swap_bytes_32(lpix);
#endif
            memcpy(&fpix, &lpix, 4);
            *p_dest++ = (float)((double)fpix * bscale + bzero) ;
        }
        break ;

        case -64:
        for (i=0 ; i<npix ; i++) {
            XLpix[0] = *p_source ++ ;
            XLpix[1] = *p_source ++ ;
            XLpix[2] = *p_source ++ ;
            XLpix[3] = *p_source ++ ;
            XLpix[4] = *p_source ++ ;
            XLpix[5] = *p_source ++ ;
            XLpix[6] = *p_source ++ ;
            XLpix[7] = *p_source ++ ;
#ifndef WORDS_BIGENDIAN
            qfits_swap_bytes(XLpix, 8);
#endif
            dpix = *((double*)XLpix) ;
            *p_dest ++ = (float)(bscale * dpix + bzero);
        }
        break ;
    }
    return baseptr ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a pixel buffer as ints.
  @param    p_source    Pointer to source buffer (as bytes).
  @param    npix        Number of pixels to load.
  @param    bitpix      FITS BITPIX in original file.
  @param    bscale      FITS BSCALE in original file.
  @param    bzero       FITS BZERO in original file.
  @return   1 pointer to a newly allocated buffer of npix ints.

  This function takes in input a pointer to a byte buffer as given
  in the original FITS file (big-endian format). It converts the
  buffer to an array of int (whatever representation is used for
  int by this platform is used) and returns the newly allocated
  buffer, or NULL if an error occurred.

  The returned buffer must be deallocated using qfits_free().
 */
/*----------------------------------------------------------------------------*/
static int * qfits_pixin_int(
        byte    *   p_source,
        int         npix,
        int         bitpix,
        double      bscale,
        double      bzero)
{
    int         i ;
    int     *   p_dest ;
    int     *   baseptr ;
    double      dpix ;
    short       spix ;
    int         lpix ;
    float       fpix ;
    byte        XLpix[8] ;

    baseptr = p_dest = qfits_malloc(npix * sizeof(int));
    switch (bitpix) {

        case 8:
        /* No swapping needed */
        for (i=0 ; i<npix ; i++) {
            p_dest[i] = (int)((double)p_source[i] * bscale + bzero) ;
        }
        break ;

        case 16:
        for (i=0 ; i<npix ; i++) {
            memcpy(&spix, p_source, 2);
            p_source += 2 ;
#ifndef WORDS_BIGENDIAN
            spix = qfits_swap_bytes_16(spix);
#endif
            *p_dest++ = (int)(bscale * (double)spix + bzero) ;
        }
        break ;

        case 32:
        for (i=0 ; i<npix ; i++) {
            memcpy(&lpix, p_source, 4);
            p_source += 4 ;
#ifndef WORDS_BIGENDIAN
            lpix = qfits_swap_bytes_32(lpix);
#endif
            *p_dest++ = (int)(bscale * (double)lpix + bzero);
        }
        break ;

        case -32:
        for (i=0 ; i<npix ; i++) {
            memcpy(&lpix, p_source, 4);
            p_source += 4 ;
#ifndef WORDS_BIGENDIAN
            lpix = qfits_swap_bytes_32(lpix);
#endif
            memcpy(&fpix, &lpix, 4);
            *p_dest++ = (int)((double)fpix * bscale + bzero) ;
        }
        break ;

        case -64:
        for (i=0 ; i<npix ; i++) {
            XLpix[0] = *p_source ++ ;
            XLpix[1] = *p_source ++ ;
            XLpix[2] = *p_source ++ ;
            XLpix[3] = *p_source ++ ;
            XLpix[4] = *p_source ++ ;
            XLpix[5] = *p_source ++ ;
            XLpix[6] = *p_source ++ ;
            XLpix[7] = *p_source ++ ;
#ifndef WORDS_BIGENDIAN
            qfits_swap_bytes(XLpix, 8);
#endif
            dpix = *((double*)XLpix) ;
            *p_dest ++ = (int)(bscale * dpix + bzero);
        }
        break ;
    }
    return baseptr ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a pixel buffer as doubles.
  @param    p_source    Pointer to source buffer (as bytes).
  @param    npix        Number of pixels to load.
  @param    bitpix      FITS BITPIX in original file.
  @param    bscale      FITS BSCALE in original file.
  @param    bzero       FITS BZERO in original file.
  @return   1 pointer to a newly allocated buffer of npix doubles.

  This function takes in input a pointer to a byte buffer as given
  in the original FITS file (big-endian format). It converts the
  buffer to an array of double (whatever representation is used for
  int by this platform is used) and returns the newly allocated
  buffer, or NULL if an error occurred.

  The returned buffer must be deallocated using qfits_free().
 */
/*----------------------------------------------------------------------------*/
static double * qfits_pixin_double(
        byte    *   p_source,
        int         npix,
        int         bitpix,
        double      bscale,
        double      bzero)
{
    int         i ;
    double  *   p_dest ;
    double  *   baseptr ;
    double      dpix ;
    short       spix ;
    int         lpix ;
    float       fpix ;
    byte        XLpix[8] ;


    baseptr = p_dest = qfits_malloc(npix * sizeof(double));
    switch (bitpix) {

        case 8:
        /* No swapping needed */
        for (i=0 ; i<npix ; i++) {
            p_dest[i] = (double)((double)p_source[i] * bscale + bzero) ;
        }
        break ;

        case 16:
        for (i=0 ; i<npix ; i++) {
            memcpy(&spix, p_source, 2);
            p_source += 2 ;
#ifndef WORDS_BIGENDIAN
            spix = qfits_swap_bytes_16(spix);
#endif
            *p_dest++ = (double)(bscale * (double)spix + bzero) ;
        }
        break ;

        case 32:
        for (i=0 ; i<npix ; i++) {
            memcpy(&lpix, p_source, 4);
            p_source += 4 ;
#ifndef WORDS_BIGENDIAN
            lpix = qfits_swap_bytes_32(lpix);
#endif
            *p_dest++ = (double)(bscale * (double)lpix + bzero);
        }
        break ;

        case -32:
        for (i=0 ; i<npix ; i++) {
            memcpy(&lpix, p_source, 4);
            p_source += 4 ;
#ifndef WORDS_BIGENDIAN
            lpix = qfits_swap_bytes_32(lpix);
#endif
            memcpy(&fpix, &lpix, 4);
            *p_dest++ = (double)((double)fpix * bscale + bzero) ;
        }
        break ;

        case -64:
        for (i=0 ; i<npix ; i++) {
            XLpix[0] = *p_source ++ ;
            XLpix[1] = *p_source ++ ;
            XLpix[2] = *p_source ++ ;
            XLpix[3] = *p_source ++ ;
            XLpix[4] = *p_source ++ ;
            XLpix[5] = *p_source ++ ;
            XLpix[6] = *p_source ++ ;
            XLpix[7] = *p_source ++ ;
#ifndef WORDS_BIGENDIAN
            qfits_swap_bytes(XLpix, 8);
#endif
            dpix = *((double*)XLpix) ;
            *p_dest ++ = (double)(bscale * dpix + bzero);
        }
        break ;
    }
    return baseptr ;
}

/* Test code */
#ifdef TESTPIXIO
static void qfitsloader_dump(qfitsloader * ql)
{
    fprintf(stderr,
            "file      : %s\n"
            "xtnum     : %d\n"
            "pnum      : %d\n"
            "ptype     : %d\n"
            "lx        : %d\n"
            "ly        : %d\n"
            "np        : %d\n"
            "bitpix    : %d\n"
            "seg_start : %d\n"
            "bscale    : %g\n"
            "bzero     : %g\n"
            "ibuf      : %p\n"
            "fbuf      : %p\n"
            "dbuf      : %p\n",
            ql->filename,
            ql->xtnum,
            ql->pnum,
            ql->ptype,
            ql->lx,
            ql->ly,
            ql->np,
            ql->bitpix,
            ql->seg_start,
            ql->bscale,
            ql->bzero,
            ql->ibuf,
            ql->fbuf,
            ql->dbuf);
}

int main (int argc, char * argv[])
{
    qfitsloader ql ;

    if (argc<2) {
        printf("use: %s <FITS>\n", argv[0]);
        return 1 ;
    }

    ql.filename = argv[1] ;
    ql.xtnum    = 0 ;
    ql.pnum     = 0 ;
    ql.ptype    = PTYPE_FLOAT ;

    if (qfits_loadpix(&ql)!=0) {
        printf("error occurred during loading: abort\n");
        return -1 ;
    }
    qfitsloader_dump(&ql);
    printf("pix[0]=%g\n"
           "pix[100]=%g\n"
           "pix[10000]=%g\n",
           ql.fbuf[0],
           ql.fbuf[100],
           ql.fbuf[10000]);
    qfits_free(ql.fbuf);

    ql.ptype   = PTYPE_INT ;
    if (qfits_loadpix(&ql)!=0) {
        printf("error occurred during loading: abort\n");
        return -1 ;
    }
    qfitsloader_dump(&ql);
    printf("pix[0]=%d\n"
           "pix[100]=%d\n"
           "pix[10000]=%d\n",
           ql.ibuf[0],
           ql.ibuf[100],
           ql.ibuf[10000]);
    qfits_free(ql.ibuf);


    ql.ptype   = PTYPE_DOUBLE ;
    if (qfits_loadpix(&ql)!=0) {
        printf("error occurred during loading: abort\n");
        return -1 ;
    }
    qfitsloader_dump(&ql);
    printf("pix[0]=%g\n"
           "pix[100]=%g\n"
           "pix[10000]=%g\n",
           ql.dbuf[0],
           ql.dbuf[100],
           ql.dbuf[10000]);
    qfits_free(ql.dbuf);

    return 0 ;
}
#endif
#endif
#ifndef QFITS_RW_H
#define QFITS_RW_H
/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/








/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

/* FITS magic number */
#define FITS_MAGIC            "SIMPLE"
/* Size of the FITS magic number */
#define FITS_MAGIC_SZ        6

/*-----------------------------------------------------------------------------
                        Private to this module
 -----------------------------------------------------------------------------*/

static int is_blank_line(const char *) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_rw    FITS header reading/writing
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                            Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Read a FITS header from a file to an internal structure.
  @param    filename    Name of the file to be read
  @return   Pointer to newly allocated qfits_header or NULL in error case.

  This function parses a FITS (main) header, and returns an allocated
  qfits_header object. The qfits_header object contains a linked-list of
  key "tuples". A key tuple contains:

  - A keyword
  - A value
  - A comment
  - An original FITS line (as read from the input file)

  Direct access to the structure is not foreseen, use accessor
  functions in fits_h.h

  Value, comment, and original line might be NULL pointers.
 */
/*----------------------------------------------------------------------------*/
qfits_header * qfits_header_read(const char * filename)
{
    /* Forward job to readext */
    return qfits_header_readext(filename, 0);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Read a FITS header from a 'hdr' file.
  @param    filename    Name of the file to be read
  @return   Pointer to newly allocated qfits_header or NULL in error case

  This function parses a 'hdr' file, and returns an allocated qfits_header
  object. A hdr file is an ASCII format were the header is written with a
  carriage return after each line. The command dfits typically displays
  a hdr file.
 */
/*----------------------------------------------------------------------------*/
qfits_header * qfits_header_read_hdr(const char * filename)
{
    qfits_header    *   hdr ;
    FILE            *   in ;
    char                line[81];
    char            *   key,
                    *   val,
                    *   com ;
    int                 i, j ;

    /* Check input */
    if (filename==NULL) return NULL ;

    /* Initialise */
    key = val = com = NULL ;

    /* Open the file */
    if ((in=fopen(filename, "r"))==NULL) {
        qfits_error("cannot read [%s]", filename) ;
        return NULL ;
    }

    /* Create the header */
    hdr = qfits_header_new() ;

    /* Go through the file */
    while (fgets(line, 81, in)!=NULL) {
        for (i=0 ; i<81 ; i++) {
            if (line[i] == '\n') {
                for (j=i ; j<81 ; j++) line[j] = ' ' ;
                line[80] = (char)0 ;
                break ;
            }
        }
        if (!strcmp(line, "END")) {
            line[3] = ' ';
            line[4] = (char)0 ;
        }

        /* Rule out blank lines */
        if (!is_blank_line(line)) {

            /* Get key, value, comment for the current line */
            key = qfits_getkey(line);
            val = qfits_getvalue(line);
            com = qfits_getcomment(line);

            /* If key or value cannot be found, trigger an error */
            if (key==NULL) {
                qfits_header_destroy(hdr);
                fclose(in) ;
                return NULL ;
            }
            /* Append card to linked-list */
            qfits_header_append(hdr, key, val, com, NULL);
        }
    }
    fclose(in) ;

    /* The last key should be 'END' */
    if (strlen(key)!=3) {
        qfits_header_destroy(hdr);
        return NULL ;
    }
    if (key[0]!='E' || key[1]!='N' || key[2]!='D') {
        qfits_header_destroy(hdr);
        return NULL ;
    }

    return hdr ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Read a FITS header from a 'hdr' string
  @param    hdr_str         String containing the hdr file
  @param    nb_char         Number of characters in the string
  @return   Pointer to newly allocated qfits_header or NULL in error case

  This function parses a 'hdr' string, and returns an allocated qfits_header
  object.
 */
/*----------------------------------------------------------------------------*/
qfits_header * qfits_header_read_hdr_string(
        const unsigned char *   hdr_str,
        int                     nb_char)
{
    qfits_header    *   hdr ;
    char                line[81];
    char            *   key,
                    *   val,
                    *   com ;
    int                 ind ;
    int                 i, j ;

    /* Check input */
    if (hdr_str==NULL) return NULL ;

    /* Initialise */
    key = val = com = NULL ;

    /* Create the header */
    hdr = qfits_header_new() ;

    /* Go through the file */
    ind = 0 ;
    while (ind <= nb_char - 80) {
        strncpy(line, (char*)hdr_str + ind, 80) ;
        line[80] = (char)0 ;
        for (i=0 ; i<81 ; i++) {
            if (line[i] == '\n') {
                for (j=i ; j<81 ; j++) line[j] = ' ' ;
                line[80] = (char)0 ;
                break ;
            }
        }
        if (!strcmp(line, "END")) {
            line[3] = ' ';
            line[4] = (char)0 ;
        }

        /* Rule out blank lines */
        if (!is_blank_line(line)) {

            /* Get key, value, comment for the current line */
            key = qfits_getkey(line);
            val = qfits_getvalue(line);
            com = qfits_getcomment(line);

            /* If key or value cannot be found, trigger an error */
            if (key==NULL) {
                qfits_header_destroy(hdr);
                return NULL ;
            }
            /* Append card to linked-list */
            qfits_header_append(hdr, key, val, com, NULL);
        }
        ind += 80 ;
    }

    /* The last key should be 'END' */
    if (strlen(key)!=3) {
        qfits_header_destroy(hdr);
        return NULL ;
    }
    if (key[0]!='E' || key[1]!='N' || key[2]!='D') {
        qfits_header_destroy(hdr);
        return NULL ;
    }

    return hdr ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Read an extension header from a FITS file.
  @param    filename    Name of the FITS file to read
  @param    xtnum       Extension number to read, starting from 0.
  @return   Newly allocated qfits_header structure.

  Strictly similar to qfits_header_read() but reads headers from
  extensions instead. If the requested xtension is 0, this function
  returns the main header.

  Returns NULL in case of error.
 */
/*----------------------------------------------------------------------------*/
qfits_header * qfits_header_readext(const char * filename, int xtnum)
{
    qfits_header*   hdr ;
    int             n_ext ;
    char            line[81];
    char        *   where ;
    char        *   start ;
    char        *   key,
                *   val,
                *   com ;
    int             seg_start ;
    int             seg_size ;
    size_t          size ;

    /* Check input */
    if (filename==NULL || xtnum<0)
        return NULL ;

    /* Check that there are enough extensions */
    if (xtnum>0) {
        n_ext = qfits_query_n_ext(filename);
        if (xtnum>n_ext) {
            return NULL ;
        }
    }

    /* Get offset to the extension header */
    if (qfits_get_hdrinfo(filename, xtnum, &seg_start, &seg_size)!=0) {
        return NULL ;
    }

    /* Memory-map the input file */
    start = qfits_falloc((char *)filename, seg_start, &size) ;
    if (start==NULL) return NULL ;

    hdr   = qfits_header_new() ;
    where = start ;
    while (1) {
        memcpy(line, where, 80);
        line[80] = (char)0;

        /* Rule out blank lines */
        if (!is_blank_line(line)) {

            /* Get key, value, comment for the current line */
            key = qfits_getkey(line);
            val = qfits_getvalue(line);
            com = qfits_getcomment(line);

            /* If key or value cannot be found, trigger an error */
            if (key==NULL) {
                qfits_header_destroy(hdr);
                hdr = NULL ;
                break ;
            }
            /* Append card to linked-list */
            qfits_header_append(hdr, key, val, com, line);
            /* Check for END keyword */
            if (strlen(key)==3)
                if (key[0]=='E' &&
                    key[1]=='N' &&
                    key[2]=='D')
                    break ;
        }
        where += 80 ;
        /* If reaching the end of file, trigger an error */
        if ((int)(where-start)>=(int)(seg_size+80)) {
            qfits_header_destroy(hdr);
            hdr = NULL ;
            break ;
        }
    }
    qfits_fdealloc(start, seg_start, size) ;
    return hdr ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Pad an existing file with zeros to a multiple of 2880.
  @param    filename    Name of the file to pad.
  @return   void

  This function simply pads an existing file on disk with enough zeros
  for the file size to reach a multiple of 2880, as required by FITS.
 */
/*----------------------------------------------------------------------------*/
void qfits_zeropad(const char * filename)
{
    struct stat sta ;
    int         size ;
    int         remaining;
    FILE    *   out ;
    char    *   buf;

    if (filename==NULL) return ;

    /* Get file size in bytes */
    if (stat(filename, &sta)!=0) {
        return ;
    }
    size = (int)sta.st_size ;
    /* Compute number of zeros to pad */
    remaining = size % FITS_BLOCK_SIZE ;
    if (remaining==0) return ;
    remaining = FITS_BLOCK_SIZE - remaining ;

    /* Open file, dump zeros, exit */
    if ((out=fopen(filename, "a"))==NULL)
        return ;
    buf = qfits_calloc(remaining, sizeof(char));
    fwrite(buf, 1, remaining, out);
    fclose(out);
    qfits_free(buf);
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Identify if a file is a FITS file.
  @param    filename name of the file to check
  @return   int 0, 1, or -1

  Returns 1 if the file name looks like a valid FITS file. Returns
  0 else. If the file does not exist, returns -1.
 */
/*----------------------------------------------------------------------------*/
int qfits_is_fits(const char * filename)
{
    FILE  *   fp ;
    char  *   magic ;
    int       isfits ;

    if (filename==NULL) return -1 ;
    if ((fp = fopen(filename, "r"))==NULL) {
        qfits_error("cannot open file [%s]", filename) ;
        return -1 ;
    }

    magic = qfits_calloc(FITS_MAGIC_SZ+1, sizeof(char)) ;
    fread(magic, 1, FITS_MAGIC_SZ, fp) ;
    fclose(fp) ;
    magic[FITS_MAGIC_SZ] = (char)0 ;
    if (strstr(magic, FITS_MAGIC)!=NULL)
        isfits = 1 ;
    else
        isfits = 0 ;
    qfits_free(magic) ;
    return isfits ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Retrieve offset to start and size of a header in a FITS file.
  @param    filename    Name of the file to examine
  @param    xtnum       Extension number (0 for main)
  @param    seg_start   Segment start in bytes (output)
  @param    seg_size    Segment size in bytes (output)
  @return   int 0 if Ok, -1 otherwise.

  This function retrieves the two most important informations about
  a header in a FITS file: the offset to its beginning, and the size
  of the header in bytes. Both values are returned in the passed
  pointers to ints. It is Ok to pass NULL for any pointer if you do
  not want to retrieve the associated value.

  You must provide an extension number for the header, 0 meaning the
  main header in the file.
 */
/*----------------------------------------------------------------------------*/
int qfits_get_hdrinfo(
        const char  *   filename,
        int             xtnum,
        int         *   seg_start,
        int         *   seg_size)
{
    if (filename==NULL || xtnum<0 || (seg_start==NULL && seg_size==NULL)) {
        return -1 ;
    }
    if (seg_start!=NULL) {
        *seg_start = qfits_query(filename, QFITS_QUERY_HDR_START | xtnum);
        if (*seg_start<0)
            return -1 ;
    }
    if (seg_size!=NULL) {
        *seg_size = qfits_query(filename, QFITS_QUERY_HDR_SIZE | xtnum);
        if (*seg_size<0)
            return -1 ;
    }
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Retrieve offset to start and size of a data section in a file.
  @param    filename    Name of the file to examine.
  @param    xtnum       Extension number (0 for main).
  @param    seg_start   Segment start in bytes (output).
  @param    seg_size    Segment size in bytes (output).
  @return   int 0 if Ok, -1 otherwise.

  This function retrieves the two most important informations about
  a data section in a FITS file: the offset to its beginning, and the size
  of the section in bytes. Both values are returned in the passed
  pointers to ints. It is Ok to pass NULL for any pointer if you do
  not want to retrieve the associated value.

  You must provide an extension number for the header, 0 meaning the
  main header in the file.
 */
/*----------------------------------------------------------------------------*/
int qfits_get_datinfo(
        const char  *   filename,
        int             xtnum,
        int         *   seg_start,
        int         *   seg_size)
{
    if (filename==NULL || xtnum<0 || (seg_start==NULL && seg_size==NULL)) {
        return -1 ;
    }
    if (seg_start!=NULL) {
        *seg_start = qfits_query(filename, QFITS_QUERY_DAT_START | xtnum);
        if (*seg_start<0)
            return -1 ;
    }
    if (seg_size!=NULL) {
        *seg_size = qfits_query(filename, QFITS_QUERY_DAT_SIZE | xtnum);
        if (*seg_size<0)
            return -1 ;
    }
    return 0 ;
}

/**@}*/

static int is_blank_line(const char * s)
{
    int     i ;

    for (i=0 ; i<(int)strlen(s) ; i++) {
        if (s[i]!=' ') return 0 ;
    }
    return 1 ;
}
#endif
#ifndef QFITS_TIME_H
#define QFITS_TIME_H

/*-----------------------------------------------------------------------------
                                   Macros
 -----------------------------------------------------------------------------*/

/* Get century from a date in long format */
#define GET_CENTURY(d)      (int) ( (d) / 1000000L)
/* Get century year from a date in long format */
#define GET_CCYEAR(d)       (int) ( (d) / 10000L)
/* Get year from a date in long format */
#define GET_YEAR(d)         (int) (((d) % 1000000L) / 10000L)
/* Get month from a date in long format */
#define GET_MONTH(d)        (int) (((d) % 10000L) / 100)
/* Get day from a date in long format */
#define GET_DAY(d)          (int) ( (d) % 100)

/* Get hours from a date in long format */
#define GET_HOUR(t)         (int) ( (t) / 1000000L)
/* Get minutes from a date in long format */
#define GET_MINUTE(t)       (int) (((t) % 1000000L) / 10000L)
/* Get seconds from a date in long format */
#define GET_SECOND(t)       (int) (((t) % 10000L) / 100)
/* Get centi-seconds from a date in long format */
#define GET_CENTI(t)        (int) ( (t) % 100)

/* Make date in long format from its components */
#define MAKE_DATE(c,y,m,d)  (long) (c) * 1000000L +                          \
                            (long) (y) * 10000L +                            \
                            (long) (m) * 100 + (d)
/* Make time in long format from its components */
#define MAKE_TIME(h,m,s,c)  (long) (h) * 1000000L +                          \
                            (long) (m) * 10000L +                            \
                            (long) (s) * 100 + (c)

/*  Interval values, specified in centiseconds */
#define INTERVAL_CENTI      1
#define INTERVAL_SEC        100
#define INTERVAL_MIN        6000
#define INTERVAL_HOUR       360000L
#define INTERVAL_DAY        8640000L

/*-----------------------------------------------------------------------------
                            Private to this module
 -----------------------------------------------------------------------------*/

static long timer_to_date(time_t time_secs) ;
static long timer_to_time(time_t time_secs) ;
static long qfits_time_now(void) ;
static long qfits_date_now (void) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_time  Get date/time, possibly in ISO8601 format
 *
 * This module contains various utilities to get the current date/time,
 * and possibly format it according to the ISO 8601 format.
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Returns the current date and time as a static string.
  @return   Pointer to statically allocated string

  Build and return a string containing the date of today and the
  current time in ISO8601 format. The returned pointer points to a
  statically allocated string in the function, so no need to free it.
 */
/*----------------------------------------------------------------------------*/
char * qfits_get_datetime_iso8601(void)
{
    static char date_iso8601[20] ;
    long        curdate ;
    long        curtime ;

    curdate  = qfits_date_now() ;
    curtime  = qfits_time_now() ;

    sprintf(date_iso8601, "%04d-%02d-%02dT%02d:%02d:%02d",
            GET_CCYEAR(curdate),
            GET_MONTH(curdate),
            GET_DAY(curdate),
            GET_HOUR(curtime),
            GET_MINUTE(curtime),
            GET_SECOND(curtime));
    return date_iso8601 ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Returns the current date as a long (CCYYMMDD).
  @return    The current date as a long number.

  Returns the current date as a long value (CCYYMMDD). Since most
  system clocks do not return a century, this function assumes that
  all years 80 and above are in the 20th century, and all years 00 to
  79 are in the 21st century.  For best results, consume before 1 Jan
  2080.
  Example:    19 Oct 2000 is returned as 20001019
 */
/*----------------------------------------------------------------------------*/
static long qfits_date_now (void)
{
    return (timer_to_date (time (NULL)));
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Returns the current time as a long (HHMMSSCC).
  @return    The current time as a long number.

  Returns the current time as a long value (HHMMSSCC). If the system
  clock does not return centiseconds, these are set to zero.

  Example: 15:36:12.84 is returned as 15361284
 */
/*----------------------------------------------------------------------------*/
static long qfits_time_now(void)
{
    struct timeval time_struct;

    gettimeofday (&time_struct, 0);
    return (timer_to_time (time_struct.tv_sec)
                         + time_struct.tv_usec / 10000);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Converts a timer value to a date.
  @param    time_secs    Current time definition in seconds.
  @return    Current date as a long (CCYYMMDD).

  Converts the supplied timer value into a long date value. Dates are
  stored as long values: CCYYMMDD. If the supplied value is zero,
  returns zero.  If the supplied value is out of range, returns 1
  January, 1970 (19700101). The timer value is assumed to be UTC
  (GMT).
 */
/*----------------------------------------------------------------------------*/
static long timer_to_date(time_t time_secs)
{
    struct tm *time_struct;

    if (time_secs == 0) {
        return 0;
    } else {
        /*  Convert into a long value CCYYMMDD */
        time_struct = localtime (&time_secs);
        if (time_struct) {
            time_struct-> tm_year += 1900;
            return (MAKE_DATE (    time_struct-> tm_year / 100,
                                time_struct-> tm_year % 100,
                                time_struct-> tm_mon + 1,
                                time_struct-> tm_mday));
        } else {
            return (19700101);
        }
    }
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert a timer value to a time.
  @param    time_secs    Current time definition in seconds.
  @return    Current time as a long.

  Converts the supplied timer value into a long time value.  Times are
  stored as long values: HHMMSS00.  Since the timer value does not
  hold centiseconds, these are set to zero.  If the supplied value was
  zero or invalid, returns zero.  The timer value is assumed to be UTC
  (GMT).
 */
/*----------------------------------------------------------------------------*/
static long timer_to_time(time_t time_secs)
{
    struct tm *time_struct;

    if (time_secs == 0) {
        return 0;
    } else {
        /*  Convert into a long value HHMMSS00 */
        time_struct = localtime (&time_secs);
        if (time_struct) {
            return (MAKE_TIME (time_struct-> tm_hour,
                               time_struct-> tm_min,
                               time_struct-> tm_sec,
                               0));
        } else {
            return 0;
        }
    }
}
#endif
#ifndef QFITS_TABLE_H
#define QFITS_TABLE_H


/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/











/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

#define ELEMENT_MAX_DISPLAY_SIZE    50


/*----------------------------------------------------------------------------*/
/**
 * @defgroup    qfits_table FITS table handling
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Identify a file as containing a FITS table in extension.
  @param    filename    Name of the FITS file to examine.
  @param    xtnum        Extension number to check (starting from 1).
  @return    int 1 if the extension contains a table, 0 else.
  Examines the requested extension and identifies the presence of a FITS table.
 */
/*----------------------------------------------------------------------------*/
int qfits_is_table(const char * filename, int xtnum)
{
    char    *    value ;
    int            ttype ;

    ttype = QFITS_INVALIDTABLE ;
    value = qfits_query_ext(filename, "XTENSION", xtnum);
    if (value==NULL) return ttype ;

    value = qfits_pretty_string(value);
    if (!strcmp(value, "TABLE")) {
        ttype = QFITS_ASCIITABLE;
    } else if (!strcmp(value, "BINTABLE")) {
        ttype = QFITS_BINTABLE;
    }
    return ttype;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Generate a default primary header to store tables
  @return    the header object
 */
/*----------------------------------------------------------------------------*/
qfits_header * qfits_table_prim_header_default(void)
{
    qfits_header    *    fh ;

    fh = qfits_header_new() ;

    qfits_header_append(fh, "SIMPLE", "T", "Standard FITS file", NULL) ;
    qfits_header_append(fh, "BITPIX", "8", "ASCII or bytes array", NULL) ;
    qfits_header_append(fh, "NAXIS", "0", "Minimal header", NULL) ;
    qfits_header_append(fh, "EXTEND", "T", "There may be FITS ext", NULL);
    qfits_header_append(fh, "END", NULL, NULL, NULL) ;

    return fh ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Generate a default extension header to store tables
  @return   the header object
 */
/*----------------------------------------------------------------------------*/
qfits_header * qfits_table_ext_header_default(const qfits_table * t)
{
    qfits_header    *   fh ;
    qfits_col       *   curr_col ;
    char                str_val[FITS_LINESZ] ;
    char                str_val2[FITS_LINESZ] ;
    char            *   date ;
    int                 tab_width ;
    int                 col_pos ;
    int                 i ;

    /* Compute the table width   */
    if ((tab_width = qfits_compute_table_width(t)) == -1) {
        qfits_error("cannot get the table width") ;
        return NULL ;
    }

    /* Create fits header */
    if ((fh=qfits_header_new()) == NULL) {
        qfits_error("cannot create new fits header") ;
        return NULL ;
    }

    /* Check the kind of table */
    if (t->tab_t == QFITS_BINTABLE) {

        /* Write extension header */
        qfits_header_append(fh, "XTENSION", "BINTABLE",
                "FITS Binary Table Extension", NULL) ;
        qfits_header_append(fh, "BITPIX", "8", "8-bits character format", NULL);
        qfits_header_append(fh, "NAXIS", "2","Tables are 2-D char. array",NULL);
        sprintf(str_val, "%d", tab_width) ;
        qfits_header_append(fh, "NAXIS1", str_val, "Bytes in row", NULL) ;
        sprintf(str_val, "%d", (int)(t->nr)) ;
        qfits_header_append(fh, "NAXIS2", str_val, "No. of rows in table",NULL);
        qfits_header_append(fh, "PCOUNT", "0", "Parameter count always 0",NULL);
        qfits_header_append(fh, "GCOUNT", "1", "Group count always 1", NULL);
        sprintf(str_val, "%d", (int)(t->nc)) ;
        qfits_header_append(fh, "TFIELDS", str_val, "No. of col in table",NULL);
        /* Columns descriptors */
        curr_col = t->col ;
        for (i=0 ; i<t->nc ; i++) {
            sprintf(str_val, "TFORM%d", i+1) ;
            sprintf(str_val2, "'%s'", qfits_build_format(curr_col)) ;
            qfits_header_append(fh, str_val, str_val2, "Format of field", NULL);

            sprintf(str_val, "TTYPE%d", i+1) ;
            sprintf(str_val2, "%s", curr_col->tlabel) ;
            qfits_header_append(fh, str_val, str_val2, "Field label", NULL) ;

            sprintf(str_val, "TUNIT%d", i+1) ;
            sprintf(str_val2, "%s", curr_col->tunit) ;
            qfits_header_append(fh, str_val, str_val2, "Physical unit of field",
                    NULL) ;
            if (curr_col->zero_present) {
                sprintf(str_val, "TZERO%d", i+1) ;
                sprintf(str_val2, "%f", curr_col->zero) ;
                qfits_header_append(fh, str_val, str_val2,
                        "NULL value is defined", NULL) ;
            }
            if (curr_col->scale_present) {
                sprintf(str_val, "TSCAL%d", i+1) ;
                sprintf(str_val2, "%f", curr_col->scale) ;
                qfits_header_append(fh, str_val, str_val2, "Scaling applied",
                        NULL);
            }
            curr_col++ ;
        }
        qfits_header_append(fh,"ORIGIN","ESO-QFITS", "Written by QFITS", NULL);

        date = qfits_get_datetime_iso8601() ;
        sprintf(str_val, "'%s'", date) ;
        qfits_header_append(fh, "DATE", str_val, "[UTC] Date of writing", NULL);
        qfits_header_append(fh, "END", NULL, NULL, NULL);

    } else if (t->tab_t == QFITS_ASCIITABLE) {

        /* Write extension header */
        qfits_header_append(fh, "XTENSION", "TABLE",
                        "FITS ASCII Table Extension", NULL) ;
        qfits_header_append(fh, "BITPIX", "8", "8-bits character format", NULL);
        qfits_header_append(fh, "NAXIS", "2", "ASCII table has 2 axes", NULL) ;

        /* Fill the header  */
        sprintf(str_val, "%d", tab_width) ;
        qfits_header_append(fh, "NAXIS1", str_val, "Characters in a row", NULL);
        sprintf(str_val, "%d", (int)(t->nr)) ;
        qfits_header_append(fh, "NAXIS2", str_val, "No. of rows in table",NULL);
        qfits_header_append(fh, "PCOUNT", "0", "No group parameters", NULL) ;
        qfits_header_append(fh, "GCOUNT", "1", "Only one group", NULL);
        sprintf(str_val, "%d", (int)(t->nc)) ;
        qfits_header_append(fh, "TFIELDS", str_val, "No. of col in table",NULL);
        qfits_header_append(fh, "ORIGIN","ESO-QFITS","Written by QFITS",NULL);
        date = qfits_get_datetime_iso8601() ;
        sprintf(str_val, "'%s'", date) ;
        qfits_header_append(fh, "DATE", str_val, "[UTC] Date of writing", NULL);

        /* Columns descriptors */
        curr_col = t->col ;
        col_pos = 1 ;
        for (i=0 ; i<t->nc ; i++) {
            sprintf(str_val, "TTYPE%d", i+1) ;
            sprintf(str_val2, "%s", curr_col->tlabel) ;
            qfits_header_append(fh, str_val, str_val2, "Field label", NULL) ;

            sprintf(str_val, "TFORM%d", i+1) ;
            sprintf(str_val2, "'%s'", qfits_build_format(curr_col)) ;
            qfits_header_append(fh, str_val, str_val2, "Format of field", NULL);

            sprintf(str_val, "TBCOL%d", i+1) ;
            sprintf(str_val2, "%d", col_pos) ;
            qfits_header_append(fh, str_val, str_val2,"Start column of field",
                    NULL);
            col_pos += curr_col->atom_nb ;

            sprintf(str_val, "TUNIT%d", i+1) ;
            sprintf(str_val2, "%s", curr_col->tunit) ;
            qfits_header_append(fh, str_val, str_val2, "Physical unit of field",
                    NULL) ;
            if (curr_col->zero_present) {
                sprintf(str_val, "TZERO%d", i+1) ;
                sprintf(str_val2, "%f", curr_col->zero) ;
                qfits_header_append(fh, str_val, str_val2,
                        "NULL value is defined", NULL) ;
            }
            if (curr_col->scale_present) {
                sprintf(str_val, "TSCAL%d", i+1) ;
                sprintf(str_val2, "%f", curr_col->scale) ;
                qfits_header_append(fh, str_val, str_val2, "Scaling applied",
                        NULL);
            }
            curr_col++ ;
        }
        qfits_header_append(fh, "END", NULL, NULL, NULL);

    } else {
        qfits_error("Table type not known") ;
        qfits_header_destroy(fh) ;
        return NULL ;
    }
    return fh ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Table object constructor
  @param    filename    Name of the FITS file associated to the table
  @param    table_type    Type of the table (QFITS_ASCIITABLE or QFITS_BINTABLE)
  @param    table_width Width in bytes of the table
  @param    nb_cols        Number of columns
  @param    nb_raws        Number of raws
  @return    The table object
  The columns are also allocated. The object has to be deallocated with
  qfits_table_close()
 */
/*----------------------------------------------------------------------------*/
qfits_table * qfits_table_new(
        const char  *   filename,
        int             table_type,
        int             table_width,
        int             nb_cols,
        int             nb_raws)
{
    qfits_table    *    qt ;
    qt = qfits_malloc(sizeof(qfits_table)) ;
    (void)strcpy(qt->filename, filename) ;
    qt->tab_t = table_type ;
    qt->nc = nb_cols ;
    qt->nr = nb_raws ;
    qt->col = qfits_calloc(qt->nc, sizeof(qfits_col)) ;
    qt->tab_w = table_width ;

    return qt ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Fill a column object with some provided informations
  @param    qc      Pointer to the column that has to be filled
  @param    unit    Unit of the data
  @param    label   Label of the column
  @param    disp    Way to display the data
  @param    nullval Null value
  @param    atom_nb Number of atoms per field. According to the type, an atom
                      is a double, an int, a char, ...
  @param    atom_dec_nb Number of decimals as specified in TFORM
  @param    atom_size    Size in bytes of the field for ASCII tables, and of
                          an atom for BIN tables. ASCII tables only contain 1
                        atom per field (except for A type where you can of
                        course have more than one char per field)
  @param    atom_type    Type of data (11 types for BIN, 5 for ASCII)
  @param    zero_present    Flag to use or not zero
  @param    zero            Zero value
  @param    scale_present   Flag to use or not scale
  @param    scale           Scale value
  @param    offset_beg  Gives the position of the column
  @return     -1 in error case, 0 otherwise
 */
/*----------------------------------------------------------------------------*/
int qfits_col_fill(
        qfits_col   *   qc,
        int             atom_nb,
        int             atom_dec_nb,
        int             atom_size,
        tfits_type      atom_type,
        const char  *   label,
        const char  *   unit,
        const char  *   nullval,
        const char  *   disp,
        int             zero_present,
        float           zero,
        int             scale_present,
        float           scale,
        int             offset_beg)
{
    /* Number of atoms per column */
    qc->atom_nb = atom_nb ;

    /* Number of decimals in a field in ASCII table (0 in BINTABLE) */
    qc->atom_dec_nb = atom_dec_nb ;

    /* Size in bytes of an atom  */
    qc->atom_size = atom_size ;

    /* Data type in the column */
    qc->atom_type = atom_type ;

    /* Label of the column */
    (void)strcpy(qc->tlabel, label) ;

    /* Unit of the column data */
    (void)strcpy(qc->tunit, unit) ;

    /* Null value*/
    (void)strcpy(qc->nullval, nullval) ;

    /* How to display the data */
    (void)strcpy(qc->tdisp, disp) ;

    /* Default values for zero and scales */
    qc->zero_present = zero_present ;
    qc->scale_present = scale_present ;
    qc->zero = zero ;
    qc->scale = scale ;

    /* Number of bytes between two consecutive fields of the same column */
    qc->off_beg = offset_beg ;

    /* A column is a priori readable */
    qc->readable = 1 ;

    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Read a FITS extension.
  @param    filename    Name of the FITS file to examine.
  @param    xtnum        Extension number to read (starting from 1).
  @return    Pointer to newly allocated qfits_table structure.

  Read a FITS table from a given file name and extension, and return a
  newly allocated qfits_table structure.
 */
/*----------------------------------------------------------------------------*/
qfits_table * qfits_table_open(
        const char  *   filename,
        int             xtnum)
{
    qfits_table     *   tload ;
    qfits_col       *   curr_col ;
    char            *   str_val ;
    char                keyword[FITSVALSZ] ;
    /* Table infos  */
    int                 table_type ;
    int                 nb_col ;
    int                 table_width ;
    int                 nb_raws ;
    /* Column infos */
    char                label[FITSVALSZ] ;
    char                unit[FITSVALSZ] ;
    char                disp[FITSVALSZ] ;
    char                nullval[FITSVALSZ] ;
    int                 atom_nb ;
    int                 atom_dec_nb ;
    int                 atom_size ;
    tfits_type          atom_type ;
    int                 offset_beg ;
    int                 data_size ;
    int                 theory_size ;
    int                 zero_present ;
    int                 scale_present ;
    float               zero ;
    float               scale ;

    /* For ASCII tables */
    int                    col_pos ;
    int                    next_col_pos ;

    /* For X type */
    int                    nb_bits ;

    int                    i ;

     /* See if 'filename' is a fits file  */
    if (qfits_is_fits(filename) != 1) {
        qfits_error("[%s] is not FITS", filename) ;
        return NULL ;
    }

    /* Identify a table and get the table type : ASCII or BIN */
    if ((table_type = qfits_is_table(filename, xtnum))==0) {
        qfits_error("[%s] extension %d is not a table", filename, xtnum) ;
        return NULL ;
    }

    /* Get number of columns and allocate them: nc <-> TFIELDS */
    if ((str_val = qfits_query_ext(filename, "TFIELDS", xtnum)) == NULL) {
        qfits_error("cannot read TFIELDS in [%s]:[%d]", filename, xtnum) ;
        return NULL ;
    }
    nb_col = atoi(str_val) ;

    /* Get the width in bytes of the table */
    if ((str_val = qfits_query_ext(filename, "NAXIS1", xtnum)) == NULL) {
        qfits_error("cannot read NAXIS1 in [%s]:[%d]", filename, xtnum) ;
        return NULL ;
    }
    table_width = atoi(str_val) ;

    /* Get the number of raws */
    if ((str_val = qfits_query_ext(filename, "NAXIS2", xtnum)) == NULL) {
        qfits_error("cannot read NAXIS2 in [%s]:[%d]", filename, xtnum) ;
        return NULL ;
    }
    nb_raws = atoi(str_val) ;

    /* Create the table object */
    tload = qfits_table_new(filename, table_type, table_width, nb_col, nb_raws);

    /* Initialize offset_beg */
    if (qfits_get_datinfo(filename, xtnum, &offset_beg, &data_size)!=0) {
        qfits_error("cannot find data start in [%s]:[%d]", filename, xtnum);
        qfits_table_close(tload);
        return NULL ;
    }

    /* Loop on all columns and get column descriptions  */
    curr_col = tload->col ;
    for (i=0 ; i<tload->nc ; i++) {
        /* label <-> TTYPE     */
        sprintf(keyword, "TTYPE%d", i+1) ;
        if ((str_val=qfits_query_ext(filename, keyword, xtnum)) == NULL) {
            label[0] = (char)0 ;
        } else strcpy(label, qfits_pretty_string(str_val)) ;

        /* unit <-> TUNIT */
        sprintf(keyword, "TUNIT%d", i+1) ;
        if ((str_val=qfits_query_ext(filename, keyword, xtnum)) == NULL) {
            unit[0] = (char)0 ;
        } else strcpy(unit, qfits_pretty_string(str_val)) ;

        /* disp <-> TDISP */
        sprintf(keyword, "TDISP%d", i+1) ;
        if ((str_val=qfits_query_ext(filename, keyword, xtnum)) == NULL) {
            disp[0] = (char)0 ;
        } else strcpy(disp, qfits_pretty_string(str_val)) ;

        /* nullval <-> TNULL */
        sprintf(keyword, "TNULL%d", i+1) ;
        if ((str_val=qfits_query_ext(filename, keyword, xtnum)) == NULL) {
            nullval[0] = (char)0 ;
        } else strcpy(nullval, qfits_pretty_string(str_val)) ;

        /* atom_size, atom_nb, atom_dec_nb, atom_type    <-> TFORM */
        sprintf(keyword, "TFORM%d", i+1) ;
        if ((str_val=qfits_query_ext(filename, keyword, xtnum))==NULL) {
            qfits_error("cannot read [%s] in [%s]:[%d]", keyword, filename,
                    xtnum);
            qfits_table_close(tload);
            return NULL ;
        }
        /* Interpret the type in header */
        if (qfits_table_interpret_type(qfits_pretty_string(str_val),
                        &(atom_nb),
                        &(atom_dec_nb),
                        &(atom_type),
                        table_type) == -1) {
            qfits_error("cannot interpret the type: %s", str_val) ;
            qfits_table_close(tload) ;
            return NULL ;
        }

        /* Set atom_size */
        switch (atom_type) {
            case TFITS_BIN_TYPE_A:
            case TFITS_BIN_TYPE_L:
            case TFITS_BIN_TYPE_B:
                atom_size = 1 ;
                break ;
            case TFITS_BIN_TYPE_I:
                atom_size = 2 ;
                break ;
            case TFITS_BIN_TYPE_J:
            case TFITS_BIN_TYPE_E:
            case TFITS_ASCII_TYPE_I:
            case TFITS_ASCII_TYPE_E:
            case TFITS_ASCII_TYPE_F:
                atom_size = 4 ;
                break ;
            case TFITS_BIN_TYPE_C:
            case TFITS_BIN_TYPE_P:
                atom_size = 4 ;
                atom_nb *= 2 ;
                break ;
            case TFITS_BIN_TYPE_D:
            case TFITS_ASCII_TYPE_D:
                atom_size = 8 ;
                break ;
            case TFITS_BIN_TYPE_M:
                atom_size = 8 ;
                atom_nb *= 2 ;
                break ;
            case TFITS_BIN_TYPE_X:
                atom_size = 1 ;
                nb_bits = atom_nb ;
                atom_nb = (int)((nb_bits - 1)/ 8) + 1 ;
                break ;
            case TFITS_ASCII_TYPE_A:
                atom_size = atom_nb ;
                break ;
            default:
                qfits_error("unrecognized type") ;
                qfits_table_close(tload) ;
                return NULL ;
                break ;
        }

        /* zero <-> TZERO */
        sprintf(keyword, "TZERO%d", i+1) ;
        if ((str_val=qfits_query_ext(filename, keyword, xtnum)) != NULL) {
            zero = (float)atof(str_val) ;
            zero_present = 1 ;
        } else {
            zero = (float)0.0 ;
            zero_present = 0 ;
        }

        /* scale <-> TSCAL */
        sprintf(keyword, "TSCAL%d", i+1) ;
        if ((str_val=qfits_query_ext(filename, keyword, xtnum)) != NULL) {
            scale = (float)atof(str_val) ;
            scale_present = 1 ;
        } else {
            scale = (float)1.0 ;
            scale_present = 0 ;
        }

        /* Fill the current column object */
        qfits_col_fill(curr_col, atom_nb, atom_dec_nb, atom_size, atom_type,
                label, unit, nullval, disp, zero_present, zero, scale_present,
                scale, offset_beg) ;

        /* Compute offset_beg but for the last column */
        if (i < tload->nc - 1) {
            if (table_type == QFITS_ASCIITABLE) {
                /* column width <-> TBCOLi and TBCOLi+1 */
                sprintf(keyword, "TBCOL%d", i+1) ;
                if ((str_val=qfits_query_ext(filename, keyword, xtnum))==NULL) {
                    qfits_error("cannot read [%s] in [%s]", keyword, filename);
                    qfits_table_close(tload);
                    return NULL ;
                }
                col_pos = atoi(qfits_pretty_string(str_val)) ;

                sprintf(keyword, "TBCOL%d", i+2) ;
                if ((str_val=qfits_query_ext(filename, keyword, xtnum))==NULL){
                    qfits_error("cannot read [%s] in [%s]", keyword, filename) ;
                    qfits_table_close(tload) ;
                    return NULL ;
                }
                next_col_pos = atoi(qfits_pretty_string(str_val)) ;
                offset_beg += (int)(next_col_pos - col_pos) ;
            } else if (table_type == QFITS_BINTABLE) {
                offset_beg += atom_nb * atom_size ;
            }
        }
        curr_col++ ;
    }

    /* Check that the theoritical data size is not far from the measured */
    /* one by more than 2880 */
    theory_size = qfits_compute_table_width(tload)*tload->nr ;
    if (data_size < theory_size) {
        qfits_error("Uncoherent data sizes") ;
        qfits_table_close(tload) ;
        return NULL ;
    }

    /* Return  */
    return tload ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Free a FITS table and associated pointers
  @param    t qfits_table to free
  @return    void
  Frees all memory associated to a qfits_table structure.
 */
/*----------------------------------------------------------------------------*/
void qfits_table_close(qfits_table * t)
{
    if (t==NULL) return ;
    if (t->nc>0) if (t->col!=NULL) qfits_free(t->col) ;
    qfits_free(t);
    return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Extract data from a column in a FITS table
  @param    th        Allocated qfits_table
  @param    colnum    Number of the column to extract (from 0 to colnum-1)
  @param    selection  boolean array to define the selected rows
  @return    unsigned char array

  If selection is NULL, select the complete column.

  Extract a column from a FITS table and return the data as a bytes
  array. The returned array type and size are determined by the
  column object in the qfits_table and by the selection parameter.

  Returned array size in bytes is:
  nbselected * col->natoms * col->atom_size

  Numeric types are correctly understood and byte-swapped if needed,
  to be converted to the local machine type.

  NULL values have to be handled by the caller.

  The returned object must be deallocated with qfits_free().
 */
/*----------------------------------------------------------------------------*/
unsigned char * qfits_query_column(
        const qfits_table   *   th,
        int                     colnum,
        const int           *   selection)
{
    char            *    start ;
    qfits_col       *   col ;
    int                    field_size ;
    unsigned char   *   array ;
    unsigned char   *   r ;
    unsigned char   *   inbuf ;
    int                 table_width ;
    int                 nb_rows ;
    size_t              size ;
    int                 i ;

    if (th->tab_w == -1) {
        /* Compute the table width in bytes */
        if ((table_width = qfits_compute_table_width(th)) == -1) {
            qfits_error("cannot compute the table width") ;
            return NULL ;
        }
    } else table_width = th->tab_w ;

    /* Compute the number of selected rows */
    nb_rows = 0 ;
    if (selection == NULL) {
        nb_rows = th->nr ;
    } else {
        for (i=0 ; i<th->nr ; i++) if (selection[i] == 1) nb_rows++ ;
    }

    /* Pointer to requested column */
    col = th->col + colnum ;

    /* Test if column is empty */
    if (nb_rows * col->atom_size * col->atom_nb == 0) col->readable = 0 ;

    /* Test if column is readable */
    if (col->readable == 0)  return NULL ;

    /* Compute the size in bytes of one field stored in the file */
    if ((field_size=qfits_table_get_field_size(th->tab_t,col))==-1) return NULL;

    /* Load input file */
    if ((start=qfits_falloc((char *)(th->filename), 0, &size))==NULL) {
        qfits_error("cannot open table for query [%s]", th->filename);
        return NULL ;
    }

    /* Allocate data array */
    array = qfits_malloc(nb_rows * field_size * sizeof(char)) ;

    /* Position the input pointer at the begining of the column data */
    r = array ;
    inbuf = (unsigned char*)start + col->off_beg ;

    /* Copy the values in array */
    if (selection == NULL) {
        /* No selection : get the complete column */
        for (i=0 ; i<th->nr ; i++) {
            /* Copy all atoms on this field into array */
            memcpy(r, inbuf, field_size);
            r += field_size ;
            /* Jump to next line */
            inbuf += table_width ;
        }
    } else {
        /* Get only the selected rows */
        for (i=0 ; i<th->nr ; i++) {
            if (selection[i] == 1) {
                /* Copy all atoms on this field into array */
                memcpy(r, inbuf, field_size);
                r += field_size ;
            }
            /* Jump to next line */
            inbuf += table_width ;
        }
    }
    qfits_fdealloc(start, 0, size) ;

    /* SWAP the bytes if necessary */
#ifndef WORDS_BIGENDIAN
    if (th->tab_t == QFITS_BINTABLE) {
        r = array ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            qfits_swap_bytes(r, col->atom_size);
            r += col->atom_size ;
        }
    }
#endif

     /* Return allocated and converted array */
    return array ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Extract consequtive values from a column in a FITS table
  @param    th        Allocated qfits_table
  @param    colnum    Number of the column to extract (from 0 to colnum-1)
  @param    start_ind   Index of the first row (0 for the first)
  @param    nb_rows     Number of rows to extract
  @return    unsigned char array
  Does the same as qfits_query_column() but on a consequtive sequence of rows
  Spares the overhead of the selection object allocation
  The returned object must be deallocated with qfits_free().
 */
/*----------------------------------------------------------------------------*/
unsigned char * qfits_query_column_seq(
        const qfits_table   *   th,
        int                     colnum,
        int                     start_ind,
        int                     nb_rows)
{
    char            *   start ;
    qfits_col       *   col ;
    int                 field_size ;
    unsigned char   *   array ;
    unsigned char   *   r ;
    unsigned char   *   inbuf ;
    int                 table_width ;
    size_t              size ;
    int                 i ;

    if (th->tab_w == -1) {
        /* Compute the table width in bytes */
        if ((table_width = qfits_compute_table_width(th)) == -1) {
            qfits_error("cannot compute the table width") ;
            return NULL ;
        }
    } else table_width = th->tab_w ;

    /* Check the validity of start_ind and nb_rows */
    if ((start_ind<0) || (start_ind+nb_rows>th->nr)) {
        qfits_error("bad start index and number of rows") ;
        return NULL ;
    }

    /* Pointer to requested column */
    col = th->col + colnum ;

    /* Test if column is empty */
    if (nb_rows * col->atom_size * col->atom_nb == 0) col->readable = 0 ;

    /* Test if column is readable */
    if (col->readable == 0)  return NULL ;

    /* Compute the size in bytes of one field stored in the file */
    if ((field_size=qfits_table_get_field_size(th->tab_t,col))==-1) return NULL;

    /* Load input file */
    if ((start=qfits_falloc((char *)(th->filename), 0, &size))==NULL) {
        qfits_error("cannot open table for query [%s]", th->filename);
        return NULL ;
    }

    /* Allocate data array */
    array = qfits_malloc(nb_rows * field_size * sizeof(char)) ;

    /* Position the input pointer at the begining of the column data */
    r = array ;
    inbuf = (unsigned char*)start + col->off_beg + table_width * start_ind ;

    /* Copy the values in array */
    /* Get only the selected rows */
    for (i=0 ; i<nb_rows ; i++) {
        /* Copy all atoms on this field into array */
        memcpy(r, inbuf, field_size);
        r += field_size ;
        /* Jump to next line */
        inbuf += table_width ;
    }
    qfits_fdealloc(start, 0, size) ;

    /* SWAP the bytes if necessary */
#ifndef WORDS_BIGENDIAN
    if (th->tab_t == QFITS_BINTABLE) {
        r = array ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            qfits_swap_bytes(r, col->atom_size);
            r += col->atom_size ;
        }
    }
#endif

     /* Return allocated and converted array */
    return array ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Extract binary data from a column in a FITS table
  @param    th        Allocated qfits_table
  @param    colnum    Number of the column to extract (from 0 to colnum-1)
  @param    selection  bollean array to identify selected rows
  @param    null_value    Value to return when a NULL value comes
  @return    Pointer to void *

  Extract a column from a FITS table and return the data as a generic
  void* array. The returned array type and size are determined by the
  column object in the qfits_table.

  Returned array size in bytes is:
  nb_selected * col->atom_nb * col->atom_size

  NULL values are recognized and replaced by the specified value.
  The returned object must be deallocated with qfits_free().
 */
/*----------------------------------------------------------------------------*/
void * qfits_query_column_data(
        const qfits_table   *   th,
        int                     colnum,
        const int           *   selection,
        const void          *    null_value)
{
    void            *    out_array ;
    qfits_col       *   col ;
    int                 nb_rows ;
    unsigned char    *    in_array ;
    char            *    field ;

    unsigned char        ucnull ;
    short                snull ;
    int                 inull ;
    double                dnull ;
    float                fnull ;

    int                 i ;

    /* Initialize */
    if (null_value == NULL) {
        inull  = (int)0 ;
        snull  = (short)0 ;
        ucnull = (unsigned char)0 ;
        fnull  = (float)0.0 ;
        dnull  = (double)0.0 ;
    } else {
        inull  = *(int*)null_value ;
        snull  = *(short*)null_value ;
        ucnull = *(unsigned char*)null_value ;
        fnull  = *(float*)null_value ;
        dnull  = *(double*)null_value ;
    }

    /* Get the number of selected rows */
    nb_rows = 0 ;
    if (selection == NULL) {
        nb_rows = th->nr ;
    } else {
        for (i=0 ; i<th->nr ; i++) if (selection[i] == 1) nb_rows++ ;
    }

    /* Pointer to requested column */
    col = th->col+colnum ;

    /* Test if column is readable */
    if (col->readable == 0) return NULL ;

    /* Handle each type separately */
    switch(col->atom_type) {
        case TFITS_ASCII_TYPE_A:
        out_array = (char*)qfits_query_column(th, colnum, selection) ;
        break ;

        case TFITS_ASCII_TYPE_I:
        in_array = (unsigned char*)qfits_query_column(th, colnum, selection) ;
        out_array = qfits_malloc(nb_rows*col->atom_size);
        field = qfits_malloc((col->atom_nb+1)*sizeof(char)) ;
        for (i=0 ; i<nb_rows ; i++) {
            /* Copy all atoms of the field into 'field' */
            memcpy(field, &in_array[i*col->atom_nb], col->atom_nb);
            field[col->atom_nb]=(char)0 ;
            /* Write the data in out_array */
            /* Test if a NULL val is encoutered */
            if (!strcmp(col->nullval, qfits_strstrip(field))) {
                ((int*)out_array)[i] = inull ;
            } else {
                ((int*)out_array)[i] = (int)atoi(field) ;
            }
        }
        qfits_free(field) ;
        qfits_free(in_array) ;
        break ;

        case TFITS_ASCII_TYPE_E:
        case TFITS_ASCII_TYPE_F:
        in_array = (unsigned char*)qfits_query_column(th, colnum, selection) ;
        out_array = qfits_malloc(nb_rows*col->atom_size);
        field = qfits_malloc((col->atom_nb+1)*sizeof(char)) ;
        for (i=0 ; i<nb_rows ; i++) {
            /* Copy all atoms of the field into 'field' */
            memcpy(field, &in_array[i*col->atom_nb], col->atom_nb);
            field[col->atom_nb]=(char)0 ;
            /* Write the data in out_array */
            /* Test if a NULL val is encoutered */
            if (!strcmp(col->nullval, qfits_strstrip(field))) {
                ((float*)out_array)[i] = fnull ;
            } else {
                /* Add the decimal handling */
                ((float*)out_array)[i]=(float)qfits_str2dec(field,
                                                             col->atom_dec_nb) ;
            }
        }
        qfits_free(field) ;
        qfits_free(in_array) ;
        break ;

        case TFITS_ASCII_TYPE_D:
        in_array = (unsigned char*)qfits_query_column(th, colnum, selection) ;
        out_array = qfits_malloc(nb_rows*col->atom_size);
        field = qfits_malloc((col->atom_nb+1)*sizeof(char)) ;
        for (i=0 ; i<nb_rows ; i++) {
            /* Copy all atoms of the field into 'field' */
            memcpy(field, &in_array[i*col->atom_nb], col->atom_nb);
            field[col->atom_nb]=(char)0 ;
            /* Write the data in out_array */
            /* Test if a NULL val is encoutered */
            if (!strcmp(col->nullval, field)) {
                ((double*)out_array)[i] = dnull ;
            } else {
                /* Add the decimal handling */
                ((double*)out_array)[i]=qfits_str2dec(field, col->atom_dec_nb) ;
            }
        }
        qfits_free(field) ;
        qfits_free(in_array) ;

        break ;
        case TFITS_BIN_TYPE_A:
        case TFITS_BIN_TYPE_L:
        out_array = (char*)qfits_query_column(th, colnum, selection) ;
        break ;

        case TFITS_BIN_TYPE_D:
        case TFITS_BIN_TYPE_M:
        out_array = (double*)qfits_query_column(th, colnum, selection) ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (qfits_isnan(((double*)out_array)[i]) ||
                    qfits_isinf(((double*)out_array)[i])) {
                ((double*)out_array)[i] = dnull ;
            }
        }
        break ;

        case TFITS_BIN_TYPE_E:
        case TFITS_BIN_TYPE_C:
        out_array = (float*)qfits_query_column(th, colnum, selection) ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (qfits_isnan(((float*)out_array)[i]) ||
                    qfits_isinf(((float*)out_array)[i])) {
                ((float*)out_array)[i] = fnull ;
            }
        }
        break ;

        case TFITS_BIN_TYPE_X:
        out_array = (unsigned char*)qfits_query_column(th, colnum, selection) ;
        break ;

        case TFITS_BIN_TYPE_B:
        out_array = (unsigned char*)qfits_query_column(th, colnum, selection) ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (((col->nullval)[0] != (char)0) &&
                (atoi(col->nullval) == (int)((unsigned char*)out_array)[i])) {
                ((unsigned char*)out_array)[i] = ucnull ;
            }
        }
        break ;

        case TFITS_BIN_TYPE_I:
        out_array = (short*)qfits_query_column(th, colnum, selection) ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (((col->nullval)[0] != (char)0) &&
                (atoi(col->nullval)==(int)((short*)out_array)[i])) {
                ((short*)out_array)[i] = snull ;
            }
        }
        break ;

        case TFITS_BIN_TYPE_J:
        out_array = (int*)qfits_query_column(th, colnum, selection) ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (((col->nullval)[0] != (char)0) &&
                (atoi(col->nullval)==((int*)out_array)[i])) {
                ((int*)out_array)[i] = inull ;
            }
        }
        break ;

        case TFITS_BIN_TYPE_P:
        out_array = (int*)qfits_query_column(th, colnum, selection) ;
        break ;

        default:
        qfits_error("unrecognized data type") ;
        return NULL ;
    }
    return out_array ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Extract binary data from a column in a FITS table
  @param    th        Allocated qfits_table
  @param    colnum    Number of the column to extract (from 0 to colnum-1)
  @param    start_ind   Index of the first row (0 for the first)
  @param    nb_rows     Number of rows to extract
  @param    null_value    Value to return when a NULL value comes
  @return    Pointer to void *
  Does the same as qfits_query_column_data() but on a consequtive sequence
  of rows.  Spares the overhead of the selection object allocation
  The returned object must be deallocated with qfits_free().
 */
/*----------------------------------------------------------------------------*/
void * qfits_query_column_seq_data(
        const qfits_table   *   th,
        int                     colnum,
        int                     start_ind,
        int                     nb_rows,
        const void          *   null_value)
{
    void            *    out_array ;
    qfits_col       *   col ;
    unsigned char    *    in_array ;
    char            *    field ;

    unsigned char        ucnull ;
    short                snull ;
    int                 inull ;
    double                dnull ;
    float                fnull ;

    int                 i ;

    /* Initialize */
    if (null_value == NULL) {
        inull  = (int)0 ;
        snull  = (short)0 ;
        ucnull = (unsigned char)0 ;
        fnull  = (float)0.0 ;
        dnull  = (double)0.0 ;
    } else {
        inull  = *(int*)null_value ;
        snull  = *(short*)null_value ;
        ucnull = *(unsigned char*)null_value ;
        fnull  = *(float*)null_value ;
        dnull  = *(double*)null_value ;
    }

    /* Pointer to requested column */
    col = th->col+colnum ;

    /* Test if column is readable */
    if (col->readable == 0) return NULL ;

    /* Handle each type separately */
    switch(col->atom_type) {
        case TFITS_ASCII_TYPE_A:
        out_array = (char*)qfits_query_column_seq(th, colnum, start_ind,
                nb_rows) ;
        break ;

        case TFITS_ASCII_TYPE_I:
        in_array = (unsigned char*)qfits_query_column_seq(th, colnum,
                start_ind, nb_rows) ;
        out_array = qfits_malloc(nb_rows*col->atom_size);
        field = qfits_malloc((col->atom_nb+1)*sizeof(char)) ;
        for (i=0 ; i<nb_rows ; i++) {
            /* Copy all atoms of the field into 'field' */
            memcpy(field, &in_array[i*col->atom_nb], col->atom_nb);
            field[col->atom_nb]=(char)0 ;
            /* Write the data in out_array */
            /* Test if a NULL val is encoutered */
            if (!strcmp(col->nullval, qfits_strstrip(field))) {
                ((int*)out_array)[i] = inull ;
            } else {
                ((int*)out_array)[i] = (int)atoi(field) ;
            }
        }
        qfits_free(field) ;
        qfits_free(in_array) ;
        break ;

        case TFITS_ASCII_TYPE_E:
        case TFITS_ASCII_TYPE_F:
        in_array = (unsigned char*)qfits_query_column_seq(th, colnum,
                start_ind, nb_rows) ;
        out_array = qfits_malloc(nb_rows*col->atom_size);
        field = qfits_malloc((col->atom_nb+1)*sizeof(char)) ;
        for (i=0 ; i<nb_rows ; i++) {
            /* Copy all atoms of the field into 'field' */
            memcpy(field, &in_array[i*col->atom_nb], col->atom_nb);
            field[col->atom_nb]=(char)0 ;
            /* Write the data in out_array */
            /* Test if a NULL val is encoutered */
            if (!strcmp(col->nullval, qfits_strstrip(field))) {
                ((float*)out_array)[i] = fnull ;
            } else {
                /* Add the decimal handling */
                ((float*)out_array)[i]=(float)qfits_str2dec(field,
                                                            col->atom_dec_nb) ;
            }
        }
        qfits_free(field) ;
        qfits_free(in_array) ;
        break ;

        case TFITS_ASCII_TYPE_D:
        in_array = (unsigned char*)qfits_query_column_seq(th, colnum,
                start_ind, nb_rows) ;
        out_array = qfits_malloc(nb_rows*col->atom_size);
        field = qfits_malloc((col->atom_nb+1)*sizeof(char)) ;
        for (i=0 ; i<nb_rows ; i++) {
            /* Copy all atoms of the field into 'field' */
            memcpy(field, &in_array[i*col->atom_nb], col->atom_nb);
            field[col->atom_nb]=(char)0 ;
            /* Write the data in out_array */
            /* Test if a NULL val is encoutered */
            if (!strcmp(col->nullval, qfits_strstrip(field))) {
                ((double*)out_array)[i] = dnull ;
            } else {
                /* Add the decimal handling */
                ((double*)out_array)[i]=qfits_str2dec(field, col->atom_dec_nb) ;
            }
        }
        qfits_free(field) ;
        qfits_free(in_array) ;
        break ;

        case TFITS_BIN_TYPE_A:
        case TFITS_BIN_TYPE_L:
        out_array = (char*)qfits_query_column_seq(th, colnum,
                start_ind, nb_rows) ;
        break ;

        case TFITS_BIN_TYPE_D:
        case TFITS_BIN_TYPE_M:
        out_array = (double*)qfits_query_column_seq(th, colnum,
                start_ind, nb_rows) ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (qfits_isnan(((double*)out_array)[i]) ||
                    qfits_isinf(((double*)out_array)[i])) {
                ((double*)out_array)[i] = dnull ;
            }
        }
        break ;

        case TFITS_BIN_TYPE_E:
        case TFITS_BIN_TYPE_C:
        out_array = (float*)qfits_query_column_seq(th, colnum,
                start_ind, nb_rows) ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (qfits_isnan(((float*)out_array)[i]) ||
                    qfits_isinf(((float*)out_array)[i])) {
                ((float*)out_array)[i] = fnull ;
            }
        }
        break ;

        case TFITS_BIN_TYPE_X:
        out_array = (unsigned char*)qfits_query_column_seq(th, colnum,
                start_ind, nb_rows) ;
        break ;

        case TFITS_BIN_TYPE_B:
        out_array = (unsigned char*)qfits_query_column_seq(th, colnum,
                start_ind, nb_rows) ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (((col->nullval)[0] != (char)0) &&
                (atoi(col->nullval)== (int)((unsigned char*)out_array)[i])) {
                ((unsigned char*)out_array)[i] = ucnull ;
            }
        }
        break ;

        case TFITS_BIN_TYPE_I:
        out_array = (short*)qfits_query_column_seq(th, colnum,
                start_ind, nb_rows) ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (((col->nullval)[0] != (char)0) &&
                (atoi(col->nullval)==(int)((short*)out_array)[i])) {
                ((short*)out_array)[i] = snull ;
            }
        }
        break ;

        case TFITS_BIN_TYPE_J:
        out_array = (int*)qfits_query_column_seq(th, colnum,
                start_ind, nb_rows) ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (((col->nullval)[0] != (char)0) &&
                (atoi(col->nullval)==((int*)out_array)[i])) {
                ((int*)out_array)[i] = inull ;
            }
        }
        break ;

        case TFITS_BIN_TYPE_P:
        out_array = (int*)qfits_query_column_seq(th, colnum,
                start_ind, nb_rows) ;
        break ;

        default:
        qfits_error("unrecognized data type") ;
        return NULL ;
    }
    return out_array ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Detect NULL values in a column
  @param    th        Allocated qfits_table
  @param    colnum    Number of the column to check (from 0 to colnum-1)
  @param    selection Array to identify selected rows
  @param    nb_vals Gives the size of the output array
  @param    nb_nulls Gives the number of detected null values
  @return     array with 1 for NULLs and 0 for non-NULLs
  The returned object must be deallocated with qfits_free().
 */
/*----------------------------------------------------------------------------*/
int * qfits_query_column_nulls(
        const qfits_table   *   th,
        int                     colnum,
        const int           *   selection,
        int                 *   nb_vals,
        int                 *   nb_nulls)
{
    int             *    out_array ;
    qfits_col       *   col ;
    unsigned char    *    in_array ;
    void            *    tmp_array ;
    char            *    field ;
    int                 nb_rows ;
    int                 i ;

    /* Initialize */
    *nb_nulls = 0 ;
    *nb_vals = 0 ;

    /* Get the number of selected rows */
    nb_rows = 0 ;
    if (selection == NULL) {
        nb_rows = th->nr ;
    } else {
       for (i=0 ; i<th->nr ; i++) if (selection[i] == 1) nb_rows++ ;
    }

    /* Pointer to requested column */
    col = th->col+colnum ;

    /* Test if column is readable */
    if (col->readable == 0) return NULL ;

    /* Handle each type separately */
    switch(col->atom_type) {
        case TFITS_ASCII_TYPE_A:
        case TFITS_ASCII_TYPE_D:
        case TFITS_ASCII_TYPE_E:
        case TFITS_ASCII_TYPE_F:
        case TFITS_ASCII_TYPE_I:
        in_array = (unsigned char*)qfits_query_column(th, colnum, selection) ;
        out_array = qfits_calloc(nb_rows, sizeof(int));
        *nb_vals = nb_rows ;
        field = qfits_malloc((col->atom_nb+1)*sizeof(char)) ;
        for (i=0 ; i<nb_rows ; i++) {
            /* Copy all atoms of the field into 'field' */
            memcpy(field, &in_array[i*col->atom_nb], col->atom_nb);
            field[col->atom_nb]=(char)0 ;
            /* Test if a NULL val is encoutered */
            if (!strcmp(col->nullval, qfits_strstrip(field))) {
                out_array[i] = 1 ;
                (*nb_nulls)++ ;
            }
        }
        qfits_free(field) ;
        if (in_array != NULL) qfits_free(in_array) ;
        break ;

        case TFITS_BIN_TYPE_A:
        /* No NULL values */
        out_array = qfits_calloc(nb_rows * col->atom_nb, sizeof(int)) ;
        *nb_vals = nb_rows * col->atom_nb ;
        break ;

        case TFITS_BIN_TYPE_L:
        case TFITS_BIN_TYPE_X:
        case TFITS_BIN_TYPE_P:
        /* No NULL values */
        out_array = qfits_calloc(nb_rows * col->atom_nb, sizeof(int)) ;
        *nb_vals = nb_rows * col->atom_nb ;
        break ;

        case TFITS_BIN_TYPE_D:
        case TFITS_BIN_TYPE_M:
        tmp_array = (double*)qfits_query_column(th, colnum, selection) ;
        out_array = qfits_calloc(nb_rows * col->atom_nb, sizeof(int)) ;
        *nb_vals = nb_rows * col->atom_nb ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (qfits_isnan(((double*)tmp_array)[i]) ||
                qfits_isinf(((double*)tmp_array)[i])) {
                out_array[i] = 1 ;
                (*nb_nulls)++ ;
            }
        }
        if (tmp_array != NULL) qfits_free(tmp_array) ;
        break ;

        case TFITS_BIN_TYPE_E:
        case TFITS_BIN_TYPE_C:
        tmp_array = (float*)qfits_query_column(th, colnum, selection) ;
        out_array = qfits_calloc(nb_rows * col->atom_nb, sizeof(int)) ;
        *nb_vals = nb_rows * col->atom_nb ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (qfits_isnan(((float*)tmp_array)[i]) ||
                qfits_isinf(((float*)tmp_array)[i])) {
                out_array[i] = 1 ;
                (*nb_nulls)++ ;
            }
        }
        if (tmp_array != NULL) qfits_free(tmp_array) ;
        break ;

        case TFITS_BIN_TYPE_B:
        tmp_array = (unsigned char*)qfits_query_column(th, colnum, selection) ;
        out_array = qfits_calloc(nb_rows * col->atom_nb, sizeof(int)) ;
        *nb_vals = nb_rows * col->atom_nb ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (((col->nullval)[0] != (char)0) &&
                (atoi(col->nullval)==(int)((unsigned char*)tmp_array)[i])) {
                out_array[i] = 1 ;
                (*nb_nulls)++ ;
            }
        }
        if (tmp_array != NULL) qfits_free(tmp_array) ;
        break ;

        case TFITS_BIN_TYPE_I:
        tmp_array = (short*)qfits_query_column(th, colnum, selection) ;
        out_array = qfits_calloc(nb_rows * col->atom_nb, sizeof(int)) ;
        *nb_vals = nb_rows * col->atom_nb ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (((col->nullval)[0] != (char)0) &&
                (atoi(col->nullval)==(int)((short*)tmp_array)[i])) {
                out_array[i] = 1 ;
                (*nb_nulls)++ ;
            }
        }
        if (tmp_array != NULL) qfits_free(tmp_array) ;
        break ;

        case TFITS_BIN_TYPE_J:
        tmp_array = (int*)qfits_query_column(th, colnum, selection) ;
        out_array = qfits_calloc(nb_rows * col->atom_nb, sizeof(int)) ;
        *nb_vals = nb_rows * col->atom_nb ;
        for (i=0 ; i<nb_rows * col->atom_nb ; i++) {
            if (((col->nullval)[0] != (char)0) &&
                (atoi(col->nullval)==((int*)tmp_array)[i])) {
                out_array[i] = 1 ;
                (*nb_nulls)++ ;
            }
        }
        if (tmp_array != NULL) qfits_free(tmp_array) ;
        break ;

        default:
        qfits_error("unrecognized data type") ;
        return NULL ;
    }
    return out_array ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a table to a FITS file with a given FITS header.
  @param    array            Data array.
  @param    table            table
  @param    fh              FITS header to insert in the output file.
  @return   -1 in error case, 0 otherwise
 */
/*----------------------------------------------------------------------------*/
int qfits_save_table_hdrdump(
        const void          **  array,
        const qfits_table   *   table,
        const qfits_header  *   fh)
{
    FILE        *   outfile ;
    const char  *   md5hash ;
    char            md5card[81];

    /* Open the destination file */
    if ((outfile = fopen(table->filename, "w")) == NULL) {
        qfits_error("cannot open file [%s]", table->filename) ;
        return -1 ;
    }
    /* Write the fits header in the file 'outname' */
    if (qfits_header_dump(fh, outfile) == -1) {
        qfits_error("cannot dump header in file") ;
        fclose(outfile) ;
        return -1 ;
    }
    /* Append the extension */
    if (table->tab_t == QFITS_BINTABLE) {
        if (qfits_table_append_bin_xtension(outfile, table, array) == -1) {
            qfits_error("in writing fits table") ;
            fclose(outfile) ;
            return -1 ;
        }
    } else if (table->tab_t == QFITS_ASCIITABLE) {
        if (qfits_table_append_ascii_xtension(outfile, table, array) == -1) {
            qfits_error("in writing fits table") ;
            fclose(outfile) ;
            return -1 ;
        }
    } else {
        qfits_error("Unrecognized table type") ;
        fclose(outfile) ;
        return -1 ;
    }
    fclose(outfile) ;
    /* Update MD5 keyword */
    if (strcmp(table->filename, "STDOUT")) {
        md5hash = qfits_datamd5(table->filename);
        if (md5hash==NULL) {
            qfits_error("computing MD5 signature for output file %s",
                    table->filename);
            return -1 ;
        }
        sprintf(md5card, "DATAMD5 = '%s' / MD5 checksum", md5hash);
        qfits_replace_card(table->filename, "DATAMD5", md5card);
    }
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Appends a std extension header + data to a FITS table file.
  @param    outfile        Pointer to (opened) file ready for writing.
  @param    t            Pointer to qfits_table
  @param    data        Table data to write
  @return    int 0 if Ok, -1 otherwise

  Dumps a FITS table to a file. The whole table described by qfits_table, and
  the data arrays contained in 'data' are dumped to the file. An extension
  header is produced with all keywords needed to describe the table, then the
  data is dumped to the file.
  The output is then padded to reach a multiple of 2880 bytes in size.
  Notice that no main header is produced, only the extension part.
 */
/*----------------------------------------------------------------------------*/
int qfits_table_append_xtension(
        FILE                *   outfile,
        const qfits_table   *   t,
        const void          **  data)
{
    /* Append the extension */
    if (t->tab_t == QFITS_BINTABLE) {
        if (qfits_table_append_bin_xtension(outfile, t, data) == -1) {
            qfits_error("in writing fits table") ;
            return -1 ;
        }
    } else if (t->tab_t == QFITS_ASCIITABLE) {
        if (qfits_table_append_ascii_xtension(outfile, t, data) == -1) {
            qfits_error("in writing fits table") ;
            return -1 ;
        }
    } else {
        qfits_error("Unrecognized table type") ;
        return -1 ;
    }
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Appends a specified extension header + data to a FITS table file.
  @param    outfile        Pointer to (opened) file ready for writing.
  @param    t            Pointer to qfits_table
  @param    data        Table data to write
  @param    hdr         Specified extension header
  @return    int 0 if Ok, -1 otherwise

  Dumps a FITS table to a file. The whole table described by qfits_table, and
  the data arrays contained in 'data' are dumped to the file following the
  specified fits header.
  The output is then padded to reach a multiple of 2880 bytes in size.
  Notice that no main header is produced, only the extension part.
 */
/*----------------------------------------------------------------------------*/
int qfits_table_append_xtension_hdr(
        FILE                *   outfile,
        const qfits_table   *   t,
        const void          **  data,
        const qfits_header  *   hdr)
{
    /* Write the fits header in the file  */
    if (qfits_header_dump(hdr, outfile) == -1) {
        qfits_error("cannot dump header in file") ;
        return -1 ;
    }

    /* Append the data to the file */
    return qfits_table_append_data(outfile, t, data) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    given a col and a row, find out the string to write for display
  @param    table    table structure
  @param    col_id    col id (0 -> nbcol-1)
  @param    row_id    row id (0 -> nrow-1)
  @param    use_zero_scale    Flag to use or not zero and scale
  @return    the string

  This function is highly inefficient, it should not be used in loops to
  display a complete table. It is more to get one field from time to
  time, or for debugging puposes.
  The returned object must be deallocated with qfits_free().
 */
/*----------------------------------------------------------------------------*/
char * qfits_table_field_to_string(
        const qfits_table   *   table,
        int                     col_id,
        int                     row_id,
        int                     use_zero_scale)
{
    char    *    str ;

    switch (table->tab_t) {
        case QFITS_BINTABLE:
            str = qfits_bintable_field_to_string(table, col_id, row_id,
                    use_zero_scale) ;
            break ;

        case QFITS_ASCIITABLE:
            str = qfits_asciitable_field_to_string(table, col_id, row_id,
                    use_zero_scale) ;
            break ;
        default:
            qfits_error("Table type not recognized") ;
            return NULL ;
            break ;
    }
    return str ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the table width in bytes from the columns infos
  @param    th        Allocated qfits_table
  @return    the width (-1 in error case)
 */
/*----------------------------------------------------------------------------*/
static int qfits_compute_table_width(const qfits_table * th)
{
    int             width ;
    qfits_col   *   curr_col ;
    int             i ;

    /* Initialize */
    width = 0 ;

    /* Loop on all columns and get column descriptions  */
    curr_col = th->col ;
    for (i=0 ; i<th->nc ; i++) {
        if (th->tab_t == QFITS_ASCIITABLE) {
            width += curr_col->atom_nb ;
        } else if (th->tab_t == QFITS_BINTABLE) {
            width += curr_col->atom_nb * curr_col->atom_size ;
        }
        curr_col++ ;
    }
    return width ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    given a col and a row, find out the string to write for display
  @param    table    table structure
  @param    col_id    col id (0 -> nbcol-1)
  @param    row_id    row id (0 -> nrow-1)
  @param    use_zero_scale    Flag to use or not zero and scale
  @return    the string to write

  The returned object must be deallocated with qfits_free().
  ASCII tables specific
 */
/*----------------------------------------------------------------------------*/
static char * qfits_asciitable_field_to_string(
        const qfits_table   *   table,
        int                     col_id,
        int                     row_id,
        int                     use_zero_scale)
{
    qfits_col       *       col ;
    char            *       ccol ;
    int             *       icol ;
    float           *       fcol ;
    double          *       dcol ;
    char                    ctmp[512];
    char            *       stmp ;
    int                     field_size ;
    void            *        field ;
    int             *       selection ;

    /* Test inputs */
    if (table->tab_t != QFITS_ASCIITABLE) return NULL ;

    /* Initialize */
    ctmp[0] = (char)0 ;

    /* Set selection to select the requested row */
    selection = qfits_calloc(table->nr, sizeof(int)) ;
    selection[row_id] = 1 ;

    /* Load the column data */
    if ((field = qfits_query_column_data(table, col_id, selection,
                    NULL)) == NULL) return NULL ;
    qfits_free(selection) ;

    /* Set reference to the column */
    col = table->col + col_id ;

    /* Compute field size and allocate stmp */
    if (col->atom_nb > ELEMENT_MAX_DISPLAY_SIZE) field_size = col->atom_nb + 1 ;
    else field_size = ELEMENT_MAX_DISPLAY_SIZE ;
    stmp = qfits_malloc(field_size * sizeof(char)) ;
    stmp[0] = (char)0 ;

    /* Get the string to write according to the type */
    switch(col->atom_type) {
        case TFITS_ASCII_TYPE_A:
            ccol = (char*)field ;
            strncpy(ctmp, ccol, col->atom_nb);
            ctmp[col->atom_nb] = (char)0;
            strcpy(stmp, ctmp);
            break ;

        case TFITS_ASCII_TYPE_I:
            icol = (int*)field ;
            /* Two cases: use col->zero and col->scale or not */
            if (col->zero_present && col->scale_present && use_zero_scale) {
                sprintf(stmp, "%f", (float)(col->zero +
                            (float)icol[0] * col->scale)) ;
            } else {
                sprintf(stmp, "%d", icol[0]);
            }
            break ;

        case TFITS_ASCII_TYPE_E:
        case TFITS_ASCII_TYPE_F:
            fcol = (float*)field ;
            /* Two cases: use col->zero and col->scale or not */
            if (col->zero_present && col->scale_present && use_zero_scale) {
                sprintf(stmp, "%f", (float)(col->zero +
                            fcol[0] * col->scale)) ;
            } else {
                sprintf(stmp, "%f", fcol[0]);
            }
            break ;

        case TFITS_ASCII_TYPE_D:
            dcol = (double*)field ;
            /* Two cases: use col->zero and col->scale or not */
            if (col->zero_present && col->scale_present && use_zero_scale) {
                sprintf(stmp, "%f", (float)(col->zero +
                        (float)dcol[0] * col->scale)) ;
            } else {
                sprintf(stmp, "%g", dcol[0]) ;
            }
            break ;
        default:
            qfits_warning("Type not recognized") ;
            break ;
    }

    /* Free and return */
    qfits_free(field) ;
    return stmp ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Given a col and a row, find out the string to write for display
  @param    table    table structure
  @param    col_id    col id (0 -> nbcol-1)
  @param    row_id    row id (0 -> nrow-1)
  @param    use_zero_scale    Flag to use or not zero and scale
  @return    the allocated string to write

  The returned object must be deallocated with qfits_free().
  BIN tables specific
 */
/*----------------------------------------------------------------------------*/
static char * qfits_bintable_field_to_string(
        const qfits_table   *   table,
        int                     col_id,
        int                     row_id,
        int                     use_zero_scale)
{
    qfits_col       *       col ;
    unsigned char   *       uccol ;
    char            *       ccol ;
    int             *       icol ;
    short           *       scol ;
    float           *       fcol ;
    double          *       dcol ;
    char                    ctmp[512];
    char            *       stmp ;
    int                     field_size ;
    void            *        field ;
    int             *       selection ;

    int                     i ;

    /* Test inputs */
    if (table->tab_t != QFITS_BINTABLE) return NULL ;

   /* Initialize */
    ctmp[0] = (char)0 ;

    /* Set selection to select the requested row */
    selection = qfits_calloc(table->nr, sizeof(int)) ;
    selection[row_id] = 1 ;

    /* Load the data column */
    if ((field = qfits_query_column_data(table, col_id, selection,
                    NULL)) == NULL) {
        qfits_free(selection) ;
        return NULL ;
    }
    qfits_free(selection) ;

    /* Set reference to the column */
    col = table->col + col_id ;

    /* Compute field size and allocate stmp */
    field_size = col->atom_nb * ELEMENT_MAX_DISPLAY_SIZE ;
    stmp = qfits_malloc(field_size * sizeof(char)) ;
    stmp[0] = (char)0 ;

    /* Get the string to write according to the type */
    switch(col->atom_type) {
        case TFITS_BIN_TYPE_A:
        ccol = (char*)field ;
        strncpy(ctmp, ccol, col->atom_size * col->atom_nb) ;
        ctmp[col->atom_size*col->atom_nb] = (char)0 ;
        strcpy(stmp, ctmp) ;
        break ;

        case TFITS_BIN_TYPE_B:
        uccol = (unsigned char*)field ;
        /* Two cases: use col->zero and col->scale or not */
        if (col->zero_present && col->scale_present && use_zero_scale) {
            /* For each atom of the column */
            for (i=0 ; i<col->atom_nb-1 ; i++) {
                sprintf(ctmp, "%f, ", (float)(col->zero +
                        (float)uccol[i] * col->scale)) ;
                strcat(stmp, ctmp) ;
            }
            /* Handle the last atom differently: no ',' */
            sprintf(ctmp, "%f", (float)(col->zero +
                (float)uccol[col->atom_nb-1]*col->scale)) ;
            strcat(stmp, ctmp) ;
        } else {
            /* For each atom of the column */
            for (i=0 ; i<col->atom_nb-1 ; i++) {
                sprintf(ctmp, "%d, ", (int)uccol[i]) ;
                strcat(stmp, ctmp) ;
            }
            /* Handle the last atom differently: no ',' */
            sprintf(ctmp,"%d",(int)uccol[col->atom_nb-1]);
            strcat(stmp, ctmp) ;
        }
        break ;

        case TFITS_BIN_TYPE_D:
        case TFITS_BIN_TYPE_M:
        dcol = (double*)field ;
        /* Two cases: use col->zero and col->scale or not */
        if (col->zero_present && col->scale_present && use_zero_scale) {
            /* For each atom of the column */
            for (i=0 ; i<col->atom_nb-1 ; i++) {
                sprintf(ctmp, "%g, ", (double)((double)col->zero +
                        dcol[i] * (double)col->scale)) ;
                strcat(stmp, ctmp) ;
            }
            /* Handle the last atom differently: no ',' */
            sprintf(ctmp, "%g", (double)((double)col->zero +
                dcol[col->atom_nb-1] * (double)col->scale));
            strcat(stmp, ctmp) ;
        } else {
            /* For each atom of the column */
            for (i=0 ; i<col->atom_nb-1 ; i++) {
                sprintf(ctmp, "%g, ", dcol[i]) ;
                strcat(stmp, ctmp) ;
            }
            /* Handle the last atom differently: no ',' */
            sprintf(ctmp, "%g", dcol[col->atom_nb-1]) ;
            strcat(stmp, ctmp) ;
        }
        break ;

        case TFITS_BIN_TYPE_E:
        case TFITS_BIN_TYPE_C:
        fcol = (float*)field ;
        /* Two cases: use col->zero and col->scale or not */
        if (col->zero_present && col->scale_present && use_zero_scale) {
            /* For each atom of the column */
            for (i=0 ; i<col->atom_nb-1 ; i++) {
                sprintf(ctmp, "%f, ", (float)(col->zero +
                        (float)fcol[i] * col->scale)) ;
                strcat(stmp, ctmp) ;
            }
            /* Handle the last atom differently: no ',' */
            sprintf(ctmp, "%f", (float)(col->zero +
                (float)fcol[col->atom_nb-1] * col->scale)) ;
            strcat(stmp, ctmp) ;
        } else {
            /* For each atom of the column */
            for (i=0 ; i<col->atom_nb-1 ; i++) {
                sprintf(ctmp, "%f, ", fcol[i]) ;
                strcat(stmp, ctmp) ;
            }
            /* Handle the last atom differently: no ',' */
            sprintf(ctmp, "%f", fcol[col->atom_nb-1]) ;
            strcat(stmp, ctmp) ;
        }
        break ;

        case TFITS_BIN_TYPE_I:
        scol = (short*)field ;
        /* Two cases: use col->zero and col->scale or not */
        if (col->zero_present && col->scale_present && use_zero_scale) {
            /* For each atom of the column */
            for (i=0 ; i<col->atom_nb-1 ; i++) {
                sprintf(ctmp, "%f, ", (float)(col->zero +
                        (float)scol[i] * col->scale)) ;
                strcat(stmp, ctmp) ;
            }
            /* Handle the last atom differently: no ',' */
            sprintf(ctmp, "%f", (float)(col->zero +
                (float)scol[col->atom_nb-1] * col->scale)) ;
            strcat(stmp, ctmp) ;
        } else {
            /* For each atom of the column */
            for (i=0 ; i<col->atom_nb-1 ; i++) {
                sprintf(ctmp, "%d, ", (int)scol[i]) ;
                strcat(stmp, ctmp) ;
            }
            /* Handle the last atom differently: no ',' */
            sprintf(ctmp, "%d",(int)scol[col->atom_nb-1]);
            strcat(stmp, ctmp) ;
        }
        break ;

        case TFITS_BIN_TYPE_J:
        icol = (int*)field ;
        /* Two cases: use col->zero and col->scale or not */
        if (col->zero_present && col->scale_present && use_zero_scale) {
            /* For each atom of the column */
            for (i=0 ; i<col->atom_nb-1 ; i++) {
                sprintf(ctmp, "%f, ", (float)(col->zero +
                        (float)icol[i] * col->scale)) ;
                strcat(stmp, ctmp) ;
            }
            /* Handle the last atom differently: no ',' */
            sprintf(ctmp, "%f", (float)(col->zero +
                (float)icol[col->atom_nb-1] * col->scale)) ;
            strcat(stmp, ctmp) ;
        } else {
            /* For each atom of the column */
            for (i=0 ; i<col->atom_nb-1 ; i++) {
                sprintf(ctmp, "%d, ", (int)icol[i]) ;
                strcat(stmp, ctmp) ;
            }
            /* Handle the last atom differently: no ',' */
            sprintf(ctmp, "%d",(int)icol[col->atom_nb-1]);
            strcat(stmp, ctmp) ;
        }
        break ;

        case TFITS_BIN_TYPE_L:
        ccol = (char*)field ;
        /* For each atom of the column */
        for (i=0 ; i<col->atom_nb-1 ; i++) {
            sprintf(ctmp, "%c, ", ccol[i]) ;
            strcat(stmp, ctmp) ;
        }
        /* Handle the last atom differently: no ',' */
        sprintf(ctmp, "%c", ccol[col->atom_nb-1]) ;
        strcat(stmp, ctmp) ;
        break ;

        case TFITS_BIN_TYPE_X:
        uccol = (unsigned char*)field ;
        /* For each atom of the column */
        for (i=0 ; i<col->atom_nb-1 ; i++) {
            sprintf(ctmp, "%d, ", uccol[i]) ;
            strcat(stmp, ctmp) ;
        }
        /* Handle the last atom differently: no ',' */
        sprintf(ctmp, "%d", uccol[col->atom_nb-1]) ;
        strcat(stmp, ctmp) ;
        break ;

        case TFITS_BIN_TYPE_P:
        icol = (int*)field ;
        /* For each atom of the column */
        for (i=0 ; i<col->atom_nb-1 ; i++) {
            sprintf(ctmp, "%d, ", (int)icol[i]) ;
            strcat(stmp, ctmp) ;
        }
        /* Handle the last atom differently: no ',' */
        sprintf(ctmp, "%d",(int)icol[col->atom_nb-1]);
        strcat(stmp, ctmp) ;
        break ;

        default:
        qfits_warning("Type not recognized") ;
        break ;
    }
    qfits_free(field) ;
    return stmp ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Remove blanks at the beginning and the end of a string.
  @param    s   String to parse.
  @return   ptr to statically allocated string.

  This function returns a pointer to a statically allocated string,
  which is identical to the input string, except that all blank
  characters at the end and the beg. of the string have been removed.
  Do not free or modify the returned string!
  Since the returned string is statically allocated, it will be modified at
  each function call (not re-entrant).
 */
/*----------------------------------------------------------------------------*/
static char * qfits_strstrip(const char * s)
{
    static char l[1024+1];
    char * last ;

    if (s==NULL) return NULL ;

    while (isspace((int)*s) && *s) s++;

    memset(l, 0, 1024+1);
    strcpy(l, s);
    last = l + strlen(l);
    while (last > l) {
        if (!isspace((int)*(last-1)))
            break ;
        last -- ;
    }
    *last = (char)0;

    return (char*)l ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Make a double out of a string and a number of decimals
  @param    to_format   the string to convert
  @param    nb_dec      the number of decimals
  @return   the double
  A field with 123 of type F3.1 actually contains 12.3
  This is handled by this function.
 */
/*----------------------------------------------------------------------------*/
static double qfits_str2dec(
        const char  *   to_format,
        int             nb_dec)
{
    double      val ;
    int         i ;

    /* Test entries */
    if (to_format == NULL) return 0.00 ;

    val = (double)atof(to_format) ;
    /* First handle case where there are no decimals or the dot is there */
    if ((strstr(to_format, ".") == NULL) && (nb_dec > 0)) {
        for (i=0 ; i<nb_dec ; i++) val /=10 ;
    }
    return val ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Parse a FITS type
  @param    str        string read in the FITS header (e.g. TFORM value)
  @param    nb        pointer to the number
  @param    dec_nb  pointer to the number of decimals
  @param    type    pointer to the type
  @param    table_type    Table type (BIN, ASCII, ...)
  @return    0 if ok, -1 otherwise

  This functions reads the input string and uses it to update nb and type
 */
/*----------------------------------------------------------------------------*/
static int qfits_table_interpret_type(
        const char  *   str,
        int         *   nb,
        int         *   dec_nb,
        tfits_type  *   type,
        int             table_type)
{
    char        type_c ;

    *dec_nb = 0 ;
    if (table_type == QFITS_BINTABLE) {
        if (sscanf(str, "%d%c", nb, &type_c) == 0) {
            /* nb is 1 by default */
            if (sscanf(str, "%c", &type_c) == 0) {
                qfits_error("cannot interpret this type: %s", str) ;
                return -1 ;
            }
            *nb = 1 ;
        }
        switch(type_c) {
            case 'A': *type = TFITS_BIN_TYPE_A ; break ;
            case 'B': *type = TFITS_BIN_TYPE_B ; break ;
            case 'C': *type = TFITS_BIN_TYPE_C ; break ;
            case 'D': *type = TFITS_BIN_TYPE_D ; break ;
            case 'E': *type = TFITS_BIN_TYPE_E ; break ;
            case 'I': *type = TFITS_BIN_TYPE_I ; break ;
            case 'J': *type = TFITS_BIN_TYPE_J ; break ;
            case 'L': *type = TFITS_BIN_TYPE_L ; break ;
            case 'M': *type = TFITS_BIN_TYPE_M ; break ;
            case 'P': *type = TFITS_BIN_TYPE_P ; break ;
            case 'X': *type = TFITS_BIN_TYPE_X ; break ;
            default: return -1 ;
        }
    } else if (table_type == QFITS_ASCIITABLE) {
        if (sscanf(str, "%c%d.%d", &type_c, nb, dec_nb) == 0) {
            qfits_error("cannot interpret this type: %s", str) ;
            return -1 ;
        }
        switch(type_c) {
            case 'A': *type = TFITS_ASCII_TYPE_A ; break ;
            case 'D': *type = TFITS_ASCII_TYPE_D ; break ;
            case 'E': *type = TFITS_ASCII_TYPE_E ; break ;
            case 'F': *type = TFITS_ASCII_TYPE_F ; break ;
            case 'I': *type = TFITS_ASCII_TYPE_I ; break ;
            default: return -1 ;
        }
    } else {
        qfits_error("unrecognized table type") ;
        return -1 ;
    }
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Generate a FITS type string
  @param    col        input column
  @return    The string to write to TFORM
 */
/*----------------------------------------------------------------------------*/
static char * qfits_build_format(const qfits_col * col)
{
    static char sval[10] ;
    int         nb ;

    switch (col->atom_type) {
        case TFITS_ASCII_TYPE_A:
            nb=sprintf(sval, "A%d.%d", col->atom_nb, col->atom_dec_nb) ; break ;
        case TFITS_ASCII_TYPE_D:
            nb=sprintf(sval, "D%d.%d", col->atom_nb, col->atom_dec_nb) ; break ;
        case TFITS_ASCII_TYPE_E:
            nb=sprintf(sval, "E%d.%d", col->atom_nb, col->atom_dec_nb) ; break ;
        case TFITS_ASCII_TYPE_I:
            nb=sprintf(sval, "I%d.%d", col->atom_nb, col->atom_dec_nb) ; break ;
        case TFITS_ASCII_TYPE_F:
            nb=sprintf(sval, "F%d.%d", col->atom_nb, col->atom_dec_nb) ; break ;
        case TFITS_BIN_TYPE_D: nb=sprintf(sval, "%dD", col->atom_nb) ; break ;
        case TFITS_BIN_TYPE_E: nb=sprintf(sval, "%dE", col->atom_nb) ; break ;
        case TFITS_BIN_TYPE_I: nb=sprintf(sval, "%dI", col->atom_nb) ; break ;
        case TFITS_BIN_TYPE_A: nb=sprintf(sval, "%dA", col->atom_nb) ; break ;
        case TFITS_BIN_TYPE_B: nb=sprintf(sval, "%dB", col->atom_nb) ; break ;
        case TFITS_BIN_TYPE_C: nb=sprintf(sval, "%dC",
                                       (int)(col->atom_nb/2)) ; break ;
        case TFITS_BIN_TYPE_J: nb=sprintf(sval, "%dJ", col->atom_nb) ; break ;
        case TFITS_BIN_TYPE_L: nb=sprintf(sval, "%dL", col->atom_nb) ; break ;
        case TFITS_BIN_TYPE_M: nb=sprintf(sval, "%dM",
                                       (int)(col->atom_nb/2)) ; break ;
        case TFITS_BIN_TYPE_P: nb=sprintf(sval, "%dP",
                                       (int)(col->atom_nb/2)) ; break ;
        case TFITS_BIN_TYPE_X: nb=sprintf(sval, "%dX",
                                       8*col->atom_nb) ; break ;
        default: return NULL ;
    }
    sval[nb] = (char)0 ;
    return sval ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Appends a std extension header + data to a FITS BIN table file.
  @param    outfile        Pointer to (opened) file ready for writing.
  @param    t            Pointer to qfits_table
  @param    data        Table data to write
  @return    int 0 if Ok, -1 otherwise

  Dumps a FITS table to a file. The whole table described by qfits_table, and
  the data arrays contained in 'data' are dumped to the file. An extension
  header is produced with all keywords needed to describe the table, then the
  data is dumped to the file.
  The output is then padded to reach a multiple of 2880 bytes in size.
  Notice that no main header is produced, only the extension part.
 */
/*----------------------------------------------------------------------------*/
static int qfits_table_append_bin_xtension(
        FILE                *   outfile,
        const qfits_table   *   t,
        const void          **  data)
{
    qfits_header    *    fh ;

    if ((fh=qfits_table_ext_header_default(t)) == NULL) {
        qfits_error("cannot create new fits header") ;
        return -1 ;
    }

    /* Write the fits header in the file  */
    if (qfits_header_dump(fh, outfile) == -1) {
        qfits_error("cannot dump header in file") ;
        qfits_header_destroy(fh) ;
        fclose(outfile) ;
        return -1 ;
    }
    qfits_header_destroy(fh) ;

    /* Append the data to the file */
    return qfits_table_append_data(outfile, t, data) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Appends an extension header + data to a FITS ASCII table file.
  @param    outfile        Pointer to (opened) file ready for writing.
  @param    t            Pointer to qfits_table
  @param    data        Table data to write
  @return    int 0 if Ok, -1 otherwise

  Dumps a FITS table to a file. The whole table described by
  qfits_table, and the data arrays contained in 'data' are dumped to
  the file. An extension header is produced with all keywords needed
  to describe the table, then the data is dumped to the file.

  The output is then padded to reach a multiple of 2880 bytes in size.

  Notice that no main header is produced, only the extension part.
 */
/*----------------------------------------------------------------------------*/
static int qfits_table_append_ascii_xtension(
        FILE                *   outfile,
        const qfits_table   *   t,
        const void          **  data)
{
    qfits_header    *    fh ;

    if ((fh=qfits_table_ext_header_default(t)) == NULL) {
        qfits_error("cannot create new fits header") ;
        return -1 ;
    }

    /* Write the fits header in the file  */
    if (qfits_header_dump(fh, outfile) == -1) {
        qfits_error("cannot dump header in file") ;
        qfits_header_destroy(fh) ;
        return -1 ;
    }
    qfits_header_destroy(fh) ;

    /* Append the data to the file */
    return qfits_table_append_data(outfile, t, data) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Appends data to a FITS table file.
  @param    outfile        Pointer to (opened) file ready for writing.
  @param    t            Pointer to qfits_table
  @param    data        Table data to write
  @return    int 0 if Ok, -1 otherwise

  Dumps the data part of a FITS table to a file. The primary header, as well as
  the extension header are supposed to be already there (and properly padded).
  The output is then padded to reach a multiple of 2880 bytes in size.
 */
/*----------------------------------------------------------------------------*/
static int qfits_table_append_data(
        FILE                *   outfile,
        const qfits_table   *   t,
        const void          **  data)
{
    qfits_col       *   curr_col ;
    char                field[1024] ;
    char            *   line ;
    unsigned char   **  array ;
    unsigned char   *   r ;
    unsigned char   *   inbuf ;
    int                 writt_char ;
    int                 nb_blanks ;
    int                 field_size ;
    int                 i, j ;

    /* Write DATA */
    array = qfits_malloc(t->nc*sizeof(unsigned char *)) ;

    curr_col = t->col ;
    for (i=0 ; i<t->nc ; i++) {
        /* Compute the field size */
        field_size = qfits_table_get_field_size(t->tab_t, curr_col) ;

        /* Copy data from data to array (unsigned char) */
        array[i] = qfits_malloc(t->nr * field_size) ;
        r = (unsigned char *)array[i] ;
        inbuf = (unsigned char *)(data[i]) ;

        /* Copy the data */
        if (t->tab_t == QFITS_ASCIITABLE) {
            /* ASCII table */
            for (j=0 ; j<t->nr ; j++) {
                switch(curr_col->atom_type) {
                    case TFITS_ASCII_TYPE_A :
                        strncpy(field, (char*)inbuf, curr_col->atom_nb) ;
                        field[curr_col->atom_nb] = (char)0 ;
                        inbuf += curr_col->atom_nb ;
                        break ;
                    case TFITS_ASCII_TYPE_D :
                        memset(field, ' ', curr_col->atom_nb) ;
                        sprintf(field, "%g", ((double*)data[i])[j]) ;
                        field[curr_col->atom_nb] = (char)0 ;
                        break ;
                    case TFITS_ASCII_TYPE_E :
                    case TFITS_ASCII_TYPE_F :
                        memset(field, ' ', curr_col->atom_nb) ;
                        sprintf(field, "%f", ((float*)data[i])[j]) ;
                        field[curr_col->atom_nb] = (char)0 ;
                        break ;
                    case TFITS_ASCII_TYPE_I :
                        memset(field, ' ', curr_col->atom_nb) ;
                        sprintf(field, "%d", ((int*)data[i])[j]) ;
                        field[curr_col->atom_nb] = (char)0 ;
                        break ;
                    default:
                        break ;
                }
                memcpy(r, field, curr_col->atom_nb) ;
                r += (curr_col->atom_nb) ;
            }
        } else if (t->tab_t == QFITS_BINTABLE) {
            /* BINARY table */
            for (j=0 ; j<t->nr ; j++) {
                memcpy(r, inbuf, field_size) ;
                inbuf += field_size ;
                r += field_size ;
            }

            /* Byte swapping needed if on a little-endian machine */
#ifndef WORDS_BIGENDIAN
            r = array[i] ;
            for (j=0 ; j<t->nr * curr_col->atom_nb ; j++) {
                qfits_swap_bytes(r, curr_col->atom_size);
                r += curr_col->atom_size ;
            }
#endif
        } else return -1 ;
        curr_col++ ;
    }

    /* Write to the outfile */
    writt_char = 0 ;
    for (i=0 ; i<t->nr ; i++) {
        curr_col = t->col ;
        for (j=0 ; j<t->nc ; j++) {
            field_size = qfits_table_get_field_size(t->tab_t, curr_col) ;
            r = array[j] + field_size * i ;
            line = (char *)qfits_calloc (field_size+1, sizeof (char)) ;
            memcpy(line, r, field_size) ;
            line[field_size] = (char)0 ;
            fwrite(line, 1, field_size, outfile) ;
            writt_char += field_size ;
            curr_col++ ;
            qfits_free(line) ;
        }
    }

    /* Complete with blanks to FITS_BLOCK_SIZE characters */
    if (writt_char % FITS_BLOCK_SIZE) {
        nb_blanks = FITS_BLOCK_SIZE - (writt_char%FITS_BLOCK_SIZE) ;
        for (i=1 ; i<=nb_blanks ; i++) fwrite(" ", 1, 1, outfile) ;
    }

    /* Free and return  */
    for(i=0 ; i<t->nc ; i++) {
        if (array[i] != NULL) qfits_free(array[i]) ;
    }
    qfits_free(array) ;
    return  0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the size in bytes of a field
  @param    type    table type
  @param    col     a column
  @return   the size
 */
/*----------------------------------------------------------------------------*/
static int qfits_table_get_field_size(
        int                 type,
        const qfits_col *   col)
{
    int     field_size ;

    switch (type) {
        case QFITS_BINTABLE:
            field_size = col->atom_nb * col->atom_size ;
            break ;
        case QFITS_ASCIITABLE:
            field_size = col->atom_nb ;
            break ;
        default:
            qfits_warning("unrecognized table type") ;
            field_size = -1 ;
    }
    return field_size ;
}
#endif
