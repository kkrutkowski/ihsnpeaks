/* Copyright 2019 Herman ten Brugge
 *
 * Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
 * http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
 * <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
 * option. This file may not be copied, modified, or distributed
 * except according to those terms.
 */

#ifndef __FAST_STDIO_H
#define __FAST_STDIO_H

#include <string.h>
#include <ctype.h>
#include <math.h>
#include <locale.h>
#include <inttypes.h>

#if defined (__cplusplus)
extern "C"
{
#endif

/* See: https://en.wikipedia.org/wiki/IEEE_754#2019 */
#define PREC_FLT_NR	9
#define PREC_DBL_NR	17
#define PREC_FLT	"9"
#define PREC_DBL	"17"

/** \brief fast_sint32
 *
 * \b Description
 *
 * Convert signed integer to string
 *
 * \param v Integer
 * \param str Buffer to print to
 * \returns lenght string
 */
  extern unsigned int fast_sint32 (int32_t v, char *str);

/** \brief fast_sint64
 *
 * \b Description
 *
 * Convert signed long integer to string
 *
 * \param v Long integer
 * \param str Buffer to print to
 * \returns lenght string
 */
  extern unsigned int fast_sint64 (int64_t v, char *str);

/** \brief fast_uint32
 *
 * \b Description
 *
 * Convert unsigned integer to string
 *
 * \param v Integer
 * \param str Buffer to print to
 * \returns lenght string
 */
  extern unsigned int fast_uint32 (uint32_t v, char *str);

/** \brief fast_uint64
 *
 * \b Description
 *
 * Convert unsigned long integer to string
 *
 * \param v Long integer
 * \param str Buffer to print to
 * \returns lenght string
 */
  extern unsigned int fast_uint64 (uint64_t v, char *str);

/** \brief fast_base_sint32
 *
 * \b Description
 *
 * Convert signed integer to string with base and upper case
 *
 * \param v Integer
 * \param str Buffer to print to
 * \param base Base to use
 * \param upper Use uppecase
 * \returns lenght string
 */
  extern unsigned int fast_base_sint32 (int32_t v, char *str, int base,
					int upper);

/** \brief fast_base_sint64
 *
 * \b Description
 *
 * Convert signed long integer to string with base and upper case
 *
 * \param v Long integer
 * \param str Buffer to print to
 * \param base Base to use
 * \param upper Use uppecase
 * \returns lenght string
 */
  extern unsigned int fast_base_sint64 (int64_t v, char *str, int base,
					int upper);

/** \brief fast_uint32
 *
 * \b Description
 *
 * Convert unsigned base_integer to string with base and upper case
 *
 * \param v Integer
 * \param str Buffer to print to
 * \param base Base to use
 * \param upper Use uppecase
 * \returns lenght string
 */
  extern unsigned int fast_base_uint32 (uint32_t v, char *str, int base,
					int upper);

/** \brief fast_base_uint64
 *
 * \b Description
 *
 * Convert unsigned long integer to string with base and upper case
 *
 * \param v Long integer
 * \param str Buffer to print to
 * \param base Base to use
 * \param upper Use uppecase
 * \returns lenght string
 */
  extern unsigned int fast_base_uint64 (uint64_t v, char *str, int base,
					int upper);

/** \brief fast_strtos32
 *
 * \b Description
 *
 * Convert string to signed integer
 *
 * \param str String to convert from
 * \param endptr optional endptr
 * \param base base of converting.
 * \returns converted string
 */
  extern int32_t fast_strtos32 (const char *str, char **endptr, int base);

/** \brief fast_strtos64
 *
 * \b Description
 *
 * Convert string to signed long
 *
 * \param str String to convert from
 * \param endptr optional endptr
 * \param base base of converting.
 * \returns converted string
 */
  extern int64_t fast_strtos64 (const char *str, char **endptr, int base);

/** \brief fast_strtou32
 *
 * \b Description
 *
 * Convert string to unsigned int
 *
 * \param str String to convert from
 * \param endptr optional endptr
 * \param base base of converting.
 * \returns converted string
 */
  extern uint32_t fast_strtou32 (const char *str, char **endptr, int base);

/** \brief fast_strtou64
 *
 * \b Description
 *
 * Convert string to unsigned long
 *
 * \param str String to convert from
 * \param endptr optional endptr
 * \param base base of converting.
 * \returns converted string
 */
  extern uint64_t fast_strtou64 (const char *str, char **endptr, int base);

/** \brief fast_ftoa
 *
 * \b Description
 *
 * Convert float to ascii
 *
 * \param v float value
 * \param size precision
 * \param line pointer to result
 * \returns lenght string
 */
  extern unsigned int fast_ftoa (float v, int size, char *line);

/** \brief fast_dtoa
 *
 * \b Description
 *
 * Convert double to ascii
 *
 * \param v float value
 * \param size precision
 * \param line pointer to result
 * \returns lenght string
 */
  extern unsigned int fast_dtoa (double v, int size, char *line);

/** \brief fast_strtof
 *
 * \b Description
 *
 * Convert string to float
 *
 * \param str string to convert
 * \param endptr optional endptr
 * \returns converted float value
 */
  extern float fast_strtof (const char *str, char **endptr);

/** \brief fast_strtod
 *
 * \b Description
 *
 * Convert string to double
 *
 * \param str string to convert
 * \param endptr optional endptr
 * \returns converted double value
 */
  extern double fast_strtod (const char *str, char **endptr);

#if defined (__cplusplus)
}
#endif

#ifndef __WORDSIZE
#define	__WORDSIZE	64
#endif

#ifdef __GNUC__
#define LIKELY(x)               __builtin_expect ((x), 1)
#define UNLIKELY(x)             __builtin_expect ((x), 0)
#else
#define LIKELY(x)               (x)
#define UNLIKELY(x)             (x)
#endif

#define	DO_BASE(b) case b: do { *--p = d[u % b]; u /= b; } while (u); break

#ifdef __attribute__
#undef __attribute__
#endif

#if 0
/* gcc -g -O3 -Wall a.c -o a -lmpfr -lgmp */
#include <stdio.h>
#include <mpfr.h>

#define	PRINT_VALUE	0

int
main (void)
{
  int i;
  int cmp;
  unsigned long value;
  unsigned long value1;
  unsigned long value2;
  long exp;
  mpfr_t p64;
  mpfr_t p32e;
  mpfr_t p32;
  mpfr_t p24;
  mpfr_t t1;
  mpfr_t t2;
  mpfr_t t3;
  mpfr_t t4;

  mpfr_init2 (p64, 128);
  mpfr_init2 (p32e, 128);
  mpfr_init2 (p32, 128);
  mpfr_init2 (p24, 128);
  mpfr_init2 (t1, 128);
  mpfr_init2 (t2, 128);
  mpfr_init2 (t3, 128);
  mpfr_init2 (t4, 128);
  mpfr_ui_pow_ui (p64, 2, 64, MPFR_RNDN);
  mpfr_ui_pow_ui (p32e, 10, 12, MPFR_RNDN);
  mpfr_ui_pow_ui (p32, 2, 32, MPFR_RNDN);
  mpfr_ui_pow_ui (p24, 2, 24, MPFR_RNDN);

  printf ("static const struct\n");
  printf ("{\n");
  printf ("  uint64_t mul1;\n");
  printf ("  uint32_t mul2;\n");
  printf ("  int32_t exp;\n");
  printf ("} dpowers2[] = {\n");
  for (i = -1022 - 53 + 64 - 11 - 54; i <= 1025 - 53 + 64 - 11; i++) {
    /* t1 = powl (2.0, i); */
    mpfr_set_ui (t3, 2, MPFR_RNDN);
    mpfr_pow_si (t1, t3, i, MPFR_RNDN);
    if (i < 0) {
      /* t2 = log10l (1.0 / t1); */
      mpfr_set_ui (t3, 1, MPFR_RNDN);
      mpfr_div (t3, t3, t1, MPFR_RNDN);
      mpfr_log10 (t2, t3, MPFR_RNDN);
      /* t2 = floorl (t2 + 22); */
      mpfr_add_ui (t3, t2, 22, MPFR_RNDN);
      mpfr_rint_floor (t2, t3, MPFR_RNDN);
      /* while (powl (10.0, t2) * t1 >= p64)  */
      mpfr_set_ui (t3, 10, MPFR_RNDN);
      mpfr_pow (t3, t3, t2, MPFR_RNDN);
      mpfr_mul (t3, t3, t1, MPFR_RNDN);
      cmp = mpfr_cmp (t3, p64);
      while (cmp >= 0) {
	/* t2 -= 1.0; */
	mpfr_sub_ui (t2, t2, 1, MPFR_RNDN);
	mpfr_set_ui (t3, 10, MPFR_RNDN);
	mpfr_pow (t3, t3, t2, MPFR_RNDN);
	mpfr_mul (t3, t3, t1, MPFR_RNDN);
	cmp = mpfr_cmp (t3, p64);
      }
      mpfr_mul (t3, t3, p32, MPFR_RNDN);
      mpfr_round (t3, t3);
      mpfr_div (t4, t3, p32, MPFR_RNDN);
      mpfr_floor (t4, t4);
      value1 = mpfr_get_ui (t4, MPFR_RNDN);
      mpfr_mul (t4, t4, p32, MPFR_RNDN);
      mpfr_sub (t4, t3, t4, MPFR_RNDN);
      value2 = mpfr_get_ui (t4, MPFR_RNDN);
      exp = mpfr_get_si (t2, MPFR_RNDN);
#if PRINT_VALUE
      mpfr_printf (" /* %40.16Rf */ ", t3);
#endif
      printf ("  { /* %5d */ UINT64_C (0x%016lX), 0x%08lX, %4ld },\n", i,
	      value1, value2, -exp);
    }
    else {
      /* t2 = log10l (t1); */
      mpfr_log10 (t2, t1, MPFR_RNDN);
      /* t2 = floorl (t2 - 22); */
      mpfr_sub_ui (t3, t2, 22, MPFR_RNDN);
      mpfr_rint_floor (t2, t3, MPFR_RNDN);
      /* while (t1 / powl (10.0, t2) >= p64) */
      mpfr_set_ui (t3, 10, MPFR_RNDN);
      mpfr_pow (t3, t3, t2, MPFR_RNDN);
      mpfr_div (t3, t1, t3, MPFR_RNDN);
      cmp = mpfr_cmp (t3, p64);
      while (cmp >= 0) {
	/* t2 += 1.0; */
	mpfr_add_ui (t2, t2, 1, MPFR_RNDN);
	mpfr_set_ui (t3, 10, MPFR_RNDN);
	mpfr_pow (t3, t3, t2, MPFR_RNDN);
	mpfr_div (t3, t1, t3, MPFR_RNDN);
	cmp = mpfr_cmp (t3, p64);
      }
      mpfr_mul (t3, t3, p32, MPFR_RNDN);
      mpfr_round (t3, t3);
      mpfr_div (t4, t3, p32, MPFR_RNDN);
      mpfr_floor (t4, t4);
      value1 = mpfr_get_ui (t4, MPFR_RNDN);
      mpfr_mul (t4, t4, p32, MPFR_RNDN);
      mpfr_sub (t4, t3, t4, MPFR_RNDN);
      value2 = mpfr_get_ui (t4, MPFR_RNDN);
      exp = mpfr_get_si (t2, MPFR_RNDN);
#if PRINT_VALUE
      mpfr_printf (" /* %40.16Rf */ ", t3);
#endif
      printf ("  { /* %5d */ UINT64_C (0x%016lX), 0x%08lX, %4ld },\n", i,
	      value1, value2, exp);
    }
  }
  printf ("};\n\n");
  printf ("static const struct\n");
  printf ("{\n");
  printf ("  uint64_t mul;\n");
  printf ("  int32_t exp;\n");
  printf ("} fpowers2[] = {\n");
  for (i = -126 - 24 + 32 - 8 - 25; i <= 129 - 24 + 32 - 8; i++) {
    /* t1 = powl (2.0, i); */
    mpfr_set_ui (t3, 2, MPFR_RNDN);
    mpfr_pow_si (t1, t3, i, MPFR_RNDN);
    if (i < 0) {
      /* t2 = log10l (1.0 / t1); */
      mpfr_set_ui (t3, 1, MPFR_RNDN);
      mpfr_div (t3, t3, t1, MPFR_RNDN);
      mpfr_log10 (t2, t3, MPFR_RNDN);
      /* t2 = floorl (t2 + 16); */
      mpfr_add_ui (t3, t2, 16, MPFR_RNDN);
      mpfr_rint_floor (t2, t3, MPFR_RNDN);
      /* while (powl (10.0, t2) * t1 >= p32e)  */
      mpfr_set_ui (t3, 10, MPFR_RNDN);
      mpfr_pow (t3, t3, t2, MPFR_RNDN);
      mpfr_mul (t3, t3, t1, MPFR_RNDN);
      cmp = mpfr_cmp (t3, p32e);
      while (cmp >= 0) {
	/* t2 -= 1.0; */
	mpfr_sub_ui (t2, t2, 1, MPFR_RNDN);
	mpfr_set_ui (t3, 10, MPFR_RNDN);
	mpfr_pow (t3, t3, t2, MPFR_RNDN);
	mpfr_mul (t3, t3, t1, MPFR_RNDN);
	cmp = mpfr_cmp (t3, p32e);
      }
      mpfr_mul (t3, t3, p24, MPFR_RNDN);
      mpfr_round (t3, t3);
      value = mpfr_get_ui (t3, MPFR_RNDN);
      exp = mpfr_get_si (t2, MPFR_RNDN);
#if PRINT_VALUE
      mpfr_printf (" /* %40.16Rf */ ", t3);
#endif
      printf ("  { /* %5d */ UINT64_C (0x%016lX), %4ld },\n", i, value, -exp);
    }
    else {
      /* t2 = log10l (t1); */
      mpfr_log10 (t2, t1, MPFR_RNDN);
      /* t2 = floorl (t2 - 16); */
      mpfr_sub_ui (t3, t2, 16, MPFR_RNDN);
      mpfr_rint_floor (t2, t3, MPFR_RNDN);
      /* while (t1 / powl (10.0, t2) >= p32e) */
      mpfr_set_ui (t3, 10, MPFR_RNDN);
      mpfr_pow (t3, t3, t2, MPFR_RNDN);
      mpfr_div (t3, t1, t3, MPFR_RNDN);
      cmp = mpfr_cmp (t3, p32e);
      while (cmp >= 0) {
	/* t2 += 1.0; */
	mpfr_add_ui (t2, t2, 1, MPFR_RNDN);
	mpfr_set_ui (t3, 10, MPFR_RNDN);
	mpfr_pow (t3, t3, t2, MPFR_RNDN);
	mpfr_div (t3, t1, t3, MPFR_RNDN);
	cmp = mpfr_cmp (t3, p32e);
      }
      mpfr_mul (t3, t3, p24, MPFR_RNDN);
      mpfr_round (t3, t3);
      value = mpfr_get_ui (t3, MPFR_RNDN);
      exp = mpfr_get_si (t2, MPFR_RNDN);
#if PRINT_VALUE
      mpfr_printf (" /* %40.16Rf */ ", t3);
#endif
      printf ("  { /* %5d */ UINT64_C (0x%016lX), %4ld },\n", i, value, exp);
    }
  }
  printf ("};\n");
  mpfr_clear (p64);
  mpfr_clear (p32e);
  mpfr_clear (p32);
  mpfr_clear (p24);
  mpfr_clear (t1);
  mpfr_clear (t2);
  mpfr_clear (t3);
  mpfr_clear (t4);
  mpfr_free_cache ();
  return 0;
}
#endif

static const struct
{
  uint64_t mul1;
  uint32_t mul2;
  int32_t exp;
} dpowers2[] = {
  { /* -1076 */ UINT64_C (0xAB69D82E364948D4), 0x6A0A4F6C, -343},
  { /* -1075 */ UINT64_C (0x22485E6FA4750E90), 0xE2020FE2, -342},
  { /* -1074 */ UINT64_C (0x4490BCDF48EA1D21), 0xC4041FC5, -342},
  { /* -1073 */ UINT64_C (0x892179BE91D43A43), 0x88083F89, -342},
  { /* -1072 */ UINT64_C (0x1B6D1859505DA540), 0xB4CE731B, -341},
  { /* -1071 */ UINT64_C (0x36DA30B2A0BB4A81), 0x699CE637, -341},
  { /* -1070 */ UINT64_C (0x6DB4616541769502), 0xD339CC6E, -341},
  { /* -1069 */ UINT64_C (0xDB68C2CA82ED2A05), 0xA67398DC, -341},
  { /* -1068 */ UINT64_C (0x2BE1C08EE6FC3B9A), 0xBAE3EB5F, -340},
  { /* -1067 */ UINT64_C (0x57C3811DCDF87735), 0x75C7D6BE, -340},
  { /* -1066 */ UINT64_C (0xAF87023B9BF0EE6A), 0xEB8FAD7C, -340},
  { /* -1065 */ UINT64_C (0x231B0072526362E2), 0x2F1CBC4C, -339},
  { /* -1064 */ UINT64_C (0x463600E4A4C6C5C4), 0x5E397898, -339},
  { /* -1063 */ UINT64_C (0x8C6C01C9498D8B88), 0xBC72F130, -339},
  { /* -1062 */ UINT64_C (0x1C1599F50EB5E8B4), 0xF27D6370, -338},
  { /* -1061 */ UINT64_C (0x382B33EA1D6BD169), 0xE4FAC6E0, -338},
  { /* -1060 */ UINT64_C (0x705667D43AD7A2D3), 0xC9F58DC0, -338},
  { /* -1059 */ UINT64_C (0xE0ACCFA875AF45A7), 0x93EB1B81, -338},
  { /* -1058 */ UINT64_C (0x2CEF5CBB4ABCA787), 0xEA6238B3, -337},
  { /* -1057 */ UINT64_C (0x59DEB97695794F0F), 0xD4C47167, -337},
  { /* -1056 */ UINT64_C (0xB3BD72ED2AF29E1F), 0xA988E2CD, -337},
  { /* -1055 */ UINT64_C (0x23F2B095D563B939), 0x884E93C3, -336},
  { /* -1054 */ UINT64_C (0x47E5612BAAC77273), 0x109D2785, -336},
  { /* -1053 */ UINT64_C (0x8FCAC257558EE4E6), 0x213A4F0B, -336},
  { /* -1052 */ UINT64_C (0x1CC226DE444FC761), 0x39D87635, -335},
  { /* -1051 */ UINT64_C (0x39844DBC889F8EC2), 0x73B0EC6B, -335},
  { /* -1050 */ UINT64_C (0x73089B79113F1D84), 0xE761D8D5, -335},
  { /* -1049 */ UINT64_C (0xE61136F2227E3B09), 0xCEC3B1AB, -335},
  { /* -1048 */ UINT64_C (0x2E037163A07FA568), 0x5C8D89EF, -334},
  { /* -1047 */ UINT64_C (0x5C06E2C740FF4AD0), 0xB91B13DE, -334},
  { /* -1046 */ UINT64_C (0xB80DC58E81FE95A1), 0x723627BC, -334},
  { /* -1045 */ UINT64_C (0x24CF8DE94D32EAB9), 0xE3A46E59, -333},
  { /* -1044 */ UINT64_C (0x499F1BD29A65D573), 0xC748DCB1, -333},
  { /* -1043 */ UINT64_C (0x933E37A534CBAAE7), 0x8E91B963, -333},
  { /* -1042 */ UINT64_C (0x1D72D7EDD75BEEFB), 0x1C838B7A, -332},
  { /* -1041 */ UINT64_C (0x3AE5AFDBAEB7DDF6), 0x390716F4, -332},
  { /* -1040 */ UINT64_C (0x75CB5FB75D6FBBEC), 0x720E2DE9, -332},
  { /* -1039 */ UINT64_C (0xEB96BF6EBADF77D8), 0xE41C5BD2, -332},
  { /* -1038 */ UINT64_C (0x2F1E2649589317F8), 0x2D9F4590, -331},
  { /* -1037 */ UINT64_C (0x5E3C4C92B1262FF0), 0x5B3E8B21, -331},
  { /* -1036 */ UINT64_C (0xBC789925624C5FE0), 0xB67D1641, -331},
  { /* -1035 */ UINT64_C (0x25B1B83AAD427993), 0x57B29E0D, -330},
  { /* -1034 */ UINT64_C (0x4B6370755A84F326), 0xAF653C1A, -330},
  { /* -1033 */ UINT64_C (0x96C6E0EAB509E64D), 0x5ECA7834, -330},
  { /* -1032 */ UINT64_C (0x1E27C69557686142), 0xAC8EE4D7, -329},
  { /* -1031 */ UINT64_C (0x3C4F8D2AAED0C285), 0x591DC9AE, -329},
  { /* -1030 */ UINT64_C (0x789F1A555DA1850A), 0xB23B935D, -329},
  { /* -1029 */ UINT64_C (0xF13E34AABB430A15), 0x647726BA, -329},
  { /* -1028 */ UINT64_C (0x303FA4222573CED1), 0x1417D48C, -328},
  { /* -1027 */ UINT64_C (0x607F48444AE79DA2), 0x282FA917, -328},
  { /* -1026 */ UINT64_C (0xC0FE908895CF3B44), 0x505F522E, -328},
  { /* -1025 */ UINT64_C (0x2699501B51297240), 0xDCDFDD3C, -327},
  { /* -1024 */ UINT64_C (0x4D32A036A252E481), 0xB9BFBA79, -327},
  { /* -1023 */ UINT64_C (0x9A65406D44A5C903), 0x737F74F2, -327},
  { /* -1022 */ UINT64_C (0x1EE10CE2A7545B67), 0x17197DCA, -326},
  { /* -1021 */ UINT64_C (0x3DC219C54EA8B6CE), 0x2E32FB94, -326},
  { /* -1020 */ UINT64_C (0x7B84338A9D516D9C), 0x5C65F728, -326},
  { /* -1019 */ UINT64_C (0xF70867153AA2DB38), 0xB8CBEE50, -326},
  { /* -1018 */ UINT64_C (0x3168149DD886F8A4), 0xF1C262DD, -325},
  { /* -1017 */ UINT64_C (0x62D0293BB10DF149), 0xE384C5BA, -325},
  { /* -1016 */ UINT64_C (0xC5A05277621BE293), 0xC7098B73, -325},
  { /* -1015 */ UINT64_C (0x278676E4AD38C6EA), 0x5B01E8B1, -324},
  { /* -1014 */ UINT64_C (0x4F0CEDC95A718DD4), 0xB603D161, -324},
  { /* -1013 */ UINT64_C (0x9E19DB92B4E31BA9), 0x6C07A2C2, -324},
  { /* -1012 */ UINT64_C (0x1F9EC583BDC70588), 0x48CE53C0, -323},
  { /* -1011 */ UINT64_C (0x3F3D8B077B8E0B10), 0x919CA781, -323},
  { /* -1010 */ UINT64_C (0x7E7B160EF71C1621), 0x23394F02, -323},
  { /* -1009 */ UINT64_C (0xFCF62C1DEE382C42), 0x46729E04, -323},
  { /* -1008 */ UINT64_C (0x3297A26C62D808DA), 0x0E16EC67, -322},
  { /* -1007 */ UINT64_C (0x652F44D8C5B011B4), 0x1C2DD8CE, -322},
  { /* -1006 */ UINT64_C (0xCA5E89B18B602368), 0x385BB19D, -322},
  { /* -1005 */ UINT64_C (0x28794EBD1BE00714), 0xD81256B9, -321},
  { /* -1004 */ UINT64_C (0x50F29D7A37C00E29), 0xB024AD72, -321},
  { /* -1003 */ UINT64_C (0xA1E53AF46F801C53), 0x60495AE4, -321},
  { /* -1002 */ UINT64_C (0x20610BCA7CB338DD), 0x79A84561, -320},
  { /* -1001 */ UINT64_C (0x40C21794F96671BA), 0xF3508AC2, -320},
  { /* -1000 */ UINT64_C (0x81842F29F2CCE375), 0xE6A11583, -320},
  { /*  -999 */ UINT64_C (0x19E73CA1FD5C2D7D), 0xFAED044D, -319},
  { /*  -998 */ UINT64_C (0x33CE7943FAB85AFB), 0xF5DA089B, -319},
  { /*  -997 */ UINT64_C (0x679CF287F570B5F7), 0xEBB41136, -319},
  { /*  -996 */ UINT64_C (0xCF39E50FEAE16BEF), 0xD768226B, -319},
  { /*  -995 */ UINT64_C (0x2971FA9CC8937BFC), 0xC4AE6D49, -318},
  { /*  -994 */ UINT64_C (0x52E3F5399126F7F9), 0x895CDA91, -318},
  { /*  -993 */ UINT64_C (0xA5C7EA73224DEFF3), 0x12B9B523, -318},
  { /*  -992 */ UINT64_C (0x2127FBB0A075FCCA), 0x36F1F107, -317},
  { /*  -991 */ UINT64_C (0x424FF76140EBF994), 0x6DE3E20E, -317},
  { /*  -990 */ UINT64_C (0x849FEEC281D7F328), 0xDBC7C41C, -317},
  { /*  -989 */ UINT64_C (0x1A8662F3B3919708), 0x2BF4C0D2, -316},
  { /*  -988 */ UINT64_C (0x350CC5E767232E10), 0x57E981A5, -316},
  { /*  -987 */ UINT64_C (0x6A198BCECE465C20), 0xAFD30349, -316},
  { /*  -986 */ UINT64_C (0xD433179D9C8CB841), 0x5FA60693, -316},
  { /*  -985 */ UINT64_C (0x2A709E52B8E8F1A6), 0xACBACE1D, -315},
  { /*  -984 */ UINT64_C (0x54E13CA571D1E34D), 0x59759C3B, -315},
  { /*  -983 */ UINT64_C (0xA9C2794AE3A3C69A), 0xB2EB3875, -315},
  { /*  -982 */ UINT64_C (0x21F3B1DBC720C152), 0x23C8A4E4, -314},
  { /*  -981 */ UINT64_C (0x43E763B78E4182A4), 0x479149C9, -314},
  { /*  -980 */ UINT64_C (0x87CEC76F1C830548), 0x8F229391, -314},
  { /*  -979 */ UINT64_C (0x1B295B1638E7010E), 0x8306EA50, -313},
  { /*  -978 */ UINT64_C (0x3652B62C71CE021D), 0x060DD4A0, -313},
  { /*  -977 */ UINT64_C (0x6CA56C58E39C043A), 0x0C1BA941, -313},
  { /*  -976 */ UINT64_C (0xD94AD8B1C7380874), 0x18375282, -313},
  { /*  -975 */ UINT64_C (0x2B755E89F4A4CE7D), 0x9E7176E7, -312},
  { /*  -974 */ UINT64_C (0x56EABD13E9499CFB), 0x3CE2EDCD, -312},
  { /*  -973 */ UINT64_C (0xADD57A27D29339F6), 0x79C5DB9B, -312},
  { /*  -972 */ UINT64_C (0x22C44BA19083D864), 0x7EC12BEC, -311},
  { /*  -971 */ UINT64_C (0x458897432107B0C8), 0xFD8257D8, -311},
  { /*  -970 */ UINT64_C (0x8B112E86420F6191), 0xFB04AFAF, -311},
  { /*  -969 */ UINT64_C (0x1BD03C81406979E9), 0xFF00EFF0, -310},
  { /*  -968 */ UINT64_C (0x37A0790280D2F3D3), 0xFE01DFE0, -310},
  { /*  -967 */ UINT64_C (0x6F40F20501A5E7A7), 0xFC03BFBF, -310},
  { /*  -966 */ UINT64_C (0xDE81E40A034BCF4F), 0xF8077F7F, -310},
  { /*  -965 */ UINT64_C (0x2C8060CECD758FDC), 0xCB34B319, -309},
  { /*  -964 */ UINT64_C (0x5900C19D9AEB1FB9), 0x96696633, -309},
  { /*  -963 */ UINT64_C (0xB201833B35D63F73), 0x2CD2CC65, -309},
  { /*  -962 */ UINT64_C (0x2399E70BD7913FE3), 0xD5C3C27B, -308},
  { /*  -961 */ UINT64_C (0x4733CE17AF227FC7), 0xAB8784F5, -308},
  { /*  -960 */ UINT64_C (0x8E679C2F5E44FF8F), 0x570F09EB, -308},
  { /*  -959 */ UINT64_C (0x1C7B1F3CAC74331C), 0xAB0301FC, -307},
  { /*  -958 */ UINT64_C (0x38F63E7958E86639), 0x560603F7, -307},
  { /*  -957 */ UINT64_C (0x71EC7CF2B1D0CC72), 0xAC0C07EF, -307},
  { /*  -956 */ UINT64_C (0xE3D8F9E563A198E5), 0x58180FDE, -307},
  { /*  -955 */ UINT64_C (0x2D91CB94472051C7), 0x7804CFF9, -306},
  { /*  -954 */ UINT64_C (0x5B2397288E40A38E), 0xF0099FF2, -306},
  { /*  -953 */ UINT64_C (0xB6472E511C81471D), 0xE0133FE5, -306},
  { /*  -952 */ UINT64_C (0x2474A2DD05B3749F), 0x93370CC7, -305},
  { /*  -951 */ UINT64_C (0x48E945BA0B66E93F), 0x266E198F, -305},
  { /*  -950 */ UINT64_C (0x91D28B7416CDD27E), 0x4CDC331D, -305},
  { /*  -949 */ UINT64_C (0x1D2A1BE4048F907F), 0xA8F8D706, -304},
  { /*  -948 */ UINT64_C (0x3A5437C8091F20FF), 0x51F1AE0C, -304},
  { /*  -947 */ UINT64_C (0x74A86F90123E41FE), 0xA3E35C17, -304},
  { /*  -946 */ UINT64_C (0xE950DF20247C83FD), 0x47C6B82F, -304},
  { /*  -945 */ UINT64_C (0x2EA9C639A0E5B3FF), 0x74C15809, -303},
  { /*  -944 */ UINT64_C (0x5D538C7341CB67FE), 0xE982B013, -303},
  { /*  -943 */ UINT64_C (0xBAA718E68396CFFD), 0xD3056026, -303},
  { /*  -942 */ UINT64_C (0x25549E9480B7C332), 0xC3CDE008, -302},
  { /*  -941 */ UINT64_C (0x4AA93D29016F8665), 0x879BC00F, -302},
  { /*  -940 */ UINT64_C (0x95527A5202DF0CCB), 0x0F37801E, -302},
  { /*  -939 */ UINT64_C (0x1DDD4BAA0093028F), 0x030B19A0, -301},
  { /*  -938 */ UINT64_C (0x3BBA97540126051E), 0x0616333F, -301},
  { /*  -937 */ UINT64_C (0x77752EA8024C0A3C), 0x0C2C667E, -301},
  { /*  -936 */ UINT64_C (0xEEEA5D5004981478), 0x1858CCFD, -301},
  { /*  -935 */ UINT64_C (0x2FC8791000EB374B), 0x3811C299, -300},
  { /*  -934 */ UINT64_C (0x5F90F22001D66E96), 0x70238532, -300},
  { /*  -933 */ UINT64_C (0xBF21E44003ACDD2C), 0xE0470A64, -300},
  { /*  -932 */ UINT64_C (0x2639FA7333EF5F6F), 0x600E3547, -299},
  { /*  -931 */ UINT64_C (0x4C73F4E667DEBEDE), 0xC01C6A8E, -299},
  { /*  -930 */ UINT64_C (0x98E7E9CCCFBD7DBD), 0x8038D51D, -299},
  { /*  -929 */ UINT64_C (0x1E94C85C298C4C59), 0x19A4F76C, -298},
  { /*  -928 */ UINT64_C (0x3D2990B8531898B2), 0x3349EED8, -298},
  { /*  -927 */ UINT64_C (0x7A532170A6313164), 0x6693DDB1, -298},
  { /*  -926 */ UINT64_C (0xF4A642E14C6262C8), 0xCD27BB61, -298},
  { /*  -925 */ UINT64_C (0x30EE0D60427A13C1), 0xC2A18BE0, -297},
  { /*  -924 */ UINT64_C (0x61DC1AC084F42783), 0x854317C0, -297},
  { /*  -923 */ UINT64_C (0xC3B8358109E84F07), 0x0A862F81, -297},
  { /*  -922 */ UINT64_C (0x2724D780352E7634), 0x9BB46FE7, -296},
  { /*  -921 */ UINT64_C (0x4E49AF006A5CEC69), 0x3768DFCD, -296},
  { /*  -920 */ UINT64_C (0x9C935E00D4B9D8D2), 0x6ED1BF9A, -296},
  { /*  -919 */ UINT64_C (0x1F50AC6690F1F82A), 0x1629F31F, -295},
  { /*  -918 */ UINT64_C (0x3EA158CD21E3F054), 0x2C53E63E, -295},
  { /*  -917 */ UINT64_C (0x7D42B19A43C7E0A8), 0x58A7CC7B, -295},
  { /*  -916 */ UINT64_C (0xFA856334878FC150), 0xB14F98F7, -295},
  { /*  -915 */ UINT64_C (0x321AAD70E7E98D10), 0x237651CB, -294},
  { /*  -914 */ UINT64_C (0x64355AE1CFD31A20), 0x46ECA396, -294},
  { /*  -913 */ UINT64_C (0xC86AB5C39FA63440), 0x8DD9472C, -294},
  { /*  -912 */ UINT64_C (0x2815578D865470D9), 0xB5F8416F, -293},
  { /*  -911 */ UINT64_C (0x502AAF1B0CA8E1B3), 0x6BF082DE, -293},
  { /*  -910 */ UINT64_C (0xA0555E361951C366), 0xD7E105BD, -293},
  { /*  -909 */ UINT64_C (0x201112D79EA9F3E1), 0x5E603459, -292},
  { /*  -908 */ UINT64_C (0x402225AF3D53E7C2), 0xBCC068B2, -292},
  { /*  -907 */ UINT64_C (0x80444B5E7AA7CF85), 0x7980D164, -292},
  { /*  -906 */ UINT64_C (0x19A742461887F64D), 0xE519C37A, -291},
  { /*  -905 */ UINT64_C (0x334E848C310FEC9B), 0xCA3386F5, -291},
  { /*  -904 */ UINT64_C (0x669D0918621FD937), 0x94670DE9, -291},
  { /*  -903 */ UINT64_C (0xCD3A1230C43FB26F), 0x28CE1BD3, -291},
  { /*  -902 */ UINT64_C (0x290B9D3CF40CBD49), 0x6E8F9F2A, -290},
  { /*  -901 */ UINT64_C (0x52173A79E8197A92), 0xDD1F3E54, -290},
  { /*  -900 */ UINT64_C (0xA42E74F3D032F525), 0xBA3E7CA9, -290},
  { /*  -899 */ UINT64_C (0x20D61763F670976D), 0xF20C7F55, -289},
  { /*  -898 */ UINT64_C (0x41AC2EC7ECE12EDB), 0xE418FEAA, -289},
  { /*  -897 */ UINT64_C (0x83585D8FD9C25DB7), 0xC831FD54, -289},
  { /*  -896 */ UINT64_C (0x1A44DF832B8D45F1), 0x8E7065DE, -288},
  { /*  -895 */ UINT64_C (0x3489BF06571A8BE3), 0x1CE0CBBB, -288},
  { /*  -894 */ UINT64_C (0x69137E0CAE3517C6), 0x39C19776, -288},
  { /*  -893 */ UINT64_C (0xD226FC195C6A2F8C), 0x73832EEC, -288},
  { /*  -892 */ UINT64_C (0x2A07CC05127BA31C), 0x171A3C96, -287},
  { /*  -891 */ UINT64_C (0x540F980A24F74638), 0x2E34792B, -287},
  { /*  -890 */ UINT64_C (0xA81F301449EE8C70), 0x5C68F257, -287},
  { /*  -889 */ UINT64_C (0x219FD66A752FB5B0), 0x127B63AB, -286},
  { /*  -888 */ UINT64_C (0x433FACD4EA5F6B60), 0x24F6C756, -286},
  { /*  -887 */ UINT64_C (0x867F59A9D4BED6C0), 0x49ED8EAC, -286},
  { /*  -886 */ UINT64_C (0x1AE64521F7595E26), 0x752F82EF, -285},
  { /*  -885 */ UINT64_C (0x35CC8A43EEB2BC4C), 0xEA5F05DE, -285},
  { /*  -884 */ UINT64_C (0x6B991487DD657899), 0xD4BE0BBD, -285},
  { /*  -883 */ UINT64_C (0xD732290FBACAF133), 0xA97C1779, -285},
  { /*  -882 */ UINT64_C (0x2B0A0836588EFD0A), 0x5518D17F, -284},
  { /*  -881 */ UINT64_C (0x5614106CB11DFA14), 0xAA31A2FD, -284},
  { /*  -880 */ UINT64_C (0xAC2820D9623BF429), 0x546345FB, -284},
  { /*  -879 */ UINT64_C (0x226E6CF846D8CA6E), 0xAA7A4132, -283},
  { /*  -878 */ UINT64_C (0x44DCD9F08DB194DD), 0x54F48264, -283},
  { /*  -877 */ UINT64_C (0x89B9B3E11B6329BA), 0xA9E904C8, -283},
  { /*  -876 */ UINT64_C (0x1B8B8A6038AD6EBE), 0xEEC83428, -282},
  { /*  -875 */ UINT64_C (0x371714C0715ADD7D), 0xDD906850, -282},
  { /*  -874 */ UINT64_C (0x6E2E2980E2B5BAFB), 0xBB20D0A0, -282},
  { /*  -873 */ UINT64_C (0xDC5C5301C56B75F7), 0x7641A141, -282},
  { /*  -872 */ UINT64_C (0x2C1277005AAF1797), 0xE47386A7, -281},
  { /*  -871 */ UINT64_C (0x5824EE00B55E2F2F), 0xC8E70D4D, -281},
  { /*  -870 */ UINT64_C (0xB049DC016ABC5E5F), 0x91CE1A9A, -281},
  { /*  -869 */ UINT64_C (0x2341F8CD1558DFAC), 0xB6C2D21F, -280},
  { /*  -868 */ UINT64_C (0x4683F19A2AB1BF59), 0x6D85A43E, -280},
  { /*  -867 */ UINT64_C (0x8D07E33455637EB2), 0xDB0B487B, -280},
  { /*  -866 */ UINT64_C (0x1C34C70A777A4C8A), 0x2BCF0E7F, -279},
  { /*  -865 */ UINT64_C (0x38698E14EEF49914), 0x579E1CFE, -279},
  { /*  -864 */ UINT64_C (0x70D31C29DDE93228), 0xAF3C39FC, -279},
  { /*  -863 */ UINT64_C (0xE1A63853BBD26451), 0x5E7873F9, -279},
  { /*  -862 */ UINT64_C (0x2D213E7725907A76), 0xAC7E7D98, -278},
  { /*  -861 */ UINT64_C (0x5A427CEE4B20F4ED), 0x58FCFB30, -278},
  { /*  -860 */ UINT64_C (0xB484F9DC9641E9DA), 0xB1F9F661, -278},
  { /*  -859 */ UINT64_C (0x241A985F514061F8), 0x89FECAE0, -277},
  { /*  -858 */ UINT64_C (0x483530BEA280C3F1), 0x13FD95C0, -277},
  { /*  -857 */ UINT64_C (0x906A617D450187E2), 0x27FB2B80, -277},
  { /*  -856 */ UINT64_C (0x1CE2137F74338193), 0xA198A24D, -276},
  { /*  -855 */ UINT64_C (0x39C426FEE8670327), 0x4331449A, -276},
  { /*  -854 */ UINT64_C (0x73884DFDD0CE064E), 0x86628934, -276},
  { /*  -853 */ UINT64_C (0xE7109BFBA19C0C9D), 0x0CC51267, -276},
  { /*  -852 */ UINT64_C (0x2E368598B9EC0285), 0xCF5A9D48, -275},
  { /*  -851 */ UINT64_C (0x5C6D0B3173D8050B), 0x9EB53A90, -275},
  { /*  -850 */ UINT64_C (0xB8DA1662E7B00A17), 0x3D6A751F, -275},
  { /*  -849 */ UINT64_C (0x24F86AE094BCCED1), 0x72AEE439, -274},
  { /*  -848 */ UINT64_C (0x49F0D5C129799DA2), 0xE55DC873, -274},
  { /*  -847 */ UINT64_C (0x93E1AB8252F33B45), 0xCABB90E6, -274},
  { /*  -846 */ UINT64_C (0x1D9388B3AA30A574), 0x5BBF1CFB, -273},
  { /*  -845 */ UINT64_C (0x3B27116754614AE8), 0xB77E39F6, -273},
  { /*  -844 */ UINT64_C (0x764E22CEA8C295D1), 0x6EFC73EB, -273},
  { /*  -843 */ UINT64_C (0xEC9C459D51852BA2), 0xDDF8E7D6, -273},
  { /*  -842 */ UINT64_C (0x2F527452A9E76F20), 0x92CB6191, -272},
  { /*  -841 */ UINT64_C (0x5EA4E8A553CEDE41), 0x2596C322, -272},
  { /*  -840 */ UINT64_C (0xBD49D14AA79DBC82), 0x4B2D8645, -272},
  { /*  -839 */ UINT64_C (0x25DB90422185F280), 0x756F8141, -271},
  { /*  -838 */ UINT64_C (0x4BB72084430BE500), 0xEADF0282, -271},
  { /*  -837 */ UINT64_C (0x976E41088617CA01), 0xD5BE0504, -271},
  { /*  -836 */ UINT64_C (0x1E494034E79E5B99), 0xF78C6767, -270},
  { /*  -835 */ UINT64_C (0x3C928069CF3CB733), 0xEF18CECE, -270},
  { /*  -834 */ UINT64_C (0x792500D39E796E67), 0xDE319D9D, -270},
  { /*  -833 */ UINT64_C (0xF24A01A73CF2DCCF), 0xBC633B39, -270},
  { /*  -832 */ UINT64_C (0x30753387D8FD5F5C), 0xBF470BD8, -269},
  { /*  -831 */ UINT64_C (0x60EA670FB1FABEB9), 0x7E8E17B1, -269},
  { /*  -830 */ UINT64_C (0xC1D4CE1F63F57D72), 0xFD1C2F61, -269},
  { /*  -829 */ UINT64_C (0x26C429397A644C4A), 0x329F3CAD, -268},
  { /*  -828 */ UINT64_C (0x4D885272F4C89894), 0x653E795A, -268},
  { /*  -827 */ UINT64_C (0x9B10A4E5E9913128), 0xCA7CF2B4, -268},
  { /*  -826 */ UINT64_C (0x1F03542DFB83703B), 0x5BB296F1, -267},
  { /*  -825 */ UINT64_C (0x3E06A85BF706E076), 0xB7652DE2, -267},
  { /*  -824 */ UINT64_C (0x7C0D50B7EE0DC0ED), 0x6ECA5BC3, -267},
  { /*  -823 */ UINT64_C (0xF81AA16FDC1B81DA), 0xDD94B787, -267},
  { /*  -822 */ UINT64_C (0x319EED165F38B392), 0x2C50F181, -266},
  { /*  -821 */ UINT64_C (0x633DDA2CBE716724), 0x58A1E303, -266},
  { /*  -820 */ UINT64_C (0xC67BB4597CE2CE48), 0xB143C605, -266},
  { /*  -819 */ UINT64_C (0x27B2574518FA2941), 0xBD0D8E01, -265},
  { /*  -818 */ UINT64_C (0x4F64AE8A31F45283), 0x7A1B1C02, -265},
  { /*  -817 */ UINT64_C (0x9EC95D1463E8A506), 0xF4363804, -265},
  { /*  -816 */ UINT64_C (0x1FC1DF6A7A61BA9A), 0xFDA4719A, -264},
  { /*  -815 */ UINT64_C (0x3F83BED4F4C37535), 0xFB48E335, -264},
  { /*  -814 */ UINT64_C (0x7F077DA9E986EA6B), 0xF691C66A, -264},
  { /*  -813 */ UINT64_C (0xFE0EFB53D30DD4D7), 0xED238CD4, -264},
  { /*  -812 */ UINT64_C (0x32CFCBDD909C5DC4), 0xC9071C2A, -263},
  { /*  -811 */ UINT64_C (0x659F97BB2138BB89), 0x920E3855, -263},
  { /*  -810 */ UINT64_C (0xCB3F2F7642717713), 0x241C70A9, -263},
  { /*  -809 */ UINT64_C (0x28A63CB1407D17D0), 0xA0D27CEF, -262},
  { /*  -808 */ UINT64_C (0x514C796280FA2FA1), 0x41A4F9DD, -262},
  { /*  -807 */ UINT64_C (0xA298F2C501F45F42), 0x8349F3BB, -262},
  { /*  -806 */ UINT64_C (0x2084FD5A99FDACA6), 0xE70ECA59, -261},
  { /*  -805 */ UINT64_C (0x4109FAB533FB594D), 0xCE1D94B1, -261},
  { /*  -804 */ UINT64_C (0x8213F56A67F6B29B), 0x9C3B2962, -261},
  { /*  -803 */ UINT64_C (0x1A03FDE214CAF085), 0x85A56EAD, -260},
  { /*  -802 */ UINT64_C (0x3407FBC42995E10B), 0x0B4ADD5A, -260},
  { /*  -801 */ UINT64_C (0x680FF788532BC216), 0x1695BAB5, -260},
  { /*  -800 */ UINT64_C (0xD01FEF10A657842C), 0x2D2B756A, -260},
  { /*  -799 */ UINT64_C (0x299FFC9CEE1180D5), 0xA2A24AAF, -259},
  { /*  -798 */ UINT64_C (0x533FF939DC2301AB), 0x4544955D, -259},
  { /*  -797 */ UINT64_C (0xA67FF273B8460356), 0x8A892ABB, -259},
  { /*  -796 */ UINT64_C (0x214CCA1724DACD77), 0xB54EA225, -258},
  { /*  -795 */ UINT64_C (0x4299942E49B59AEF), 0x6A9D444B, -258},
  { /*  -794 */ UINT64_C (0x8533285C936B35DE), 0xD53A8896, -258},
  { /*  -793 */ UINT64_C (0x1AA3D4DF50AF0AC6), 0x2AA54E84, -257},
  { /*  -792 */ UINT64_C (0x3547A9BEA15E158C), 0x554A9D09, -257},
  { /*  -791 */ UINT64_C (0x6A8F537D42BC2B18), 0xAA953A11, -257},
  { /*  -790 */ UINT64_C (0xD51EA6FA85785631), 0x552A7422, -257},
  { /*  -789 */ UINT64_C (0x2A9FBAFEE77E77A3), 0x776EE407, -256},
  { /*  -788 */ UINT64_C (0x553F75FDCEFCEF46), 0xEEDDC80E, -256},
  { /*  -787 */ UINT64_C (0xAA7EEBFB9DF9DE8D), 0xDDBB901C, -256},
  { /*  -786 */ UINT64_C (0x2219626585FEC61C), 0x5F8BE99F, -255},
  { /*  -785 */ UINT64_C (0x4432C4CB0BFD8C38), 0xBF17D33E, -255},
  { /*  -784 */ UINT64_C (0x8865899617FB1871), 0x7E2FA67C, -255},
  { /*  -783 */ UINT64_C (0x1B4781EAD1989E7D), 0x193CBAE6, -254},
  { /*  -782 */ UINT64_C (0x368F03D5A3313CFA), 0x327975CB, -254},
  { /*  -781 */ UINT64_C (0x6D1E07AB466279F4), 0x64F2EB97, -254},
  { /*  -780 */ UINT64_C (0xDA3C0F568CC4F3E8), 0xC9E5D72E, -254},
  { /*  -779 */ UINT64_C (0x2BA59CAAE8F430C8), 0x28612B09, -253},
  { /*  -778 */ UINT64_C (0x574B3955D1E86190), 0x50C25612, -253},
  { /*  -777 */ UINT64_C (0xAE9672ABA3D0C320), 0xA184AC24, -253},
  { /*  -776 */ UINT64_C (0x22EAE3BBED902706), 0x86B4226E, -252},
  { /*  -775 */ UINT64_C (0x45D5C777DB204E0D), 0x0D6844DB, -252},
  { /*  -774 */ UINT64_C (0x8BAB8EEFB6409C1A), 0x1AD089B7, -252},
  { /*  -773 */ UINT64_C (0x1BEF1C9657A6859E), 0xD229B525, -251},
  { /*  -772 */ UINT64_C (0x37DE392CAF4D0B3D), 0xA4536A49, -251},
  { /*  -771 */ UINT64_C (0x6FBC72595E9A167B), 0x48A6D492, -251},
  { /*  -770 */ UINT64_C (0xDF78E4B2BD342CF6), 0x914DA924, -251},
  { /*  -769 */ UINT64_C (0x2CB1C756F2A408FE), 0x1D0F883A, -250},
  { /*  -768 */ UINT64_C (0x59638EADE54811FC), 0x3A1F1075, -250},
  { /*  -767 */ UINT64_C (0xB2C71D5BCA9023F8), 0x743E20EA, -250},
  { /*  -766 */ UINT64_C (0x23C16C458EE9A0CB), 0x4A72D362, -249},
  { /*  -765 */ UINT64_C (0x4782D88B1DD34196), 0x94E5A6C4, -249},
  { /*  -764 */ UINT64_C (0x8F05B1163BA6832D), 0x29CB4D88, -249},
  { /*  -763 */ UINT64_C (0x1C9ABD04725480A2), 0xA1F575E8, -248},
  { /*  -762 */ UINT64_C (0x39357A08E4A90145), 0x43EAEBD0, -248},
  { /*  -761 */ UINT64_C (0x726AF411C952028A), 0x87D5D7A0, -248},
  { /*  -760 */ UINT64_C (0xE4D5E82392A40515), 0x0FABAF40, -248},
  { /*  -759 */ UINT64_C (0x2DC461A0B6ED9A9D), 0xCFEF230D, -247},
  { /*  -758 */ UINT64_C (0x5B88C3416DDB353B), 0x9FDE461A, -247},
  { /*  -757 */ UINT64_C (0xB7118682DBB66A77), 0x3FBC8C33, -247},
  { /*  -756 */ UINT64_C (0x249D1AE6F8BE154B), 0x0CBF4F3D, -246},
  { /*  -755 */ UINT64_C (0x493A35CDF17C2A96), 0x197E9E7B, -246},
  { /*  -754 */ UINT64_C (0x92746B9BE2F8552C), 0x32FD3CF6, -246},
  { /*  -753 */ UINT64_C (0x1D4A7BEBFA31AAA2), 0x70990C31, -245},
  { /*  -752 */ UINT64_C (0x3A94F7D7F4635544), 0xE1321862, -245},
  { /*  -751 */ UINT64_C (0x7529EFAFE8C6AA89), 0xC26430C5, -245},
  { /*  -750 */ UINT64_C (0xEA53DF5FD18D5513), 0x84C86189, -245},
  { /*  -749 */ UINT64_C (0x2EDD931329E91103), 0xE75B46B5, -244},
  { /*  -748 */ UINT64_C (0x5DBB262653D22207), 0xCEB68D6A, -244},
  { /*  -747 */ UINT64_C (0xBB764C4CA7A4440F), 0x9D6D1AD4, -244},
  { /*  -746 */ UINT64_C (0x257E0F4287EDA736), 0x52AF6BC4, -243},
  { /*  -745 */ UINT64_C (0x4AFC1E850FDB4E6C), 0xA55ED788, -243},
  { /*  -744 */ UINT64_C (0x95F83D0A1FB69CD9), 0x4ABDAF10, -243},
  { /*  -743 */ UINT64_C (0x1DFE729B9FF15291), 0xDBBF896A, -242},
  { /*  -742 */ UINT64_C (0x3BFCE5373FE2A523), 0xB77F12D3, -242},
  { /*  -741 */ UINT64_C (0x77F9CA6E7FC54A47), 0x6EFE25A6, -242},
  { /*  -740 */ UINT64_C (0xEFF394DCFF8A948E), 0xDDFC4B4D, -242},
  { /*  -739 */ UINT64_C (0x2FFD842C331BB74F), 0xC5FF4243, -241},
  { /*  -738 */ UINT64_C (0x5FFB085866376E9F), 0x8BFE8485, -241},
  { /*  -737 */ UINT64_C (0xBFF610B0CC6EDD3F), 0x17FD090A, -241},
  { /*  -736 */ UINT64_C (0x266469BCF5AFC5D9), 0x6B329B68, -240},
  { /*  -735 */ UINT64_C (0x4CC8D379EB5F8BB2), 0xD66536D1, -240},
  { /*  -734 */ UINT64_C (0x9991A6F3D6BF1765), 0xACCA6DA2, -240},
  { /*  -733 */ UINT64_C (0x1EB6BAFD91596B14), 0x55C215ED, -239},
  { /*  -732 */ UINT64_C (0x3D6D75FB22B2D628), 0xAB842BDA, -239},
  { /*  -731 */ UINT64_C (0x7ADAEBF64565AC51), 0x570857B5, -239},
  { /*  -730 */ UINT64_C (0xF5B5D7EC8ACB58A2), 0xAE10AF69, -239},
  { /*  -729 */ UINT64_C (0x31245E628228AB53), 0xBC69BCAF, -238},
  { /*  -728 */ UINT64_C (0x6248BCC5045156A7), 0x78D3795D, -238},
  { /*  -727 */ UINT64_C (0xC491798A08A2AD4E), 0xF1A6F2BB, -238},
  { /*  -726 */ UINT64_C (0x27504B8201BA22A9), 0x6387CA25, -237},
  { /*  -725 */ UINT64_C (0x4EA0970403744552), 0xC70F944B, -237},
  { /*  -724 */ UINT64_C (0x9D412E0806E88AA5), 0x8E1F2895, -237},
  { /*  -723 */ UINT64_C (0x1F736F9B3494E887), 0x82D3081E, -236},
  { /*  -722 */ UINT64_C (0x3EE6DF366929D10F), 0x05A6103C, -236},
  { /*  -721 */ UINT64_C (0x7DCDBE6CD253A21E), 0x0B4C2078, -236},
  { /*  -720 */ UINT64_C (0xFB9B7CD9A4A7443C), 0x169840EF, -236},
  { /*  -719 */ UINT64_C (0x32524C2B8754A73F), 0x37B80CFD, -235},
  { /*  -718 */ UINT64_C (0x64A498570EA94E7E), 0x6F7019F9, -235},
  { /*  -717 */ UINT64_C (0xC94930AE1D529CFC), 0xDEE033F2, -235},
  { /*  -716 */ UINT64_C (0x2841D689391085CC), 0x2C933D97, -234},
  { /*  -715 */ UINT64_C (0x5083AD1272210B98), 0x59267B2E, -234},
  { /*  -714 */ UINT64_C (0xA1075A24E4421730), 0xB24CF65C, -234},
  { /*  -713 */ UINT64_C (0x2034ABA0FA739E3C), 0xF075CADF, -233},
  { /*  -712 */ UINT64_C (0x40695741F4E73C79), 0xE0EB95BE, -233},
  { /*  -711 */ UINT64_C (0x80D2AE83E9CE78F3), 0xC1D72B7C, -233},
  { /*  -710 */ UINT64_C (0x19C3BC80C85C7E97), 0x26C4A24C, -232},
  { /*  -709 */ UINT64_C (0x3387790190B8FD2E), 0x4D894498, -232},
  { /*  -708 */ UINT64_C (0x670EF2032171FA5C), 0x9B128930, -232},
  { /*  -707 */ UINT64_C (0xCE1DE40642E3F4B9), 0x36251261, -232},
  { /*  -706 */ UINT64_C (0x29392D9ADA2D9758), 0x3E076A13, -231},
  { /*  -705 */ UINT64_C (0x52725B35B45B2EB0), 0x7C0ED427, -231},
  { /*  -704 */ UINT64_C (0xA4E4B66B68B65D60), 0xF81DA84D, -231},
  { /*  -703 */ UINT64_C (0x20FA8AE248247913), 0x64D2BB43, -230},
  { /*  -702 */ UINT64_C (0x41F515C49048F226), 0xC9A57685, -230},
  { /*  -701 */ UINT64_C (0x83EA2B892091E44D), 0x934AED0B, -230},
  { /*  -700 */ UINT64_C (0x1A6208B50683940F), 0x83DBC902, -229},
  { /*  -699 */ UINT64_C (0x34C4116A0D07281F), 0x07B79204, -229},
  { /*  -698 */ UINT64_C (0x698822D41A0E503E), 0x0F6F2409, -229},
  { /*  -697 */ UINT64_C (0xD31045A8341CA07C), 0x1EDE4811, -229},
  { /*  -696 */ UINT64_C (0x2A367454D738ECE5), 0x9FC60E6A, -228},
  { /*  -695 */ UINT64_C (0x546CE8A9AE71D9CB), 0x3F8C1CD4, -228},
  { /*  -694 */ UINT64_C (0xA8D9D1535CE3B396), 0x7F1839A7, -228},
  { /*  -693 */ UINT64_C (0x21C529DD78FA571E), 0x196B3EBB, -227},
  { /*  -692 */ UINT64_C (0x438A53BAF1F4AE3C), 0x32D67D76, -227},
  { /*  -691 */ UINT64_C (0x8714A775E3E95C78), 0x65ACFAEC, -227},
  { /*  -690 */ UINT64_C (0x1B04217DFA61DF4B), 0x4788FEFC, -226},
  { /*  -689 */ UINT64_C (0x360842FBF4C3BE96), 0x8F11FDF8, -226},
  { /*  -688 */ UINT64_C (0x6C1085F7E9877D2D), 0x1E23FBF0, -226},
  { /*  -687 */ UINT64_C (0xD8210BEFD30EFA5A), 0x3C47F7E0, -226},
  { /*  -686 */ UINT64_C (0x2B39CF2FF702FEDE), 0xD8DB3193, -225},
  { /*  -685 */ UINT64_C (0x56739E5FEE05FDBD), 0xB1B66327, -225},
  { /*  -684 */ UINT64_C (0xACE73CBFDC0BFB7B), 0x636CC64D, -225},
  { /*  -683 */ UINT64_C (0x2294A5BFF8CF324B), 0xE0AF5ADC, -224},
  { /*  -682 */ UINT64_C (0x45294B7FF19E6497), 0xC15EB5B8, -224},
  { /*  -681 */ UINT64_C (0x8A5296FFE33CC92F), 0x82BD6B71, -224},
  { /*  -680 */ UINT64_C (0x1BAA1E332D728EA3), 0x1A25E24A, -223},
  { /*  -679 */ UINT64_C (0x37543C665AE51D46), 0x344BC494, -223},
  { /*  -678 */ UINT64_C (0x6EA878CCB5CA3A8C), 0x68978927, -223},
  { /*  -677 */ UINT64_C (0xDD50F1996B947518), 0xD12F124E, -223},
  { /*  -676 */ UINT64_C (0x2C4363851584176B), 0x5D096A10, -222},
  { /*  -675 */ UINT64_C (0x5886C70A2B082ED6), 0xBA12D41F, -222},
  { /*  -674 */ UINT64_C (0xB10D8E1456105DAD), 0x7425A83F, -222},
  { /*  -673 */ UINT64_C (0x23691C6A779CDF89), 0x173ABB40, -221},
  { /*  -672 */ UINT64_C (0x46D238D4EF39BF12), 0x2E75767F, -221},
  { /*  -671 */ UINT64_C (0x8DA471A9DE737E24), 0x5CEAECFF, -221},
  { /*  -670 */ UINT64_C (0x1C5416BB92E3E607), 0x45C895CD, -220},
  { /*  -669 */ UINT64_C (0x38A82D7725C7CC0E), 0x8B912B99, -220},
  { /*  -668 */ UINT64_C (0x71505AEE4B8F981D), 0x17225732, -220},
  { /*  -667 */ UINT64_C (0xE2A0B5DC971F303A), 0x2E44AE65, -220},
  { /*  -666 */ UINT64_C (0x2D535792849FD672), 0x0940EFAE, -219},
  { /*  -665 */ UINT64_C (0x5AA6AF25093FACE4), 0x1281DF5B, -219},
  { /*  -664 */ UINT64_C (0xB54D5E4A127F59C8), 0x2503BEB7, -219},
  { /*  -663 */ UINT64_C (0x2442AC7536E64528), 0x07672625, -218},
  { /*  -662 */ UINT64_C (0x488558EA6DCC8A50), 0x0ECE4C49, -218},
  { /*  -661 */ UINT64_C (0x910AB1D4DB9914A0), 0x1D9C9892, -218},
  { /*  -660 */ UINT64_C (0x1D022390F8B83753), 0x391F51B7, -217},
  { /*  -659 */ UINT64_C (0x3A044721F1706EA6), 0x723EA36E, -217},
  { /*  -658 */ UINT64_C (0x74088E43E2E0DD4C), 0xE47D46DB, -217},
  { /*  -657 */ UINT64_C (0xE8111C87C5C1BA99), 0xC8FA8DB7, -217},
  { /*  -656 */ UINT64_C (0x2E69D2818DF38BB8), 0x5B654F8B, -216},
  { /*  -655 */ UINT64_C (0x5CD3A5031BE71770), 0xB6CA9F16, -216},
  { /*  -654 */ UINT64_C (0xB9A74A0637CE2EE1), 0x6D953E2C, -216},
  { /*  -653 */ UINT64_C (0x25217534718FA2F9), 0xE2B772D6, -215},
  { /*  -652 */ UINT64_C (0x4A42EA68E31F45F3), 0xC56EE5AB, -215},
  { /*  -651 */ UINT64_C (0x9485D4D1C63E8BE7), 0x8ADDCB56, -215},
  { /*  -650 */ UINT64_C (0x1DB45DC38E0C8261), 0x822C5BDE, -214},
  { /*  -649 */ UINT64_C (0x3B68BB871C1904C3), 0x0458B7BC, -214},
  { /*  -648 */ UINT64_C (0x76D1770E38320986), 0x08B16F78, -214},
  { /*  -647 */ UINT64_C (0xEDA2EE1C7064130C), 0x1162DEF0, -214},
  { /*  -646 */ UINT64_C (0x2F86FC6C167A6A35), 0x9D13C630, -213},
  { /*  -645 */ UINT64_C (0x5F0DF8D82CF4D46B), 0x3A278C60, -213},
  { /*  -644 */ UINT64_C (0xBE1BF1B059E9A8D6), 0x744F18C0, -213},
  { /*  -643 */ UINT64_C (0x260596BCDEC854F7), 0xB0DC9E8D, -212},
  { /*  -642 */ UINT64_C (0x4C0B2D79BD90A9EF), 0x61B93D1A, -212},
  { /*  -641 */ UINT64_C (0x98165AF37B2153DE), 0xC3727A33, -212},
  { /*  -640 */ UINT64_C (0x1E6ADEFD7F06AA5F), 0xC0B07ED7, -211},
  { /*  -639 */ UINT64_C (0x3CD5BDFAFE0D54BF), 0x8160FDAE, -211},
  { /*  -638 */ UINT64_C (0x79AB7BF5FC1AA97F), 0x02C1FB5C, -211},
  { /*  -637 */ UINT64_C (0xF356F7EBF83552FE), 0x0583F6B9, -211},
  { /*  -636 */ UINT64_C (0x30AAFE6264D776FF), 0x9AB3FE25, -210},
  { /*  -635 */ UINT64_C (0x6155FCC4C9AEEDFF), 0x3567FC4A, -210},
  { /*  -634 */ UINT64_C (0xC2ABF989935DDBFE), 0x6ACFF894, -210},
  { /*  -633 */ UINT64_C (0x26EF31E850AC5F32), 0xE229981E, -209},
  { /*  -632 */ UINT64_C (0x4DDE63D0A158BE65), 0xC453303B, -209},
  { /*  -631 */ UINT64_C (0x9BBCC7A142B17CCB), 0x88A66076, -209},
  { /*  -630 */ UINT64_C (0x1F25C186A6F04C28), 0xB4EE134B, -208},
  { /*  -629 */ UINT64_C (0x3E4B830D4DE09851), 0x69DC2696, -208},
  { /*  -628 */ UINT64_C (0x7C97061A9BC130A2), 0xD3B84D2B, -208},
  { /*  -627 */ UINT64_C (0xF92E0C3537826145), 0xA7709A57, -208},
  { /*  -626 */ UINT64_C (0x31D602710B1A1374), 0x54B01EDE, -207},
  { /*  -625 */ UINT64_C (0x63AC04E2163426E8), 0xA9603DBC, -207},
  { /*  -624 */ UINT64_C (0xC75809C42C684DD1), 0x52C07B79, -207},
  { /*  -623 */ UINT64_C (0x27DE685A6F480F90), 0x43C018B2, -206},
  { /*  -622 */ UINT64_C (0x4FBCD0B4DE901F20), 0x87803163, -206},
  { /*  -621 */ UINT64_C (0x9F79A169BD203E41), 0x0F0062C7, -206},
  { /*  -620 */ UINT64_C (0x1FE52048590672D9), 0xCFCCE08E, -205},
  { /*  -619 */ UINT64_C (0x3FCA4090B20CE5B3), 0x9F99C11C, -205},
  { /*  -618 */ UINT64_C (0x7F9481216419CB67), 0x3F338239, -205},
  { /*  -617 */ UINT64_C (0xFF290242C83396CE), 0x7E670471, -205},
  { /*  -616 */ UINT64_C (0x330833A6F4D71E29), 0x4C7B00E3, -204},
  { /*  -615 */ UINT64_C (0x6610674DE9AE3C52), 0x98F601C7, -204},
  { /*  -614 */ UINT64_C (0xCC20CE9BD35C78A5), 0x31EC038E, -204},
  { /*  -613 */ UINT64_C (0x28D35C8590AC1821), 0x09FC00B6, -203},
  { /*  -612 */ UINT64_C (0x51A6B90B21583042), 0x13F8016C, -203},
  { /*  -611 */ UINT64_C (0xA34D721642B06084), 0x27F002D8, -203},
  { /*  -610 */ UINT64_C (0x20A916D14089ACE7), 0x3B300092, -202},
  { /*  -609 */ UINT64_C (0x41522DA2811359CE), 0x76600123, -202},
  { /*  -608 */ UINT64_C (0x82A45B450226B39C), 0xECC00246, -202},
  { /*  -607 */ UINT64_C (0x1A20DF0DCD3AF0B8), 0xFC2666DB, -201},
  { /*  -606 */ UINT64_C (0x3441BE1B9A75E171), 0xF84CCDB6, -201},
  { /*  -605 */ UINT64_C (0x68837C3734EBC2E3), 0xF0999B6C, -201},
  { /*  -604 */ UINT64_C (0xD106F86E69D785C7), 0xE13336D7, -201},
  { /*  -603 */ UINT64_C (0x29CE31AFAEC4B45B), 0x2D0A3E2B, -200},
  { /*  -602 */ UINT64_C (0x539C635F5D8968B6), 0x5A147C56, -200},
  { /*  -601 */ UINT64_C (0xA738C6BEBB12D16C), 0xB428F8AC, -200},
  { /*  -600 */ UINT64_C (0x2171C159589D5D15), 0xBDA1CB56, -199},
  { /*  -599 */ UINT64_C (0x42E382B2B13ABA2B), 0x7B4396AB, -199},
  { /*  -598 */ UINT64_C (0x85C7056562757456), 0xF6872D56, -199},
  { /*  -597 */ UINT64_C (0x1AC1677AAD4AB0DE), 0x314E3C44, -198},
  { /*  -596 */ UINT64_C (0x3582CEF55A9561BC), 0x629C7889, -198},
  { /*  -595 */ UINT64_C (0x6B059DEAB52AC378), 0xC538F112, -198},
  { /*  -594 */ UINT64_C (0xD60B3BD56A5586F1), 0x8A71E224, -198},
  { /*  -593 */ UINT64_C (0x2ACF0BF77BAAB496), 0xB549FA07, -197},
  { /*  -592 */ UINT64_C (0x559E17EEF755692D), 0x6A93F40E, -197},
  { /*  -591 */ UINT64_C (0xAB3C2FDDEEAAD25A), 0xD527E81D, -197},
  { /*  -590 */ UINT64_C (0x223F3CC5FC889078), 0x9107FB39, -196},
  { /*  -589 */ UINT64_C (0x447E798BF91120F1), 0x220FF672, -196},
  { /*  -588 */ UINT64_C (0x88FCF317F22241E2), 0x441FECE4, -196},
  { /*  -587 */ UINT64_C (0x1B65CA37FD3A0D2D), 0x40D32F61, -195},
  { /*  -586 */ UINT64_C (0x36CB946FFA741A5A), 0x81A65EC1, -195},
  { /*  -585 */ UINT64_C (0x6D9728DFF4E834B5), 0x034CBD83, -195},
  { /*  -584 */ UINT64_C (0xDB2E51BFE9D0696A), 0x06997B06, -195},
  { /*  -583 */ UINT64_C (0x2BD610599529AEAE), 0xCE1EB234, -194},
  { /*  -582 */ UINT64_C (0x57AC20B32A535D5D), 0x9C3D6469, -194},
  { /*  -581 */ UINT64_C (0xAF58416654A6BABB), 0x387AC8D2, -194},
  { /*  -580 */ UINT64_C (0x2311A6AE10EE2558), 0xA4E55B5D, -193},
  { /*  -579 */ UINT64_C (0x46234D5C21DC4AB1), 0x49CAB6BA, -193},
  { /*  -578 */ UINT64_C (0x8C469AB843B89562), 0x93956D74, -193},
  { /*  -577 */ UINT64_C (0x1C0E1EF1A724EAAD), 0x50B77C4A, -192},
  { /*  -576 */ UINT64_C (0x381C3DE34E49D55A), 0xA16EF895, -192},
  { /*  -575 */ UINT64_C (0x70387BC69C93AAB5), 0x42DDF12A, -192},
  { /*  -574 */ UINT64_C (0xE070F78D3927556A), 0x85BBE254, -192},
  { /*  -573 */ UINT64_C (0x2CE364B5D83B1115), 0x4DF26077, -191},
  { /*  -572 */ UINT64_C (0x59C6C96BB076222A), 0x9BE4C0EE, -191},
  { /*  -571 */ UINT64_C (0xB38D92D760EC4455), 0x37C981DD, -191},
  { /*  -570 */ UINT64_C (0x23E91D5E4695A744), 0x3E5B805F, -190},
  { /*  -569 */ UINT64_C (0x47D23ABC8D2B4E88), 0x7CB700BF, -190},
  { /*  -568 */ UINT64_C (0x8FA475791A569D10), 0xF96E017D, -190},
  { /*  -567 */ UINT64_C (0x1CBA7DE5054485D0), 0x31E2CD19, -189},
  { /*  -566 */ UINT64_C (0x3974FBCA0A890BA0), 0x63C59A32, -189},
  { /*  -565 */ UINT64_C (0x72E9F79415121740), 0xC78B3464, -189},
  { /*  -564 */ UINT64_C (0xE5D3EF282A242E81), 0x8F1668C9, -189},
  { /*  -563 */ UINT64_C (0x2DF72FD4D53A6FB3), 0x83047B5B, -188},
  { /*  -562 */ UINT64_C (0x5BEE5FA9AA74DF67), 0x0608F6B7, -188},
  { /*  -561 */ UINT64_C (0xB7DCBF5354E9BECE), 0x0C11ED6D, -188},
  { /*  -560 */ UINT64_C (0x24C5BFDD7761F2F6), 0x0269FC49, -187},
  { /*  -559 */ UINT64_C (0x498B7FBAEEC3E5EC), 0x04D3F892, -187},
  { /*  -558 */ UINT64_C (0x9316FF75DD87CBD8), 0x09A7F124, -187},
  { /*  -557 */ UINT64_C (0x1D6AFFE45F818F2B), 0x352196A1, -186},
  { /*  -556 */ UINT64_C (0x3AD5FFC8BF031E56), 0x6A432D42, -186},
  { /*  -555 */ UINT64_C (0x75ABFF917E063CAC), 0xD4865A83, -186},
  { /*  -554 */ UINT64_C (0xEB57FF22FC0C7959), 0xA90CB507, -186},
  { /*  -553 */ UINT64_C (0x2F11996D659C1845), 0x21CF5768, -185},
  { /*  -552 */ UINT64_C (0x5E2332DACB38308A), 0x439EAED0, -185},
  { /*  -551 */ UINT64_C (0xBC4665B596706114), 0x873D5D9F, -185},
  { /*  -550 */ UINT64_C (0x25A7ADF11E1679D0), 0xE7D912B9, -184},
  { /*  -549 */ UINT64_C (0x4B4F5BE23C2CF3A1), 0xCFB22573, -184},
  { /*  -548 */ UINT64_C (0x969EB7C47859E743), 0x9F644AE6, -184},
  { /*  -547 */ UINT64_C (0x1E1FBE5A7E786173), 0xECADA894, -183},
  { /*  -546 */ UINT64_C (0x3C3F7CB4FCF0C2E7), 0xD95B5129, -183},
  { /*  -545 */ UINT64_C (0x787EF969F9E185CF), 0xB2B6A251, -183},
  { /*  -544 */ UINT64_C (0xF0FDF2D3F3C30B9F), 0x656D44A3, -183},
  { /*  -543 */ UINT64_C (0x3032CA2A63F3CF1F), 0xE115DA87, -182},
  { /*  -542 */ UINT64_C (0x60659454C7E79E3F), 0xC22BB50E, -182},
  { /*  -541 */ UINT64_C (0xC0CB28A98FCF3C7F), 0x84576A1C, -182},
  { /*  -540 */ UINT64_C (0x268F0821E98FD8E6), 0x4DAB1539, -181},
  { /*  -539 */ UINT64_C (0x4D1E1043D31FB1CC), 0x9B562A71, -181},
  { /*  -538 */ UINT64_C (0x9A3C2087A63F6399), 0x36AC54E3, -181},
  { /*  -537 */ UINT64_C (0x1ED8D34E547313EB), 0x7155AA94, -180},
  { /*  -536 */ UINT64_C (0x3DB1A69CA8E627D6), 0xE2AB5528, -180},
  { /*  -535 */ UINT64_C (0x7B634D3951CC4FAD), 0xC556AA4F, -180},
  { /*  -534 */ UINT64_C (0xF6C69A72A3989F5B), 0x8AAD549E, -180},
  { /*  -533 */ UINT64_C (0x315AEBB0871E8645), 0x8222AA86, -179},
  { /*  -532 */ UINT64_C (0x62B5D7610E3D0C8B), 0x0445550C, -179},
  { /*  -531 */ UINT64_C (0xC56BAEC21C7A1916), 0x088AAA18, -179},
  { /*  -530 */ UINT64_C (0x277BEFC06C186B6A), 0xCE822205, -178},
  { /*  -529 */ UINT64_C (0x4EF7DF80D830D6D5), 0x9D04440A, -178},
  { /*  -528 */ UINT64_C (0x9DEFBF01B061ADAB), 0x3A088813, -178},
  { /*  -527 */ UINT64_C (0x1F965966BCE055EF), 0x0B9B4E6A, -177},
  { /*  -526 */ UINT64_C (0x3F2CB2CD79C0ABDE), 0x17369CD5, -177},
  { /*  -525 */ UINT64_C (0x7E59659AF38157BC), 0x2E6D39A9, -177},
  { /*  -524 */ UINT64_C (0xFCB2CB35E702AF78), 0x5CDA7352, -177},
  { /*  -523 */ UINT64_C (0x328A28A46166EFE4), 0xDF5EE3DD, -176},
  { /*  -522 */ UINT64_C (0x65145148C2CDDFC9), 0xBEBDC7BB, -176},
  { /*  -521 */ UINT64_C (0xCA28A291859BBF93), 0x7D7B8F75, -176},
  { /*  -520 */ UINT64_C (0x286E86E9E7858CB7), 0x1918B64B, -175},
  { /*  -519 */ UINT64_C (0x50DD0DD3CF0B196E), 0x32316C95, -175},
  { /*  -518 */ UINT64_C (0xA1BA1BA79E1632DC), 0x6462D92A, -175},
  { /*  -517 */ UINT64_C (0x20586BEE52D13D5F), 0x4746F83C, -174},
  { /*  -516 */ UINT64_C (0x40B0D7DCA5A27ABE), 0x8E8DF077, -174},
  { /*  -515 */ UINT64_C (0x8161AFB94B44F57D), 0x1D1BE0EF, -174},
  { /*  -514 */ UINT64_C (0x19E056584240FDE5), 0xD29F2CFD, -173},
  { /*  -513 */ UINT64_C (0x33C0ACB08481FBCB), 0xA53E59F9, -173},
  { /*  -512 */ UINT64_C (0x678159610903F797), 0x4A7CB3F2, -173},
  { /*  -511 */ UINT64_C (0xCF02B2C21207EF2E), 0x94F967E4, -173},
  { /*  -510 */ UINT64_C (0x2966F08D36CE6309), 0x50FEAE61, -172},
  { /*  -509 */ UINT64_C (0x52CDE11A6D9CC612), 0xA1FD5CC2, -172},
  { /*  -508 */ UINT64_C (0xA59BC234DB398C25), 0x43FAB983, -172},
  { /*  -507 */ UINT64_C (0x211F26D75F0B826D), 0xDA65584D, -171},
  { /*  -506 */ UINT64_C (0x423E4DAEBE1704DB), 0xB4CAB09B, -171},
  { /*  -505 */ UINT64_C (0x847C9B5D7C2E09B7), 0x69956136, -171},
  { /*  -504 */ UINT64_C (0x1A7F5245E5A2CEBE), 0x48511371, -170},
  { /*  -503 */ UINT64_C (0x34FEA48BCB459D7C), 0x90A226E2, -170},
  { /*  -502 */ UINT64_C (0x69FD4917968B3AF9), 0x21444DC5, -170},
  { /*  -501 */ UINT64_C (0xD3FA922F2D1675F2), 0x42889B8A, -170},
  { /*  -500 */ UINT64_C (0x2A65506FD5D14ACA), 0x0D4E8582, -169},
  { /*  -499 */ UINT64_C (0x54CAA0DFABA29594), 0x1A9D0B04, -169},
  { /*  -498 */ UINT64_C (0xA99541BF57452B28), 0x353A1608, -169},
  { /*  -497 */ UINT64_C (0x21EAA6BFDE4108A1), 0xA43ED135, -168},
  { /*  -496 */ UINT64_C (0x43D54D7FBC821143), 0x487DA269, -168},
  { /*  -495 */ UINT64_C (0x87AA9AFF79042286), 0x90FB44D3, -168},
  { /*  -494 */ UINT64_C (0x1B221EFFE500D3B4), 0x8365742A, -167},
  { /*  -493 */ UINT64_C (0x36443DFFCA01A769), 0x06CAE854, -167},
  { /*  -492 */ UINT64_C (0x6C887BFF94034ED2), 0x0D95D0A9, -167},
  { /*  -491 */ UINT64_C (0xD910F7FF28069DA4), 0x1B2BA152, -167},
  { /*  -490 */ UINT64_C (0x2B69CB33080152BA), 0x6BD586AA, -166},
  { /*  -489 */ UINT64_C (0x56D396661002A574), 0xD7AB0D54, -166},
  { /*  -488 */ UINT64_C (0xADA72CCC20054AE9), 0xAF561AA8, -166},
  { /*  -487 */ UINT64_C (0x22BB08F5A0010EFB), 0x89779EEE, -165},
  { /*  -486 */ UINT64_C (0x457611EB40021DF7), 0x12EF3DDD, -165},
  { /*  -485 */ UINT64_C (0x8AEC23D680043BEE), 0x25DE7BB9, -165},
  { /*  -484 */ UINT64_C (0x1BC8D3F7B3340BFC), 0x6DF94BF2, -164},
  { /*  -483 */ UINT64_C (0x3791A7EF666817F8), 0xDBF297E4, -164},
  { /*  -482 */ UINT64_C (0x6F234FDECCD02FF1), 0xB7E52FC7, -164},
  { /*  -481 */ UINT64_C (0xDE469FBD99A05FE3), 0x6FCA5F8F, -164},
  { /*  -480 */ UINT64_C (0x2C7486591EB9ACC7), 0x165BACB6, -163},
  { /*  -479 */ UINT64_C (0x58E90CB23D73598E), 0x2CB7596C, -163},
  { /*  -478 */ UINT64_C (0xB1D219647AE6B31C), 0x596EB2D9, -163},
  { /*  -477 */ UINT64_C (0x23906B7A7EFAF09F), 0x451623C5, -162},
  { /*  -476 */ UINT64_C (0x4720D6F4FDF5E13E), 0x8A2C478A, -162},
  { /*  -475 */ UINT64_C (0x8E41ADE9FBEBC27D), 0x14588F14, -162},
  { /*  -474 */ UINT64_C (0x1C73892ECBFBF3B2), 0x9DAB4FD1, -161},
  { /*  -473 */ UINT64_C (0x38E7125D97F7E765), 0x3B569FA1, -161},
  { /*  -472 */ UINT64_C (0x71CE24BB2FEFCECA), 0x76AD3F43, -161},
  { /*  -471 */ UINT64_C (0xE39C49765FDF9D94), 0xED5A7E86, -161},
  { /*  -470 */ UINT64_C (0x2D85A84ADFF985EA), 0x95DEE61B, -160},
  { /*  -469 */ UINT64_C (0x5B0B5095BFF30BD5), 0x2BBDCC36, -160},
  { /*  -468 */ UINT64_C (0xB616A12B7FE617AA), 0x577B986B, -160},
  { /*  -467 */ UINT64_C (0x246AED08B32E04BB), 0xAB18B815, -159},
  { /*  -466 */ UINT64_C (0x48D5DA11665C0977), 0x5631702B, -159},
  { /*  -465 */ UINT64_C (0x91ABB422CCB812EE), 0xAC62E056, -159},
  { /*  -464 */ UINT64_C (0x1D22573A28F19D62), 0xEF46F9AB, -158},
  { /*  -463 */ UINT64_C (0x3A44AE7451E33AC5), 0xDE8DF356, -158},
  { /*  -462 */ UINT64_C (0x74895CE8A3C6758B), 0xBD1BE6AB, -158},
  { /*  -461 */ UINT64_C (0xE912B9D1478CEB17), 0x7A37CD56, -158},
  { /*  -460 */ UINT64_C (0x2E9D585D0E4F6237), 0xE53E5C44, -157},
  { /*  -459 */ UINT64_C (0x5D3AB0BA1C9EC46F), 0xCA7CB889, -157},
  { /*  -458 */ UINT64_C (0xBA756174393D88DF), 0x94F97112, -157},
  { /*  -457 */ UINT64_C (0x254AAD173EA5E82C), 0xB765169D, -156},
  { /*  -456 */ UINT64_C (0x4A955A2E7D4BD059), 0x6ECA2D3A, -156},
  { /*  -455 */ UINT64_C (0x952AB45CFA97A0B2), 0xDD945A74, -156},
  { /*  -454 */ UINT64_C (0x1DD55745CBB7ECF0), 0x92B7454A, -155},
  { /*  -453 */ UINT64_C (0x3BAAAE8B976FD9E1), 0x256E8A95, -155},
  { /*  -452 */ UINT64_C (0x77555D172EDFB3C2), 0x4ADD152A, -155},
  { /*  -451 */ UINT64_C (0xEEAABA2E5DBF6784), 0x95BA2A54, -155},
  { /*  -450 */ UINT64_C (0x2FBBBED612BFE180), 0xEABED544, -154},
  { /*  -449 */ UINT64_C (0x5F777DAC257FC301), 0xD57DAA88, -154},
  { /*  -448 */ UINT64_C (0xBEEEFB584AFF8603), 0xAAFB5510, -154},
  { /*  -447 */ UINT64_C (0x262FCBDE75664E00), 0xBBCBDDD0, -153},
  { /*  -446 */ UINT64_C (0x4C5F97BCEACC9C01), 0x7797BBA0, -153},
  { /*  -445 */ UINT64_C (0x98BF2F79D5993802), 0xEF2F7740, -153},
  { /*  -444 */ UINT64_C (0x1E8CA3185DEB719A), 0x2FD64B0D, -152},
  { /*  -443 */ UINT64_C (0x3D194630BBD6E334), 0x5FAC961A, -152},
  { /*  -442 */ UINT64_C (0x7A328C6177ADC668), 0xBF592C33, -152},
  { /*  -441 */ UINT64_C (0xF46518C2EF5B8CD1), 0x7EB25866, -152},
  { /*  -440 */ UINT64_C (0x30E104F3C978B5C3), 0x7FBD44E1, -151},
  { /*  -439 */ UINT64_C (0x61C209E792F16B86), 0xFF7A89C3, -151},
  { /*  -438 */ UINT64_C (0xC38413CF25E2D70D), 0xFEF51385, -151},
  { /*  -437 */ UINT64_C (0x271A6A5CA12D5E35), 0xFFCA9D81, -150},
  { /*  -436 */ UINT64_C (0x4E34D4B9425ABC6B), 0xFF953B02, -150},
  { /*  -435 */ UINT64_C (0x9C69A97284B578D7), 0xFF2A7604, -150},
  { /*  -434 */ UINT64_C (0x1F485516E7577E91), 0x996EE467, -149},
  { /*  -433 */ UINT64_C (0x3E90AA2DCEAEFD23), 0x32DDC8CE, -149},
  { /*  -432 */ UINT64_C (0x7D21545B9D5DFA46), 0x65BB919D, -149},
  { /*  -431 */ UINT64_C (0xFA42A8B73ABBF48C), 0xCB77233A, -149},
  { /*  -430 */ UINT64_C (0x320D54F17225974F), 0x5BE4A0A5, -148},
  { /*  -429 */ UINT64_C (0x641AA9E2E44B2E9E), 0xB7C9414A, -148},
  { /*  -428 */ UINT64_C (0xC83553C5C8965D3D), 0x6F928295, -148},
  { /*  -427 */ UINT64_C (0x280AAA5AC1B7AC3F), 0x7CB6E6EB, -147},
  { /*  -426 */ UINT64_C (0x501554B5836F587E), 0xF96DCDD5, -147},
  { /*  -425 */ UINT64_C (0xA02AA96B06DEB0FD), 0xF2DB9BAA, -147},
  { /*  -424 */ UINT64_C (0x200888489AF95699), 0x30925255, -146},
  { /*  -423 */ UINT64_C (0x4011109135F2AD32), 0x6124A4AA, -146},
  { /*  -422 */ UINT64_C (0x802221226BE55A64), 0xC2494955, -146},
  { /*  -421 */ UINT64_C (0x19A06D06E2611214), 0x26DB7511, -145},
  { /*  -420 */ UINT64_C (0x3340DA0DC4C22428), 0x4DB6EA22, -145},
  { /*  -419 */ UINT64_C (0x6681B41B89844850), 0x9B6DD444, -145},
  { /*  -418 */ UINT64_C (0xCD036837130890A1), 0x36DBA888, -145},
  { /*  -417 */ UINT64_C (0x2900AE716A34E9B9), 0xD7C5881B, -144},
  { /*  -416 */ UINT64_C (0x52015CE2D469D373), 0xAF8B1036, -144},
  { /*  -415 */ UINT64_C (0xA402B9C5A8D3A6E7), 0x5F16206D, -144},
  { /*  -414 */ UINT64_C (0x20CD585ABB5D87C7), 0xDFD139AF, -143},
  { /*  -413 */ UINT64_C (0x419AB0B576BB0F8F), 0xBFA2735F, -143},
  { /*  -412 */ UINT64_C (0x8335616AED761F1F), 0x7F44E6BD, -143},
  { /*  -411 */ UINT64_C (0x1A3DE04895E46C9F), 0xE640FAF3, -142},
  { /*  -410 */ UINT64_C (0x347BC0912BC8D93F), 0xCC81F5E5, -142},
  { /*  -409 */ UINT64_C (0x68F781225791B27F), 0x9903EBCB, -142},
  { /*  -408 */ UINT64_C (0xD1EF0244AF2364FF), 0x3207D795, -142},
  { /*  -407 */ UINT64_C (0x29FC9A0DBCA0ADCC), 0xA39B2B1E, -141},
  { /*  -406 */ UINT64_C (0x53F9341B79415B99), 0x4736563C, -141},
  { /*  -405 */ UINT64_C (0xA7F26836F282B732), 0x8E6CAC77, -141},
  { /*  -404 */ UINT64_C (0x2196E1A496E6F170), 0x82E288E5, -140},
  { /*  -403 */ UINT64_C (0x432DC3492DCDE2E1), 0x05C511C9, -140},
  { /*  -402 */ UINT64_C (0x865B86925B9BC5C2), 0x0B8A2393, -140},
  { /*  -401 */ UINT64_C (0x1ADF1AEA12525AC0), 0x68B53A51, -139},
  { /*  -400 */ UINT64_C (0x35BE35D424A4B580), 0xD16A74A1, -139},
  { /*  -399 */ UINT64_C (0x6B7C6BA849496B01), 0xA2D4E942, -139},
  { /*  -398 */ UINT64_C (0xD6F8D7509292D603), 0x45A9D284, -139},
  { /*  -397 */ UINT64_C (0x2AFE917683B6F79A), 0x4121F6E7, -138},
  { /*  -396 */ UINT64_C (0x55FD22ED076DEF34), 0x8243EDCF, -138},
  { /*  -395 */ UINT64_C (0xABFA45DA0EDBDE69), 0x0487DB9D, -138},
  { /*  -394 */ UINT64_C (0x2265412B9C925FAE), 0x9A819253, -137},
  { /*  -393 */ UINT64_C (0x44CA82573924BF5D), 0x350324A5, -137},
  { /*  -392 */ UINT64_C (0x899504AE72497EBA), 0x6A06494A, -137},
  { /*  -391 */ UINT64_C (0x1B843422E3A84C8B), 0xAECE0EA8, -136},
  { /*  -390 */ UINT64_C (0x37086845C7509917), 0x5D9C1D51, -136},
  { /*  -389 */ UINT64_C (0x6E10D08B8EA1322E), 0xBB383AA2, -136},
  { /*  -388 */ UINT64_C (0xDC21A1171D42645D), 0x76707544, -136},
  { /*  -387 */ UINT64_C (0x2C06B9D16C407A79), 0x17B01774, -135},
  { /*  -386 */ UINT64_C (0x580D73A2D880F4F2), 0x2F602EE8, -135},
  { /*  -385 */ UINT64_C (0xB01AE745B101E9E4), 0x5EC05DD0, -135},
  { /*  -384 */ UINT64_C (0x233894A789CD2EC7), 0x4626792A, -134},
  { /*  -383 */ UINT64_C (0x4671294F139A5D8E), 0x8C4CF253, -134},
  { /*  -382 */ UINT64_C (0x8CE2529E2734BB1D), 0x1899E4A6, -134},
  { /*  -381 */ UINT64_C (0x1C2D43B93B0A8BD2), 0x9E852DBB, -133},
  { /*  -380 */ UINT64_C (0x385A8772761517A5), 0x3D0A5B76, -133},
  { /*  -379 */ UINT64_C (0x70B50EE4EC2A2F4A), 0x7A14B6EB, -133},
  { /*  -378 */ UINT64_C (0xE16A1DC9D8545E94), 0xF4296DD7, -133},
  { /*  -377 */ UINT64_C (0x2D1539285E77461D), 0xCA6EAF91, -132},
  { /*  -376 */ UINT64_C (0x5A2A7250BCEE8C3B), 0x94DD5F23, -132},
  { /*  -375 */ UINT64_C (0xB454E4A179DD1877), 0x29BABE46, -132},
  { /*  -374 */ UINT64_C (0x2410FA86B1F904E4), 0xA1F2260E, -131},
  { /*  -373 */ UINT64_C (0x4821F50D63F209C9), 0x43E44C1C, -131},
  { /*  -372 */ UINT64_C (0x9043EA1AC7E41392), 0x87C89838, -131},
  { /*  -371 */ UINT64_C (0x1CDA62055B2D9D83), 0xB4C1B80B, -130},
  { /*  -370 */ UINT64_C (0x39B4C40AB65B3B07), 0x69837016, -130},
  { /*  -369 */ UINT64_C (0x736988156CB6760E), 0xD306E02D, -130},
  { /*  -368 */ UINT64_C (0xE6D3102AD96CEC1D), 0xA60DC059, -130},
  { /*  -367 */ UINT64_C (0x2E2A366EF848FC05), 0xEE02C012, -129},
  { /*  -366 */ UINT64_C (0x5C546CDDF091F80B), 0xDC058024, -129},
  { /*  -365 */ UINT64_C (0xB8A8D9BBE123F017), 0xB80B0047, -129},
  { /*  -364 */ UINT64_C (0x24EE91F2603A6337), 0xF19BCCDB, -128},
  { /*  -363 */ UINT64_C (0x49DD23E4C074C66F), 0xE33799B6, -128},
  { /*  -362 */ UINT64_C (0x93BA47C980E98CDF), 0xC66F336C, -128},
  { /*  -361 */ UINT64_C (0x1D8BA7F519C84F5F), 0xF47CA3E2, -127},
  { /*  -360 */ UINT64_C (0x3B174FEA33909EBF), 0xE8F947C5, -127},
  { /*  -359 */ UINT64_C (0x762E9FD467213D7F), 0xD1F28F8A, -127},
  { /*  -358 */ UINT64_C (0xEC5D3FA8CE427AFF), 0xA3E51F14, -127},
  { /*  -357 */ UINT64_C (0x2F45D98829407EFF), 0xED94396A, -126},
  { /*  -356 */ UINT64_C (0x5E8BB3105280FDFF), 0xDB2872D5, -126},
  { /*  -355 */ UINT64_C (0xBD176620A501FBFF), 0xB650E5A9, -126},
  { /*  -354 */ UINT64_C (0x25D17AD3543398CC), 0xBE102DEF, -125},
  { /*  -353 */ UINT64_C (0x4BA2F5A6A8673199), 0x7C205BDD, -125},
  { /*  -352 */ UINT64_C (0x9745EB4D50CE6332), 0xF840B7BB, -125},
  { /*  -351 */ UINT64_C (0x1E412F0F768FAD70), 0x980CF18C, -124},
  { /*  -350 */ UINT64_C (0x3C825E1EED1F5AE1), 0x3019E317, -124},
  { /*  -349 */ UINT64_C (0x7904BC3DDA3EB5C2), 0x6033C62F, -124},
  { /*  -348 */ UINT64_C (0xF209787BB47D6B84), 0xC0678C5E, -124},
  { /*  -347 */ UINT64_C (0x30684B4BF0E5E24D), 0xC014B5AC, -123},
  { /*  -346 */ UINT64_C (0x60D09697E1CBC49B), 0x80296B59, -123},
  { /*  -345 */ UINT64_C (0xC1A12D2FC3978937), 0x0052D6B1, -123},
  { /*  -344 */ UINT64_C (0x26B9D5D65A5181D7), 0xCCDD5E23, -122},
  { /*  -343 */ UINT64_C (0x4D73ABACB4A303AF), 0x99BABC47, -122},
  { /*  -342 */ UINT64_C (0x9AE757596946075F), 0x3375788E, -122},
  { /*  -341 */ UINT64_C (0x1EFB1178484134AC), 0xA3E44B50, -121},
  { /*  -340 */ UINT64_C (0x3DF622F090826959), 0x47C8969F, -121},
  { /*  -339 */ UINT64_C (0x7BEC45E12104D2B2), 0x8F912D3E, -121},
  { /*  -338 */ UINT64_C (0xF7D88BC24209A565), 0x1F225A7D, -121},
  { /*  -337 */ UINT64_C (0x3191B58D40685447), 0x6CA0787F, -120},
  { /*  -336 */ UINT64_C (0x63236B1A80D0A88E), 0xD940F0FF, -120},
  { /*  -335 */ UINT64_C (0xC646D63501A1511D), 0xB281E1FD, -120},
  { /*  -334 */ UINT64_C (0x27A7C4710053769F), 0x8A19F9FF, -119},
  { /*  -333 */ UINT64_C (0x4F4F88E200A6ED3F), 0x1433F3FF, -119},
  { /*  -332 */ UINT64_C (0x9E9F11C4014DDA7E), 0x2867E7FE, -119},
  { /*  -331 */ UINT64_C (0x1FB969F40042C54C), 0x6E7B2E66, -118},
  { /*  -330 */ UINT64_C (0x3F72D3E800858A98), 0xDCF65CCC, -118},
  { /*  -329 */ UINT64_C (0x7EE5A7D0010B1531), 0xB9ECB998, -118},
  { /*  -328 */ UINT64_C (0xFDCB4FA002162A63), 0x73D97330, -118},
  { /*  -327 */ UINT64_C (0x32C24320006AD547), 0x172B7D70, -117},
  { /*  -326 */ UINT64_C (0x6584864000D5AA8E), 0x2E56FAE0, -117},
  { /*  -325 */ UINT64_C (0xCB090C8001AB551C), 0x5CADF5C0, -117},
  { /*  -324 */ UINT64_C (0x289B68E666BBDDD2), 0x78EF978D, -116},
  { /*  -323 */ UINT64_C (0x5136D1CCCD77BBA4), 0xF1DF2F1A, -116},
  { /*  -322 */ UINT64_C (0xA26DA3999AEF7749), 0xE3BE5E33, -116},
  { /*  -321 */ UINT64_C (0x207C53EB856317DB), 0x93F2DFA4, -115},
  { /*  -320 */ UINT64_C (0x40F8A7D70AC62FB7), 0x27E5BF48, -115},
  { /*  -319 */ UINT64_C (0x81F14FAE158C5F6E), 0x4FCB7E8F, -115},
  { /*  -318 */ UINT64_C (0x19FD0FEF9DE8DFE2), 0xDCC24C83, -114},
  { /*  -317 */ UINT64_C (0x33FA1FDF3BD1BFC5), 0xB9849906, -114},
  { /*  -316 */ UINT64_C (0x67F43FBE77A37F8B), 0x7309320C, -114},
  { /*  -315 */ UINT64_C (0xCFE87F7CEF46FF16), 0xE6126418, -114},
  { /*  -314 */ UINT64_C (0x2994E64C2FDAFFD1), 0x6136E0D2, -113},
  { /*  -313 */ UINT64_C (0x5329CC985FB5FFA2), 0xC26DC1A3, -113},
  { /*  -312 */ UINT64_C (0xA6539930BF6BFF45), 0x84DB8347, -113},
  { /*  -311 */ UINT64_C (0x2143EB702648CCA7), 0x80F8B3DB, -112},
  { /*  -310 */ UINT64_C (0x4287D6E04C91994F), 0x01F167B6, -112},
  { /*  -309 */ UINT64_C (0x850FADC09923329E), 0x03E2CF6C, -112},
  { /*  -308 */ UINT64_C (0x1A9CBC59B83A3D52), 0xCD93C316, -111},
  { /*  -307 */ UINT64_C (0x353978B370747AA5), 0x9B27862B, -111},
  { /*  -306 */ UINT64_C (0x6A72F166E0E8F54B), 0x364F0C56, -111},
  { /*  -305 */ UINT64_C (0xD4E5E2CDC1D1EA96), 0x6C9E18AC, -111},
  { /*  -304 */ UINT64_C (0x2A94608F8D29FBB7), 0xAF52D1BC, -110},
  { /*  -303 */ UINT64_C (0x5528C11F1A53F76F), 0x5EA5A378, -110},
  { /*  -302 */ UINT64_C (0xAA51823E34A7EEDE), 0xBD4B46F0, -110},
  { /*  -301 */ UINT64_C (0x22104D3FA421962C), 0x8C424163, -109},
  { /*  -300 */ UINT64_C (0x44209A7F48432C59), 0x188482C7, -109},
  { /*  -299 */ UINT64_C (0x884134FE908658B2), 0x3109058D, -109},
  { /*  -298 */ UINT64_C (0x1B403DCC834E11BD), 0x3D01CDE9, -108},
  { /*  -297 */ UINT64_C (0x36807B99069C237A), 0x7A039BD2, -108},
  { /*  -296 */ UINT64_C (0x6D00F7320D3846F4), 0xF40737A4, -108},
  { /*  -295 */ UINT64_C (0xDA01EE641A708DE9), 0xE80E6F48, -108},
  { /*  -294 */ UINT64_C (0x2B99FC7A6BB01C61), 0xFB361642, -107},
  { /*  -293 */ UINT64_C (0x5733F8F4D76038C3), 0xF66C2C83, -107},
  { /*  -292 */ UINT64_C (0xAE67F1E9AEC07187), 0xECD85907, -107},
  { /*  -291 */ UINT64_C (0x22E196C856267D1B), 0x2F5E7835, -106},
  { /*  -290 */ UINT64_C (0x45C32D90AC4CFA36), 0x5EBCF069, -106},
  { /*  -289 */ UINT64_C (0x8B865B215899F46C), 0xBD79E0D2, -106},
  { /*  -288 */ UINT64_C (0x1BE7ABD3781ECA7C), 0x25E52CF7, -105},
  { /*  -287 */ UINT64_C (0x37CF57A6F03D94F8), 0x4BCA59EE, -105},
  { /*  -286 */ UINT64_C (0x6F9EAF4DE07B29F0), 0x9794B3DB, -105},
  { /*  -285 */ UINT64_C (0xDF3D5E9BC0F653E1), 0x2F2967B6, -105},
  { /*  -284 */ UINT64_C (0x2CA5DFB8C03143F9), 0xD63B7B24, -104},
  { /*  -283 */ UINT64_C (0x594BBF71806287F3), 0xAC76F649, -104},
  { /*  -282 */ UINT64_C (0xB2977EE300C50FE7), 0x58EDEC92, -104},
  { /*  -281 */ UINT64_C (0x23B7E62D668DCFFB), 0x11C92F50, -103},
  { /*  -280 */ UINT64_C (0x476FCC5ACD1B9FF6), 0x23925EA1, -103},
  { /*  -279 */ UINT64_C (0x8EDF98B59A373FEC), 0x4724BD42, -103},
  { /*  -278 */ UINT64_C (0x1C931E8AB871732F), 0x416DBF74, -102},
  { /*  -277 */ UINT64_C (0x39263D1570E2E65E), 0x82DB7EE7, -102},
  { /*  -276 */ UINT64_C (0x724C7A2AE1C5CCBD), 0x05B6FDCE, -102},
  { /*  -275 */ UINT64_C (0xE498F455C38B997A), 0x0B6DFB9C, -102},
  { /*  -274 */ UINT64_C (0x2DB830DDF3E8B84B), 0x9BE2CBEC, -101},
  { /*  -273 */ UINT64_C (0x5B7061BBE7D17097), 0x37C597D8, -101},
  { /*  -272 */ UINT64_C (0xB6E0C377CFA2E12E), 0x6F8B2FB0, -101},
  { /*  -271 */ UINT64_C (0x24935A4B2986F9D6), 0x164F098A, -100},
  { /*  -270 */ UINT64_C (0x4926B496530DF3AC), 0x2C9E1313, -100},
  { /*  -269 */ UINT64_C (0x924D692CA61BE758), 0x593C2626, -100},
  { /*  -268 */ UINT64_C (0x1D42AEA2879F2E44), 0xDEA5A13B, -99},
  { /*  -267 */ UINT64_C (0x3A855D450F3E5C89), 0xBD4B4276, -99},
  { /*  -266 */ UINT64_C (0x750ABA8A1E7CB913), 0x7A9684EC, -99},
  { /*  -265 */ UINT64_C (0xEA1575143CF97226), 0xF52D09D7, -99},
  { /*  -264 */ UINT64_C (0x2ED1176A72984A07), 0xCAA29B91, -98},
  { /*  -263 */ UINT64_C (0x5DA22ED4E530940F), 0x95453723, -98},
  { /*  -262 */ UINT64_C (0xBB445DA9CA61281F), 0x2A8A6E46, -98},
  { /*  -261 */ UINT64_C (0x257412BB8EE03B39), 0x6EE87C74, -97},
  { /*  -260 */ UINT64_C (0x4AE825771DC07672), 0xDDD0F8E9, -97},
  { /*  -259 */ UINT64_C (0x95D04AEE3B80ECE5), 0xBBA1F1D1, -97},
  { /*  -258 */ UINT64_C (0x1DF67562D8B36294), 0x58B9FD2A, -96},
  { /*  -257 */ UINT64_C (0x3BECEAC5B166C528), 0xB173FA54, -96},
  { /*  -256 */ UINT64_C (0x77D9D58B62CD8A51), 0x62E7F4A7, -96},
  { /*  -255 */ UINT64_C (0xEFB3AB16C59B14A2), 0xC5CFE94F, -96},
  { /*  -254 */ UINT64_C (0x2FF0BBD15AB89DBA), 0x278FFB76, -95},
  { /*  -253 */ UINT64_C (0x5FE177A2B5713B74), 0x4F1FF6EC, -95},
  { /*  -252 */ UINT64_C (0xBFC2EF456AE276E8), 0x9E3FEDD9, -95},
  { /*  -251 */ UINT64_C (0x265A2FDAAEFA17C8), 0x1FA662C5, -94},
  { /*  -250 */ UINT64_C (0x4CB45FB55DF42F90), 0x3F4CC58A, -94},
  { /*  -249 */ UINT64_C (0x9968BF6ABBE85F20), 0x7E998B14, -94},
  { /*  -248 */ UINT64_C (0x1EAE8CAEF261ACA0), 0x1951E89E, -93},
  { /*  -247 */ UINT64_C (0x3D5D195DE4C35940), 0x32A3D13B, -93},
  { /*  -246 */ UINT64_C (0x7ABA32BBC986B280), 0x6547A276, -93},
  { /*  -245 */ UINT64_C (0xF5746577930D6500), 0xCA8F44EC, -93},
  { /*  -244 */ UINT64_C (0x3117477E509C4766), 0x8EE9742F, -92},
  { /*  -243 */ UINT64_C (0x622E8EFCA1388ECD), 0x1DD2E85F, -92},
  { /*  -242 */ UINT64_C (0xC45D1DF942711D9A), 0x3BA5D0BD, -92},
  { /*  -241 */ UINT64_C (0x2745D2CB73B0391E), 0xD8BAC359, -91},
  { /*  -240 */ UINT64_C (0x4E8BA596E760723D), 0xB17586B2, -91},
  { /*  -239 */ UINT64_C (0x9D174B2DCEC0E47B), 0x62EB0D64, -91},
  { /*  -238 */ UINT64_C (0x1F6B0F092959C74B), 0xE0956914, -90},
  { /*  -237 */ UINT64_C (0x3ED61E1252B38E97), 0xC12AD228, -90},
  { /*  -236 */ UINT64_C (0x7DAC3C24A5671D2F), 0x8255A450, -90},
  { /*  -235 */ UINT64_C (0xFB5878494ACE3A5F), 0x04AB48A0, -90},
  { /*  -234 */ UINT64_C (0x3244E4DB755C7213), 0x00EF0E86, -89},
  { /*  -233 */ UINT64_C (0x6489C9B6EAB8E426), 0x01DE1D0D, -89},
  { /*  -232 */ UINT64_C (0xC913936DD571C84C), 0x03BC3A1A, -89},
  { /*  -231 */ UINT64_C (0x28371D7C5DE38E75), 0x9A58D86C, -88},
  { /*  -230 */ UINT64_C (0x506E3AF8BBC71CEB), 0x34B1B0D7, -88},
  { /*  -229 */ UINT64_C (0xA0DC75F1778E39D6), 0x696361AE, -88},
  { /*  -228 */ UINT64_C (0x202C1796B182D85E), 0x1513E056, -87},
  { /*  -227 */ UINT64_C (0x40582F2D6305B0BC), 0x2A27C0AC, -87},
  { /*  -226 */ UINT64_C (0x80B05E5AC60B6178), 0x544F8158, -87},
  { /*  -225 */ UINT64_C (0x19BCDFABC13579E4), 0xDDA98045, -86},
  { /*  -224 */ UINT64_C (0x3379BF57826AF3C9), 0xBB53008A, -86},
  { /*  -223 */ UINT64_C (0x66F37EAF04D5E793), 0x76A60113, -86},
  { /*  -222 */ UINT64_C (0xCDE6FD5E09ABCF26), 0xED4C0227, -86},
  { /*  -221 */ UINT64_C (0x292E32AC68558FD4), 0x95DC006E, -85},
  { /*  -220 */ UINT64_C (0x525C6558D0AB1FA9), 0x2BB800DC, -85},
  { /*  -219 */ UINT64_C (0xA4B8CAB1A1563F52), 0x577001B9, -85},
  { /*  -218 */ UINT64_C (0x20F1C22386AAD976), 0xDE4999F2, -84},
  { /*  -217 */ UINT64_C (0x41E384470D55B2ED), 0xBC9333E3, -84},
  { /*  -216 */ UINT64_C (0x83C7088E1AAB65DB), 0x792667C7, -84},
  { /*  -215 */ UINT64_C (0x1A5B01B605557AC5), 0x7EA147F5, -83},
  { /*  -214 */ UINT64_C (0x34B6036C0AAAF58A), 0xFD428FE9, -83},
  { /*  -213 */ UINT64_C (0x696C06D81555EB15), 0xFA851FD2, -83},
  { /*  -212 */ UINT64_C (0xD2D80DB02AABD62B), 0xF50A3FA5, -83},
  { /*  -211 */ UINT64_C (0x2A2B35F00888C46F), 0x31020CBB, -82},
  { /*  -210 */ UINT64_C (0x54566BE0111188DE), 0x62041975, -82},
  { /*  -209 */ UINT64_C (0xA8ACD7C0222311BC), 0xC40832EA, -82},
  { /*  -208 */ UINT64_C (0x21BC2B266D3A36BF), 0x5A680A2F, -81},
  { /*  -207 */ UINT64_C (0x4378564CDA746D7E), 0xB4D0145E, -81},
  { /*  -206 */ UINT64_C (0x86F0AC99B4E8DAFD), 0x69A028BB, -81},
  { /*  -205 */ UINT64_C (0x1AFCEF51F0FB5EFF), 0x7B866E8C, -80},
  { /*  -204 */ UINT64_C (0x35F9DEA3E1F6BDFE), 0xF70CDD18, -80},
  { /*  -203 */ UINT64_C (0x6BF3BD47C3ED7BFD), 0xEE19BA2F, -80},
  { /*  -202 */ UINT64_C (0xD7E77A8F87DAF7FB), 0xDC33745F, -80},
  { /*  -201 */ UINT64_C (0x2B2E4BB64E5EFE65), 0x9270B0E0, -79},
  { /*  -200 */ UINT64_C (0x565C976C9CBDFCCB), 0x24E161C0, -79},
  { /*  -199 */ UINT64_C (0xACB92ED9397BF996), 0x49C2C37F, -79},
  { /*  -198 */ UINT64_C (0x228B6FC50B7F31EA), 0xDB8D5A4D, -78},
  { /*  -197 */ UINT64_C (0x4516DF8A16FE63D5), 0xB71AB499, -78},
  { /*  -196 */ UINT64_C (0x8A2DBF142DFCC7AB), 0x6E356932, -78},
  { /*  -195 */ UINT64_C (0x1BA2BFD0D5FF5B22), 0x493DE1D7, -77},
  { /*  -194 */ UINT64_C (0x37457FA1ABFEB644), 0x927BC3AE, -77},
  { /*  -193 */ UINT64_C (0x6E8AFF4357FD6C89), 0x24F7875C, -77},
  { /*  -192 */ UINT64_C (0xDD15FE86AFFAD912), 0x49EF0EB7, -77},
  { /*  -191 */ UINT64_C (0x2C37994E23322B6A), 0x0EC96958, -76},
  { /*  -190 */ UINT64_C (0x586F329C466456D4), 0x1D92D2B0, -76},
  { /*  -189 */ UINT64_C (0xB0DE65388CC8ADA8), 0x3B25A55F, -76},
  { /*  -188 */ UINT64_C (0x235FADD81C2822BB), 0x3F078779, -75},
  { /*  -187 */ UINT64_C (0x46BF5BB038504576), 0x7E0F0EF3, -75},
  { /*  -186 */ UINT64_C (0x8D7EB76070A08AEC), 0xFC1E1DE6, -75},
  { /*  -185 */ UINT64_C (0x1C4C8B1349B9B562), 0x98D2D2C8, -74},
  { /*  -184 */ UINT64_C (0x3899162693736AC5), 0x31A5A58F, -74},
  { /*  -183 */ UINT64_C (0x71322C4D26E6D58A), 0x634B4B1E, -74},
  { /*  -182 */ UINT64_C (0xE264589A4DCDAB14), 0xC696963C, -74},
  { /*  -181 */ UINT64_C (0x2D4744EBA9292237), 0x5AEAEAD9, -73},
  { /*  -180 */ UINT64_C (0x5A8E89D75252446E), 0xB5D5D5B2, -73},
  { /*  -179 */ UINT64_C (0xB51D13AEA4A488DD), 0x6BABAB64, -73},
  { /*  -178 */ UINT64_C (0x243903EFBA874E92), 0xAF22557A, -72},
  { /*  -177 */ UINT64_C (0x487207DF750E9D25), 0x5E44AAF5, -72},
  { /*  -176 */ UINT64_C (0x90E40FBEEA1D3A4A), 0xBC8955E9, -72},
  { /*  -175 */ UINT64_C (0x1CFA698C95390BA8), 0x8C1B7795, -71},
  { /*  -174 */ UINT64_C (0x39F4D3192A721751), 0x1836EF2A, -71},
  { /*  -173 */ UINT64_C (0x73E9A63254E42EA2), 0x306DDE54, -71},
  { /*  -172 */ UINT64_C (0xE7D34C64A9C85D44), 0x60DBBCA8, -71},
  { /*  -171 */ UINT64_C (0x2E5D75ADBB8E790D), 0xACF8BF55, -70},
  { /*  -170 */ UINT64_C (0x5CBAEB5B771CF21B), 0x59F17EAA, -70},
  { /*  -169 */ UINT64_C (0xB975D6B6EE39E436), 0xB3E2FD54, -70},
  { /*  -168 */ UINT64_C (0x25179157C93EC73E), 0x23FA32AA, -69},
  { /*  -167 */ UINT64_C (0x4A2F22AF927D8E7C), 0x47F46555, -69},
  { /*  -166 */ UINT64_C (0x945E455F24FB1CF8), 0x8FE8CAA9, -69},
  { /*  -165 */ UINT64_C (0x1DAC74463A989F64), 0xE994F555, -68},
  { /*  -164 */ UINT64_C (0x3B58E88C75313EC9), 0xD329EAAA, -68},
  { /*  -163 */ UINT64_C (0x76B1D118EA627D93), 0xA653D554, -68},
  { /*  -162 */ UINT64_C (0xED63A231D4C4FB27), 0x4CA7AAA8, -68},
  { /*  -161 */ UINT64_C (0x2F7A53A390F4323B), 0x0F54BBBB, -67},
  { /*  -160 */ UINT64_C (0x5EF4A74721E86476), 0x1EA97777, -67},
  { /*  -159 */ UINT64_C (0xBDE94E8E43D0C8EC), 0x3D52EEED, -67},
  { /*  -158 */ UINT64_C (0x25FB761C73F68E95), 0xA5DD62FC, -66},
  { /*  -157 */ UINT64_C (0x4BF6EC38E7ED1D2B), 0x4BBAC5F8, -66},
  { /*  -156 */ UINT64_C (0x97EDD871CFDA3A56), 0x97758BF1, -66},
  { /*  -155 */ UINT64_C (0x1E62C4E38FF87211), 0x517DE8CA, -65},
  { /*  -154 */ UINT64_C (0x3CC589C71FF0E422), 0xA2FBD194, -65},
  { /*  -153 */ UINT64_C (0x798B138E3FE1C845), 0x45F7A327, -65},
  { /*  -152 */ UINT64_C (0xF316271C7FC3908A), 0x8BEF464E, -65},
  { /*  -151 */ UINT64_C (0x309E07D27FF3E9B5), 0x4F2FDADC, -64},
  { /*  -150 */ UINT64_C (0x613C0FA4FFE7D36A), 0x9E5FB5B9, -64},
  { /*  -149 */ UINT64_C (0xC2781F49FFCFA6D5), 0x3CBF6B72, -64},
  { /*  -148 */ UINT64_C (0x26E4D30ECCC3215D), 0xD8F3157D, -63},
  { /*  -147 */ UINT64_C (0x4DC9A61D998642BB), 0xB1E62AFA, -63},
  { /*  -146 */ UINT64_C (0x9B934C3B330C8577), 0x63CC55F5, -63},
  { /*  -145 */ UINT64_C (0x1F1D75A5709C1AB1), 0x7A5C1131, -62},
  { /*  -144 */ UINT64_C (0x3E3AEB4AE1383562), 0xF4B82262, -62},
  { /*  -143 */ UINT64_C (0x7C75D695C2706AC5), 0xE97044C4, -62},
  { /*  -142 */ UINT64_C (0xF8EBAD2B84E0D58B), 0xD2E08987, -62},
  { /*  -141 */ UINT64_C (0x31C8BC3BE7602AB5), 0x90934EB5, -61},
  { /*  -140 */ UINT64_C (0x63917877CEC0556B), 0x21269D69, -61},
  { /*  -139 */ UINT64_C (0xC722F0EF9D80AAD6), 0x424D3AD3, -61},
  { /*  -138 */ UINT64_C (0x27D3C9C985E68891), 0x4075D891, -60},
  { /*  -137 */ UINT64_C (0x4FA793930BCD1122), 0x80EBB121, -60},
  { /*  -136 */ UINT64_C (0x9F4F2726179A2245), 0x01D76242, -60},
  { /*  -135 */ UINT64_C (0x1FDCA16E04B86D41), 0x005E46DA, -59},
  { /*  -134 */ UINT64_C (0x3FB942DC0970DA82), 0x00BC8DB4, -59},
  { /*  -133 */ UINT64_C (0x7F7285B812E1B504), 0x01791B68, -59},
  { /*  -132 */ UINT64_C (0xFEE50B7025C36A08), 0x02F236D0, -59},
  { /*  -131 */ UINT64_C (0x32FA9BE33AC0AECE), 0x66FD3E2A, -58},
  { /*  -130 */ UINT64_C (0x65F537C675815D9C), 0xCDFA7C53, -58},
  { /*  -129 */ UINT64_C (0xCBEA6F8CEB02BB39), 0x9BF4F8A7, -58},
  { /*  -128 */ UINT64_C (0x28C87CB5C89A2571), 0xEBFDCB55, -57},
  { /*  -127 */ UINT64_C (0x5190F96B91344AE3), 0xD7FB96A9, -57},
  { /*  -126 */ UINT64_C (0xA321F2D7226895C7), 0xAFF72D52, -57},
  { /*  -125 */ UINT64_C (0x20A063C4A07B5127), 0xEFFE3C44, -56},
  { /*  -124 */ UINT64_C (0x4140C78940F6A24F), 0xDFFC7887, -56},
  { /*  -123 */ UINT64_C (0x82818F1281ED449F), 0xBFF8F10E, -56},
  { /*  -122 */ UINT64_C (0x1A19E96A19FC40EC), 0xBFFE969C, -55},
  { /*  -121 */ UINT64_C (0x3433D2D433F881D9), 0x7FFD2D39, -55},
  { /*  -120 */ UINT64_C (0x6867A5A867F103B2), 0xFFFA5A72, -55},
  { /*  -119 */ UINT64_C (0xD0CF4B50CFE20765), 0xFFF4B4E4, -55},
  { /*  -118 */ UINT64_C (0x29C30F1029939B14), 0x6664242E, -54},
  { /*  -117 */ UINT64_C (0x53861E2053273628), 0xCCC8485B, -54},
  { /*  -116 */ UINT64_C (0xA70C3C40A64E6C51), 0x999090B6, -54},
  { /*  -115 */ UINT64_C (0x2168D8D9BADC7C10), 0x51E9B68B, -53},
  { /*  -114 */ UINT64_C (0x42D1B1B375B8F820), 0xA3D36D16, -53},
  { /*  -113 */ UINT64_C (0x85A36366EB71F041), 0x47A6DA2B, -53},
  { /*  -112 */ UINT64_C (0x1ABA4714957D300D), 0x0E549209, -52},
  { /*  -111 */ UINT64_C (0x35748E292AFA601A), 0x1CA92411, -52},
  { /*  -110 */ UINT64_C (0x6AE91C5255F4C034), 0x39524823, -52},
  { /*  -109 */ UINT64_C (0xD5D238A4ABE98068), 0x72A49046, -52},
  { /*  -108 */ UINT64_C (0x2AC3A4EDBBFB8014), 0xE3BA8341, -51},
  { /*  -107 */ UINT64_C (0x558749DB77F70029), 0xC7750682, -51},
  { /*  -106 */ UINT64_C (0xAB0E93B6EFEE0053), 0x8EEA0D04, -51},
  { /*  -105 */ UINT64_C (0x22361D8AFCC93343), 0xE962029A, -50},
  { /*  -104 */ UINT64_C (0x446C3B15F9926687), 0xD2C40535, -50},
  { /*  -103 */ UINT64_C (0x88D8762BF324CD0F), 0xA5880A6A, -50},
  { /*  -102 */ UINT64_C (0x1B5E7E08CA3A8F69), 0x87819BAF, -49},
  { /*  -101 */ UINT64_C (0x36BCFC1194751ED3), 0x0F03375E, -49},
  { /*  -100 */ UINT64_C (0x6D79F82328EA3DA6), 0x1E066EBB, -49},
  { /*   -99 */ UINT64_C (0xDAF3F04651D47B4C), 0x3C0CDD76, -49},
  { /*   -98 */ UINT64_C (0x2BCA63414390E575), 0xA59C2C4B, -48},
  { /*   -97 */ UINT64_C (0x5794C6828721CAEB), 0x4B385896, -48},
  { /*   -96 */ UINT64_C (0xAF298D050E4395D6), 0x9670B12B, -48},
  { /*   -95 */ UINT64_C (0x23084F676940B791), 0x5149BD09, -47},
  { /*   -94 */ UINT64_C (0x46109ECED2816F22), 0xA2937A11, -47},
  { /*   -93 */ UINT64_C (0x8C213D9DA502DE45), 0x4526F423, -47},
  { /*   -92 */ UINT64_C (0x1C06A5EC5433C60D), 0xDAA16407, -46},
  { /*   -91 */ UINT64_C (0x380D4BD8A8678C1B), 0xB542C80E, -46},
  { /*   -90 */ UINT64_C (0x701A97B150CF1837), 0x6A85901C, -46},
  { /*   -89 */ UINT64_C (0xE0352F62A19E306E), 0xD50B2038, -46},
  { /*   -88 */ UINT64_C (0x2CD76FE086B93CE2), 0xF768A00B, -45},
  { /*   -87 */ UINT64_C (0x59AEDFC10D7279C5), 0xEED14016, -45},
  { /*   -86 */ UINT64_C (0xB35DBF821AE4F38B), 0xDDA2802D, -45},
  { /*   -85 */ UINT64_C (0x23DF8CB39EFA971B), 0xF9208009, -44},
  { /*   -84 */ UINT64_C (0x47BF19673DF52E37), 0xF2410012, -44},
  { /*   -83 */ UINT64_C (0x8F7E32CE7BEA5C6F), 0xE4820024, -44},
  { /*   -82 */ UINT64_C (0x1CB2D6F618C878E3), 0x2DB399A1, -43},
  { /*   -81 */ UINT64_C (0x3965ADEC3190F1C6), 0x5B673341, -43},
  { /*   -80 */ UINT64_C (0x72CB5BD86321E38C), 0xB6CE6683, -43},
  { /*   -79 */ UINT64_C (0xE596B7B0C643C719), 0x6D9CCD06, -43},
  { /*   -78 */ UINT64_C (0x2DEAF189C140C16B), 0x7C528F68, -42},
  { /*   -77 */ UINT64_C (0x5BD5E313828182D6), 0xF8A51ECF, -42},
  { /*   -76 */ UINT64_C (0xB7ABC627050305AD), 0xF14A3D9E, -42},
  { /*   -75 */ UINT64_C (0x24BBF46E3433CDEF), 0x96A872B9, -41},
  { /*   -74 */ UINT64_C (0x4977E8DC68679BDF), 0x2D50E573, -41},
  { /*   -73 */ UINT64_C (0x92EFD1B8D0CF37BE), 0x5AA1CAE5, -41},
  { /*   -72 */ UINT64_C (0x1D6329F1C35CA4BF), 0xABB9F561, -40},
  { /*   -71 */ UINT64_C (0x3AC653E386B9497F), 0x5773EAC2, -40},
  { /*   -70 */ UINT64_C (0x758CA7C70D7292FE), 0xAEE7D584, -40},
  { /*   -69 */ UINT64_C (0xEB194F8E1AE525FD), 0x5DCFAB08, -40},
  { /*   -68 */ UINT64_C (0x2F050FE938943ACC), 0x45F65568, -39},
  { /*   -67 */ UINT64_C (0x5E0A1FD271287598), 0x8BECAAD0, -39},
  { /*   -66 */ UINT64_C (0xBC143FA4E250EB31), 0x17D955A0, -39},
  { /*   -65 */ UINT64_C (0x259DA6542D43623D), 0x04C51120, -38},
  { /*   -64 */ UINT64_C (0x4B3B4CA85A86C47A), 0x098A2240, -38},
  { /*   -63 */ UINT64_C (0x96769950B50D88F4), 0x13144480, -38},
  { /*   -62 */ UINT64_C (0x1E17B84357691B64), 0x03D0DA80, -37},
  { /*   -61 */ UINT64_C (0x3C2F7086AED236C8), 0x07A1B500, -37},
  { /*   -60 */ UINT64_C (0x785EE10D5DA46D90), 0x0F436A00, -37},
  { /*   -59 */ UINT64_C (0xF0BDC21ABB48DB20), 0x1E86D400, -37},
  { /*   -58 */ UINT64_C (0x3025F39EF241C56C), 0xD2E7C400, -36},
  { /*   -57 */ UINT64_C (0x604BE73DE4838AD9), 0xA5CF8800, -36},
  { /*   -56 */ UINT64_C (0xC097CE7BC90715B3), 0x4B9F1000, -36},
  { /*   -55 */ UINT64_C (0x2684C2E58E9B0457), 0x0F1FD000, -35},
  { /*   -54 */ UINT64_C (0x4D0985CB1D3608AE), 0x1E3FA000, -35},
  { /*   -53 */ UINT64_C (0x9A130B963A6C115C), 0x3C7F4000, -35},
  { /*   -52 */ UINT64_C (0x1ED09BEAD87C0378), 0xD8E64000, -34},
  { /*   -51 */ UINT64_C (0x3DA137D5B0F806F1), 0xB1CC8000, -34},
  { /*   -50 */ UINT64_C (0x7B426FAB61F00DE3), 0x63990000, -34},
  { /*   -49 */ UINT64_C (0xF684DF56C3E01BC6), 0xC7320000, -34},
  { /*   -48 */ UINT64_C (0x314DC6448D9338C1), 0x5B0A0000, -33},
  { /*   -47 */ UINT64_C (0x629B8C891B267182), 0xB6140000, -33},
  { /*   -46 */ UINT64_C (0xC5371912364CE305), 0x6C280000, -33},
  { /*   -45 */ UINT64_C (0x27716B6A0ADC2D67), 0x7C080000, -32},
  { /*   -44 */ UINT64_C (0x4EE2D6D415B85ACE), 0xF8100000, -32},
  { /*   -43 */ UINT64_C (0x9DC5ADA82B70B59D), 0xF0200000, -32},
  { /*   -42 */ UINT64_C (0x1F8DEF8808B02452), 0xC9A00000, -31},
  { /*   -41 */ UINT64_C (0x3F1BDF10116048A5), 0x93400000, -31},
  { /*   -40 */ UINT64_C (0x7E37BE2022C0914B), 0x26800000, -31},
  { /*   -39 */ UINT64_C (0xFC6F7C4045812296), 0x4D000000, -31},
  { /*   -38 */ UINT64_C (0x327CB2734119D3B7), 0xA9000000, -30},
  { /*   -37 */ UINT64_C (0x64F964E68233A76F), 0x52000000, -30},
  { /*   -36 */ UINT64_C (0xC9F2C9CD04674EDE), 0xA4000000, -30},
  { /*   -35 */ UINT64_C (0x2863C1F5CDAE42F9), 0x54000000, -29},
  { /*   -34 */ UINT64_C (0x50C783EB9B5C85F2), 0xA8000000, -29},
  { /*   -33 */ UINT64_C (0xA18F07D736B90BE5), 0x50000000, -29},
  { /*   -32 */ UINT64_C (0x204FCE5E3E250261), 0x10000000, -28},
  { /*   -31 */ UINT64_C (0x409F9CBC7C4A04C2), 0x20000000, -28},
  { /*   -30 */ UINT64_C (0x813F3978F8940984), 0x40000000, -28},
  { /*   -29 */ UINT64_C (0x19D971E4FE8401E7), 0x40000000, -27},
  { /*   -28 */ UINT64_C (0x33B2E3C9FD0803CE), 0x80000000, -27},
  { /*   -27 */ UINT64_C (0x6765C793FA10079D), 0x00000000, -27},
  { /*   -26 */ UINT64_C (0xCECB8F27F4200F3A), 0x00000000, -27},
  { /*   -25 */ UINT64_C (0x295BE96E64066972), 0x00000000, -26},
  { /*   -24 */ UINT64_C (0x52B7D2DCC80CD2E4), 0x00000000, -26},
  { /*   -23 */ UINT64_C (0xA56FA5B99019A5C8), 0x00000000, -26},
  { /*   -22 */ UINT64_C (0x2116545850052128), 0x00000000, -25},
  { /*   -21 */ UINT64_C (0x422CA8B0A00A4250), 0x00000000, -25},
  { /*   -20 */ UINT64_C (0x84595161401484A0), 0x00000000, -25},
  { /*   -19 */ UINT64_C (0x1A784379D99DB420), 0x00000000, -24},
  { /*   -18 */ UINT64_C (0x34F086F3B33B6840), 0x00000000, -24},
  { /*   -17 */ UINT64_C (0x69E10DE76676D080), 0x00000000, -24},
  { /*   -16 */ UINT64_C (0xD3C21BCECCEDA100), 0x00000000, -24},
  { /*   -15 */ UINT64_C (0x2A5A058FC295ED00), 0x00000000, -23},
  { /*   -14 */ UINT64_C (0x54B40B1F852BDA00), 0x00000000, -23},
  { /*   -13 */ UINT64_C (0xA968163F0A57B400), 0x00000000, -23},
  { /*   -12 */ UINT64_C (0x21E19E0C9BAB2400), 0x00000000, -22},
  { /*   -11 */ UINT64_C (0x43C33C1937564800), 0x00000000, -22},
  { /*   -10 */ UINT64_C (0x878678326EAC9000), 0x00000000, -22},
  { /*    -9 */ UINT64_C (0x1B1AE4D6E2EF5000), 0x00000000, -21},
  { /*    -8 */ UINT64_C (0x3635C9ADC5DEA000), 0x00000000, -21},
  { /*    -7 */ UINT64_C (0x6C6B935B8BBD4000), 0x00000000, -21},
  { /*    -6 */ UINT64_C (0xD8D726B7177A8000), 0x00000000, -21},
  { /*    -5 */ UINT64_C (0x2B5E3AF16B188000), 0x00000000, -20},
  { /*    -4 */ UINT64_C (0x56BC75E2D6310000), 0x00000000, -20},
  { /*    -3 */ UINT64_C (0xAD78EBC5AC620000), 0x00000000, -20},
  { /*    -2 */ UINT64_C (0x22B1C8C1227A0000), 0x00000000, -19},
  { /*    -1 */ UINT64_C (0x4563918244F40000), 0x00000000, -19},
  { /*     0 */ UINT64_C (0x8AC7230489E80000), 0x00000000, -19},
  { /*     1 */ UINT64_C (0x1BC16D674EC80000), 0x00000000, -18},
  { /*     2 */ UINT64_C (0x3782DACE9D900000), 0x00000000, -18},
  { /*     3 */ UINT64_C (0x6F05B59D3B200000), 0x00000000, -18},
  { /*     4 */ UINT64_C (0xDE0B6B3A76400000), 0x00000000, -18},
  { /*     5 */ UINT64_C (0x2C68AF0BB1400000), 0x00000000, -17},
  { /*     6 */ UINT64_C (0x58D15E1762800000), 0x00000000, -17},
  { /*     7 */ UINT64_C (0xB1A2BC2EC5000000), 0x00000000, -17},
  { /*     8 */ UINT64_C (0x2386F26FC1000000), 0x00000000, -16},
  { /*     9 */ UINT64_C (0x470DE4DF82000000), 0x00000000, -16},
  { /*    10 */ UINT64_C (0x8E1BC9BF04000000), 0x00000000, -16},
  { /*    11 */ UINT64_C (0x1C6BF52634000000), 0x00000000, -15},
  { /*    12 */ UINT64_C (0x38D7EA4C68000000), 0x00000000, -15},
  { /*    13 */ UINT64_C (0x71AFD498D0000000), 0x00000000, -15},
  { /*    14 */ UINT64_C (0xE35FA931A0000000), 0x00000000, -15},
  { /*    15 */ UINT64_C (0x2D79883D20000000), 0x00000000, -14},
  { /*    16 */ UINT64_C (0x5AF3107A40000000), 0x00000000, -14},
  { /*    17 */ UINT64_C (0xB5E620F480000000), 0x00000000, -14},
  { /*    18 */ UINT64_C (0x246139CA80000000), 0x00000000, -13},
  { /*    19 */ UINT64_C (0x48C2739500000000), 0x00000000, -13},
  { /*    20 */ UINT64_C (0x9184E72A00000000), 0x00000000, -13},
  { /*    21 */ UINT64_C (0x1D1A94A200000000), 0x00000000, -12},
  { /*    22 */ UINT64_C (0x3A35294400000000), 0x00000000, -12},
  { /*    23 */ UINT64_C (0x746A528800000000), 0x00000000, -12},
  { /*    24 */ UINT64_C (0xE8D4A51000000000), 0x00000000, -12},
  { /*    25 */ UINT64_C (0x2E90EDD000000000), 0x00000000, -11},
  { /*    26 */ UINT64_C (0x5D21DBA000000000), 0x00000000, -11},
  { /*    27 */ UINT64_C (0xBA43B74000000000), 0x00000000, -11},
  { /*    28 */ UINT64_C (0x2540BE4000000000), 0x00000000, -10},
  { /*    29 */ UINT64_C (0x4A817C8000000000), 0x00000000, -10},
  { /*    30 */ UINT64_C (0x9502F90000000000), 0x00000000, -10},
  { /*    31 */ UINT64_C (0x1DCD650000000000), 0x00000000, -9},
  { /*    32 */ UINT64_C (0x3B9ACA0000000000), 0x00000000, -9},
  { /*    33 */ UINT64_C (0x7735940000000000), 0x00000000, -9},
  { /*    34 */ UINT64_C (0xEE6B280000000000), 0x00000000, -9},
  { /*    35 */ UINT64_C (0x2FAF080000000000), 0x00000000, -8},
  { /*    36 */ UINT64_C (0x5F5E100000000000), 0x00000000, -8},
  { /*    37 */ UINT64_C (0xBEBC200000000000), 0x00000000, -8},
  { /*    38 */ UINT64_C (0x2625A00000000000), 0x00000000, -7},
  { /*    39 */ UINT64_C (0x4C4B400000000000), 0x00000000, -7},
  { /*    40 */ UINT64_C (0x9896800000000000), 0x00000000, -7},
  { /*    41 */ UINT64_C (0x1E84800000000000), 0x00000000, -6},
  { /*    42 */ UINT64_C (0x3D09000000000000), 0x00000000, -6},
  { /*    43 */ UINT64_C (0x7A12000000000000), 0x00000000, -6},
  { /*    44 */ UINT64_C (0xF424000000000000), 0x00000000, -6},
  { /*    45 */ UINT64_C (0x30D4000000000000), 0x00000000, -5},
  { /*    46 */ UINT64_C (0x61A8000000000000), 0x00000000, -5},
  { /*    47 */ UINT64_C (0xC350000000000000), 0x00000000, -5},
  { /*    48 */ UINT64_C (0x2710000000000000), 0x00000000, -4},
  { /*    49 */ UINT64_C (0x4E20000000000000), 0x00000000, -4},
  { /*    50 */ UINT64_C (0x9C40000000000000), 0x00000000, -4},
  { /*    51 */ UINT64_C (0x1F40000000000000), 0x00000000, -3},
  { /*    52 */ UINT64_C (0x3E80000000000000), 0x00000000, -3},
  { /*    53 */ UINT64_C (0x7D00000000000000), 0x00000000, -3},
  { /*    54 */ UINT64_C (0xFA00000000000000), 0x00000000, -3},
  { /*    55 */ UINT64_C (0x3200000000000000), 0x00000000, -2},
  { /*    56 */ UINT64_C (0x6400000000000000), 0x00000000, -2},
  { /*    57 */ UINT64_C (0xC800000000000000), 0x00000000, -2},
  { /*    58 */ UINT64_C (0x2800000000000000), 0x00000000, -1},
  { /*    59 */ UINT64_C (0x5000000000000000), 0x00000000, -1},
  { /*    60 */ UINT64_C (0xA000000000000000), 0x00000000, -1},
  { /*    61 */ UINT64_C (0x2000000000000000), 0x00000000, 0},
  { /*    62 */ UINT64_C (0x4000000000000000), 0x00000000, 0},
  { /*    63 */ UINT64_C (0x8000000000000000), 0x00000000, 0},
  { /*    64 */ UINT64_C (0x1999999999999999), 0x9999999A, 1},
  { /*    65 */ UINT64_C (0x3333333333333333), 0x33333333, 1},
  { /*    66 */ UINT64_C (0x6666666666666666), 0x66666666, 1},
  { /*    67 */ UINT64_C (0xCCCCCCCCCCCCCCCC), 0xCCCCCCCD, 1},
  { /*    68 */ UINT64_C (0x28F5C28F5C28F5C2), 0x8F5C28F6, 2},
  { /*    69 */ UINT64_C (0x51EB851EB851EB85), 0x1EB851EC, 2},
  { /*    70 */ UINT64_C (0xA3D70A3D70A3D70A), 0x3D70A3D7, 2},
  { /*    71 */ UINT64_C (0x20C49BA5E353F7CE), 0xD916872B, 3},
  { /*    72 */ UINT64_C (0x4189374BC6A7EF9D), 0xB22D0E56, 3},
  { /*    73 */ UINT64_C (0x83126E978D4FDF3B), 0x645A1CAC, 3},
  { /*    74 */ UINT64_C (0x1A36E2EB1C432CA5), 0x7A786C22, 4},
  { /*    75 */ UINT64_C (0x346DC5D63886594A), 0xF4F0D845, 4},
  { /*    76 */ UINT64_C (0x68DB8BAC710CB295), 0xE9E1B08A, 4},
  { /*    77 */ UINT64_C (0xD1B71758E219652B), 0xD3C36113, 4},
  { /*    78 */ UINT64_C (0x29F16B11C6D1E108), 0xC3F3E037, 5},
  { /*    79 */ UINT64_C (0x53E2D6238DA3C211), 0x87E7C06E, 5},
  { /*    80 */ UINT64_C (0xA7C5AC471B478423), 0x0FCF80DC, 5},
  { /*    81 */ UINT64_C (0x218DEF416BDB1A6D), 0x698FE692, 6},
  { /*    82 */ UINT64_C (0x431BDE82D7B634DA), 0xD31FCD25, 6},
  { /*    83 */ UINT64_C (0x8637BD05AF6C69B5), 0xA63F9A4A, 6},
  { /*    84 */ UINT64_C (0x1AD7F29ABCAF4857), 0x87A6520F, 7},
  { /*    85 */ UINT64_C (0x35AFE535795E90AF), 0x0F4CA41E, 7},
  { /*    86 */ UINT64_C (0x6B5FCA6AF2BD215E), 0x1E99483B, 7},
  { /*    87 */ UINT64_C (0xD6BF94D5E57A42BC), 0x3D329076, 7},
  { /*    88 */ UINT64_C (0x2AF31DC4611873BF), 0x3F70834B, 8},
  { /*    89 */ UINT64_C (0x55E63B88C230E77E), 0x7EE10696, 8},
  { /*    90 */ UINT64_C (0xABCC77118461CEFC), 0xFDC20D2B, 8},
  { /*    91 */ UINT64_C (0x225C17D04DAD2965), 0xCC5A02A2, 9},
  { /*    92 */ UINT64_C (0x44B82FA09B5A52CB), 0x98B40544, 9},
  { /*    93 */ UINT64_C (0x89705F4136B4A597), 0x31680A89, 9},
  { /*    94 */ UINT64_C (0x1B7CDFD9D7BDBAB7), 0xD6AE6882, 10},
  { /*    95 */ UINT64_C (0x36F9BFB3AF7B756F), 0xAD5CD104, 10},
  { /*    96 */ UINT64_C (0x6DF37F675EF6EADF), 0x5AB9A207, 10},
  { /*    97 */ UINT64_C (0xDBE6FECEBDEDD5BE), 0xB573440E, 10},
  { /*    98 */ UINT64_C (0x2BFAFFC2F2C92ABF), 0xBDE3DA69, 11},
  { /*    99 */ UINT64_C (0x57F5FF85E592557F), 0x7BC7B4D3, 11},
  { /*   100 */ UINT64_C (0xAFEBFF0BCB24AAFE), 0xF78F69A5, 11},
  { /*   101 */ UINT64_C (0x232F33025BD42232), 0xFE4FE1EE, 12},
  { /*   102 */ UINT64_C (0x465E6604B7A84465), 0xFC9FC3DC, 12},
  { /*   103 */ UINT64_C (0x8CBCCC096F5088CB), 0xF93F87B7, 12},
  { /*   104 */ UINT64_C (0x1C25C268497681C2), 0x650CB4BE, 13},
  { /*   105 */ UINT64_C (0x384B84D092ED0384), 0xCA19697D, 13},
  { /*   106 */ UINT64_C (0x709709A125DA0709), 0x9432D2F9, 13},
  { /*   107 */ UINT64_C (0xE12E13424BB40E13), 0x2865A5F2, 13},
  { /*   108 */ UINT64_C (0x2D09370D42573603), 0xD4E12130, 14},
  { /*   109 */ UINT64_C (0x5A126E1A84AE6C07), 0xA9C24261, 14},
  { /*   110 */ UINT64_C (0xB424DC35095CD80F), 0x538484C2, 14},
  { /*   111 */ UINT64_C (0x24075F3DCEAC2B36), 0x43E74DC0, 15},
  { /*   112 */ UINT64_C (0x480EBE7B9D58566C), 0x87CE9B81, 15},
  { /*   113 */ UINT64_C (0x901D7CF73AB0ACD9), 0x0F9D3701, 15},
  { /*   114 */ UINT64_C (0x1CD2B297D889BC2B), 0x6985D7CD, 16},
  { /*   115 */ UINT64_C (0x39A5652FB1137856), 0xD30BAF9A, 16},
  { /*   116 */ UINT64_C (0x734ACA5F6226F0AD), 0xA6175F34, 16},
  { /*   117 */ UINT64_C (0xE69594BEC44DE15B), 0x4C2EBE68, 16},
  { /*   118 */ UINT64_C (0x2E1DEA8C8DA92D12), 0x426FBFAE, 17},
  { /*   119 */ UINT64_C (0x5C3BD5191B525A24), 0x84DF7F5D, 17},
  { /*   120 */ UINT64_C (0xB877AA3236A4B449), 0x09BEFEBA, 17},
  { /*   121 */ UINT64_C (0x24E4BBA3A4875741), 0xCEBFCC8C, 18},
  { /*   122 */ UINT64_C (0x49C97747490EAE83), 0x9D7F9917, 18},
  { /*   123 */ UINT64_C (0x9392EE8E921D5D07), 0x3AFF322E, 18},
  { /*   124 */ UINT64_C (0x1D83C94FB6D2AC34), 0xA5663D3C, 19},
  { /*   125 */ UINT64_C (0x3B07929F6DA55869), 0x4ACC7A79, 19},
  { /*   126 */ UINT64_C (0x760F253EDB4AB0D2), 0x9598F4F2, 19},
  { /*   127 */ UINT64_C (0xEC1E4A7DB69561A5), 0x2B31E9E4, 19},
  { /*   128 */ UINT64_C (0x2F394219248446BA), 0xA23D2EC7, 20},
  { /*   129 */ UINT64_C (0x5E72843249088D75), 0x447A5D8E, 20},
  { /*   130 */ UINT64_C (0xBCE5086492111AEA), 0x88F4BB1D, 20},
  { /*   131 */ UINT64_C (0x25C768141D369EFB), 0xB4FDBF06, 21},
  { /*   132 */ UINT64_C (0x4B8ED0283A6D3DF7), 0x69FB7E0B, 21},
  { /*   133 */ UINT64_C (0x971DA05074DA7BEE), 0xD3F6FC17, 21},
  { /*   134 */ UINT64_C (0x1E392010175EE596), 0x2A6498D1, 22},
  { /*   135 */ UINT64_C (0x3C7240202EBDCB2C), 0x54C931A3, 22},
  { /*   136 */ UINT64_C (0x78E480405D7B9658), 0xA9926346, 22},
  { /*   137 */ UINT64_C (0xF1C90080BAF72CB1), 0x5324C68B, 22},
  { /*   138 */ UINT64_C (0x305B66802564A289), 0xDD6DC14F, 23},
  { /*   139 */ UINT64_C (0x60B6CD004AC94513), 0xBADB829E, 23},
  { /*   140 */ UINT64_C (0xC16D9A0095928A27), 0x75B7053C, 23},
  { /*   141 */ UINT64_C (0x26AF8533511D4ED4), 0xB1249AA6, 24},
  { /*   142 */ UINT64_C (0x4D5F0A66A23A9DA9), 0x6249354B, 24},
  { /*   143 */ UINT64_C (0x9ABE14CD44753B52), 0xC4926A96, 24},
  { /*   144 */ UINT64_C (0x1EF2D0F5DA7DD8AA), 0x27507BB8, 25},
  { /*   145 */ UINT64_C (0x3DE5A1EBB4FBB154), 0x4EA0F76F, 25},
  { /*   146 */ UINT64_C (0x7BCB43D769F762A8), 0x9D41EEDF, 25},
  { /*   147 */ UINT64_C (0xF79687AED3EEC551), 0x3A83DDBE, 25},
  { /*   148 */ UINT64_C (0x318481895D962776), 0xA54D92C0, 26},
  { /*   149 */ UINT64_C (0x63090312BB2C4EED), 0x4A9B257F, 26},
  { /*   150 */ UINT64_C (0xC612062576589DDA), 0x95364AFE, 26},
  { /*   151 */ UINT64_C (0x279D346DE4781F92), 0x1DD7A899, 27},
  { /*   152 */ UINT64_C (0x4F3A68DBC8F03F24), 0x3BAF5132, 27},
  { /*   153 */ UINT64_C (0x9E74D1B791E07E48), 0x775EA265, 27},
  { /*   154 */ UINT64_C (0x1FB0F6BE50601941), 0xB17953AE, 28},
  { /*   155 */ UINT64_C (0x3F61ED7CA0C03283), 0x62F2A75C, 28},
  { /*   156 */ UINT64_C (0x7EC3DAF941806506), 0xC5E54EB7, 28},
  { /*   157 */ UINT64_C (0xFD87B5F28300CA0D), 0x8BCA9D6E, 28},
  { /*   158 */ UINT64_C (0x32B4BDFD4D668ECF), 0x825BB916, 29},
  { /*   159 */ UINT64_C (0x65697BFA9ACD1D9F), 0x04B7722C, 29},
  { /*   160 */ UINT64_C (0xCAD2F7F5359A3B3E), 0x096EE458, 29},
  { /*   161 */ UINT64_C (0x289097FDD7853F0C), 0x684960DE, 30},
  { /*   162 */ UINT64_C (0x51212FFBAF0A7E18), 0xD092C1BD, 30},
  { /*   163 */ UINT64_C (0xA2425FF75E14FC31), 0xA125837A, 30},
  { /*   164 */ UINT64_C (0x2073ACCB12D0FF3D), 0x203AB3E5, 31},
  { /*   165 */ UINT64_C (0x40E7599625A1FE7A), 0x407567CA, 31},
  { /*   166 */ UINT64_C (0x81CEB32C4B43FCF4), 0x80EACF95, 31},
  { /*   167 */ UINT64_C (0x19F623D5A8A73297), 0x4CFBC31E, 32},
  { /*   168 */ UINT64_C (0x33EC47AB514E652E), 0x99F7863B, 32},
  { /*   169 */ UINT64_C (0x67D88F56A29CCA5D), 0x33EF0C77, 32},
  { /*   170 */ UINT64_C (0xCFB11EAD453994BA), 0x67DE18EE, 32},
  { /*   171 */ UINT64_C (0x2989D2EF743EB758), 0x7B2C6B63, 33},
  { /*   172 */ UINT64_C (0x5313A5DEE87D6EB0), 0xF658D6C5, 33},
  { /*   173 */ UINT64_C (0xA6274BBDD0FADD61), 0xECB1AD8B, 33},
  { /*   174 */ UINT64_C (0x213B0F25F69892AD), 0x2F56BC4F, 34},
  { /*   175 */ UINT64_C (0x42761E4BED31255A), 0x5EAD789E, 34},
  { /*   176 */ UINT64_C (0x84EC3C97DA624AB4), 0xBD5AF13C, 34},
  { /*   177 */ UINT64_C (0x1A95A5B7F87A0EF0), 0xF2ABC9D9, 35},
  { /*   178 */ UINT64_C (0x352B4B6FF0F41DE1), 0xE55793B2, 35},
  { /*   179 */ UINT64_C (0x6A5696DFE1E83BC3), 0xCAAF2763, 35},
  { /*   180 */ UINT64_C (0xD4AD2DBFC3D07787), 0x955E4EC6, 35},
  { /*   181 */ UINT64_C (0x2A8909265A5CE4B4), 0xB77942F4, 36},
  { /*   182 */ UINT64_C (0x5512124CB4B9C969), 0x6EF285E9, 36},
  { /*   183 */ UINT64_C (0xAA242499697392D2), 0xDDE50BD2, 36},
  { /*   184 */ UINT64_C (0x22073A8515171D5D), 0x5F943590, 37},
  { /*   185 */ UINT64_C (0x440E750A2A2E3ABA), 0xBF286B21, 37},
  { /*   186 */ UINT64_C (0x881CEA14545C7575), 0x7E50D641, 37},
  { /*   187 */ UINT64_C (0x1B38FB9DAA78E44A), 0xB2DCF7A7, 38},
  { /*   188 */ UINT64_C (0x3671F73B54F1C895), 0x65B9EF4D, 38},
  { /*   189 */ UINT64_C (0x6CE3EE76A9E3912A), 0xCB73DE9B, 38},
  { /*   190 */ UINT64_C (0xD9C7DCED53C72255), 0x96E7BD36, 38},
  { /*   191 */ UINT64_C (0x2B8E5F62AA5B06DD), 0xEAFB25D8, 39},
  { /*   192 */ UINT64_C (0x571CBEC554B60DBB), 0xD5F64BAF, 39},
  { /*   193 */ UINT64_C (0xAE397D8AA96C1B77), 0xABEC975E, 39},
  { /*   194 */ UINT64_C (0x22D84C4EEEAF38B1), 0x88C8EB13, 40},
  { /*   195 */ UINT64_C (0x45B0989DDD5E7163), 0x1191D626, 40},
  { /*   196 */ UINT64_C (0x8B61313BBABCE2C6), 0x2323AC4B, 40},
  { /*   197 */ UINT64_C (0x1BE03D0BF225C6F4), 0x6D6D88DC, 41},
  { /*   198 */ UINT64_C (0x37C07A17E44B8DE8), 0xDADB11B8, 41},
  { /*   199 */ UINT64_C (0x6F80F42FC8971BD1), 0xB5B6236F, 41},
  { /*   200 */ UINT64_C (0xDF01E85F912E37A3), 0x6B6C46DF, 41},
  { /*   201 */ UINT64_C (0x2C99FB46503C7187), 0x157C0E2D, 42},
  { /*   202 */ UINT64_C (0x5933F68CA078E30E), 0x2AF81C59, 42},
  { /*   203 */ UINT64_C (0xB267ED1940F1C61C), 0x55F038B2, 42},
  { /*   204 */ UINT64_C (0x23AE629EA696C138), 0xDDFCD824, 43},
  { /*   205 */ UINT64_C (0x475CC53D4D2D8271), 0xBBF9B047, 43},
  { /*   206 */ UINT64_C (0x8EB98A7A9A5B04E3), 0x77F3608F, 43},
  { /*   207 */ UINT64_C (0x1C8B821885456760), 0xB1971350, 44},
  { /*   208 */ UINT64_C (0x391704310A8ACEC1), 0x632E269F, 44},
  { /*   209 */ UINT64_C (0x722E086215159D82), 0xC65C4D3F, 44},
  { /*   210 */ UINT64_C (0xE45C10C42A2B3B05), 0x8CB89A7E, 44},
  { /*   211 */ UINT64_C (0x2DAC035A6ED57234), 0x4F581EE6, 45},
  { /*   212 */ UINT64_C (0x5B5806B4DDAAE468), 0x9EB03DCC, 45},
  { /*   213 */ UINT64_C (0xB6B00D69BB55C8D1), 0x3D607B98, 45},
  { /*   214 */ UINT64_C (0x24899C4858AAC1C3), 0x72ACE585, 46},
  { /*   215 */ UINT64_C (0x49133890B1558386), 0xE559CB0A, 46},
  { /*   216 */ UINT64_C (0x9226712162AB070D), 0xCAB39613, 46},
  { /*   217 */ UINT64_C (0x1D3AE36D13BBCE35), 0xF5571E04, 47},
  { /*   218 */ UINT64_C (0x3A75C6DA27779C6B), 0xEAAE3C08, 47},
  { /*   219 */ UINT64_C (0x74EB8DB44EEF38D7), 0xD55C780F, 47},
  { /*   220 */ UINT64_C (0xE9D71B689DDE71AF), 0xAAB8F01E, 47},
  { /*   221 */ UINT64_C (0x2EC49F14EC5FB056), 0x55583006, 48},
  { /*   222 */ UINT64_C (0x5D893E29D8BF60AC), 0xAAB0600C, 48},
  { /*   223 */ UINT64_C (0xBB127C53B17EC159), 0x5560C018, 48},
  { /*   224 */ UINT64_C (0x256A18DD89E626AB), 0x7779C005, 49},
  { /*   225 */ UINT64_C (0x4AD431BB13CC4D56), 0xEEF3800A, 49},
  { /*   226 */ UINT64_C (0x95A8637627989AAD), 0xDDE70013, 49},
  { /*   227 */ UINT64_C (0x1DEE7A4AD4B81EEF), 0x92C7CCD1, 50},
  { /*   228 */ UINT64_C (0x3BDCF495A9703DDF), 0x258F99A1, 50},
  { /*   229 */ UINT64_C (0x77B9E92B52E07BBE), 0x4B1F3343, 50},
  { /*   230 */ UINT64_C (0xEF73D256A5C0F77C), 0x963E6686, 50},
  { /*   231 */ UINT64_C (0x2FE3F6DE212697E5), 0xB7A61481, 51},
  { /*   232 */ UINT64_C (0x5FC7EDBC424D2FCB), 0x6F4C2902, 51},
  { /*   233 */ UINT64_C (0xBF8FDB78849A5F96), 0xDE985204, 51},
  { /*   234 */ UINT64_C (0x264FF8B1B41EDFEA), 0xF951AA01, 52},
  { /*   235 */ UINT64_C (0x4C9FF163683DBFD5), 0xF2A35402, 52},
  { /*   236 */ UINT64_C (0x993FE2C6D07B7FAB), 0xE546A804, 52},
  { /*   237 */ UINT64_C (0x1EA6608E29B24CBB), 0xFAA7BB34, 53},
  { /*   238 */ UINT64_C (0x3D4CC11C53649977), 0xF54F7668, 53},
  { /*   239 */ UINT64_C (0x7A998238A6C932EF), 0xEA9EECD0, 53},
  { /*   240 */ UINT64_C (0xF53304714D9265DF), 0xD53DD99F, 53},
  { /*   241 */ UINT64_C (0x310A3416A91D4793), 0x2AA5F853, 54},
  { /*   242 */ UINT64_C (0x6214682D523A8F26), 0x554BF0A6, 54},
  { /*   243 */ UINT64_C (0xC428D05AA4751E4C), 0xAA97E14C, 54},
  { /*   244 */ UINT64_C (0x273B5CDEEDB1060F), 0x55519376, 55},
  { /*   245 */ UINT64_C (0x4E76B9BDDB620C1E), 0xAAA326EB, 55},
  { /*   246 */ UINT64_C (0x9CED737BB6C4183D), 0x55464DD7, 55},
  { /*   247 */ UINT64_C (0x1F62B0B257C0D1A5), 0xDDDADC5E, 56},
  { /*   248 */ UINT64_C (0x3EC56164AF81A34B), 0xBBB5B8BC, 56},
  { /*   249 */ UINT64_C (0x7D8AC2C95F034697), 0x776B7178, 56},
  { /*   250 */ UINT64_C (0xFB158592BE068D2E), 0xEED6E2F1, 56},
  { /*   251 */ UINT64_C (0x3237811D593482A2), 0xFC916097, 57},
  { /*   252 */ UINT64_C (0x646F023AB2690545), 0xF922C12D, 57},
  { /*   253 */ UINT64_C (0xC8DE047564D20A8B), 0xF245825A, 57},
  { /*   254 */ UINT64_C (0x282C674AADC39BB5), 0x96DAB3AC, 58},
  { /*   255 */ UINT64_C (0x5058CE955B87376B), 0x2DB56757, 58},
  { /*   256 */ UINT64_C (0xA0B19D2AB70E6ED6), 0x5B6ACEAF, 58},
  { /*   257 */ UINT64_C (0x202385D557CFAFC4), 0x78AEF623, 59},
  { /*   258 */ UINT64_C (0x40470BAAAF9F5F88), 0xF15DEC46, 59},
  { /*   259 */ UINT64_C (0x808E17555F3EBF11), 0xE2BBD88C, 59},
  { /*   260 */ UINT64_C (0x19B604AAACA62636), 0xC6F25E82, 60},
  { /*   261 */ UINT64_C (0x336C0955594C4C6D), 0x8DE4BD05, 60},
  { /*   262 */ UINT64_C (0x66D812AAB29898DB), 0x1BC97A09, 60},
  { /*   263 */ UINT64_C (0xCDB02555653131B6), 0x3792F413, 60},
  { /*   264 */ UINT64_C (0x29233AAAADD6A38A), 0xD7EA30D1, 61},
  { /*   265 */ UINT64_C (0x524675555BAD4715), 0xAFD461A1, 61},
  { /*   266 */ UINT64_C (0xA48CEAAAB75A8E2B), 0x5FA8C342, 61},
  { /*   267 */ UINT64_C (0x20E8FBBBBE454FA2), 0x4654F3DA, 62},
  { /*   268 */ UINT64_C (0x41D1F7777C8A9F44), 0x8CA9E7B4, 62},
  { /*   269 */ UINT64_C (0x83A3EEEEF9153E89), 0x1953CF68, 62},
  { /*   270 */ UINT64_C (0x1A53FC9631D10C81), 0xD1DD8FE2, 63},
  { /*   271 */ UINT64_C (0x34A7F92C63A21903), 0xA3BB1FC3, 63},
  { /*   272 */ UINT64_C (0x694FF258C7443207), 0x47763F87, 63},
  { /*   273 */ UINT64_C (0xD29FE4B18E88640E), 0x8EEC7F0D, 63},
  { /*   274 */ UINT64_C (0x2A1FFA89E94E7A69), 0x4FC8E636, 64},
  { /*   275 */ UINT64_C (0x543FF513D29CF4D2), 0x9F91CC6C, 64},
  { /*   276 */ UINT64_C (0xA87FEA27A539E9A5), 0x3F2398D7, 64},
  { /*   277 */ UINT64_C (0x21B32ED4BAA52EBA), 0xA63A51C5, 65},
  { /*   278 */ UINT64_C (0x43665DA9754A5D75), 0x4C74A389, 65},
  { /*   279 */ UINT64_C (0x86CCBB52EA94BAEA), 0x98E94713, 65},
  { /*   280 */ UINT64_C (0x1AF5BF109550F22E), 0xEB61DB04, 66},
  { /*   281 */ UINT64_C (0x35EB7E212AA1E45D), 0xD6C3B607, 66},
  { /*   282 */ UINT64_C (0x6BD6FC425543C8BB), 0xAD876C0F, 66},
  { /*   283 */ UINT64_C (0xD7ADF884AA879177), 0x5B0ED81E, 66},
  { /*   284 */ UINT64_C (0x2B22CB4DBBB4B6B1), 0x789C91A0, 67},
  { /*   285 */ UINT64_C (0x5645969B77696D62), 0xF139233F, 67},
  { /*   286 */ UINT64_C (0xAC8B2D36EED2DAC5), 0xE272467E, 67},
  { /*   287 */ UINT64_C (0x22823C3E2FC3C55A), 0xC6E3A7B3, 68},
  { /*   288 */ UINT64_C (0x4504787C5F878AB5), 0x8DC74F66, 68},
  { /*   289 */ UINT64_C (0x8A08F0F8BF0F156B), 0x1B8E9ECB, 68},
  { /*   290 */ UINT64_C (0x1B9B6364F3030448), 0x9F1C8629, 69},
  { /*   291 */ UINT64_C (0x3736C6C9E6060891), 0x3E390C51, 69},
  { /*   292 */ UINT64_C (0x6E6D8D93CC0C1122), 0x7C7218A3, 69},
  { /*   293 */ UINT64_C (0xDCDB1B2798182244), 0xF8E43145, 69},
  { /*   294 */ UINT64_C (0x2C2BD23B1E6B3A0D), 0xCB60D6A7, 70},
  { /*   295 */ UINT64_C (0x5857A4763CD6741B), 0x96C1AD4F, 70},
  { /*   296 */ UINT64_C (0xB0AF48EC79ACE837), 0x2D835A9E, 70},
  { /*   297 */ UINT64_C (0x235641C8E52294D7), 0xD5E71220, 71},
  { /*   298 */ UINT64_C (0x46AC8391CA4529AF), 0xABCE243F, 71},
  { /*   299 */ UINT64_C (0x8D590723948A535F), 0x579C487E, 71},
  { /*   300 */ UINT64_C (0x1C45016D841BAA46), 0x44B8DB4C, 72},
  { /*   301 */ UINT64_C (0x388A02DB0837548C), 0x8971B699, 72},
  { /*   302 */ UINT64_C (0x711405B6106EA919), 0x12E36D32, 72},
  { /*   303 */ UINT64_C (0xE2280B6C20DD5232), 0x25C6DA64, 72},
  { /*   304 */ UINT64_C (0x2D3B357C0692AA0A), 0x078E2BAE, 73},
  { /*   305 */ UINT64_C (0x5A766AF80D255414), 0x0F1C575B, 73},
  { /*   306 */ UINT64_C (0xB4ECD5F01A4AA828), 0x1E38AEB6, 73},
  { /*   307 */ UINT64_C (0x242F5DFCD20EEE6E), 0x6C71BC8B, 74},
  { /*   308 */ UINT64_C (0x485EBBF9A41DDCDC), 0xD8E37916, 74},
  { /*   309 */ UINT64_C (0x90BD77F3483BB9B9), 0xB1C6F22B, 74},
  { /*   310 */ UINT64_C (0x1CF2B1970E725858), 0x56C163A2, 75},
  { /*   311 */ UINT64_C (0x39E5632E1CE4B0B0), 0xAD82C745, 75},
  { /*   312 */ UINT64_C (0x73CAC65C39C96161), 0x5B058E89, 75},
  { /*   313 */ UINT64_C (0xE7958CB87392C2C2), 0xB60B1D12, 75},
  { /*   314 */ UINT64_C (0x2E511C24E3EA26F3), 0xBE023904, 76},
  { /*   315 */ UINT64_C (0x5CA23849C7D44DE7), 0x7C047207, 76},
  { /*   316 */ UINT64_C (0xB94470938FA89BCE), 0xF808E40F, 76},
  { /*   317 */ UINT64_C (0x250DB01D8321B8C2), 0xFE682D9D, 77},
  { /*   318 */ UINT64_C (0x4A1B603B06437185), 0xFCD05B39, 77},
  { /*   319 */ UINT64_C (0x9436C0760C86E30B), 0xF9A0B672, 77},
  { /*   320 */ UINT64_C (0x1DA48CE468E7C702), 0x6520247D, 78},
  { /*   321 */ UINT64_C (0x3B4919C8D1CF8E04), 0xCA4048FA, 78},
  { /*   322 */ UINT64_C (0x76923391A39F1C09), 0x948091F5, 78},
  { /*   323 */ UINT64_C (0xED246723473E3813), 0x290123EA, 78},
  { /*   324 */ UINT64_C (0x2F6DAE3A4172D803), 0xD5003A62, 79},
  { /*   325 */ UINT64_C (0x5EDB5C7482E5B007), 0xAA0074C4, 79},
  { /*   326 */ UINT64_C (0xBDB6B8E905CB600F), 0x5400E988, 79},
  { /*   327 */ UINT64_C (0x25F1582E9AC24669), 0x773361E8, 80},
  { /*   328 */ UINT64_C (0x4BE2B05D35848CD2), 0xEE66C3D0, 80},
  { /*   329 */ UINT64_C (0x97C560BA6B0919A5), 0xDCCD87A0, 80},
  { /*   330 */ UINT64_C (0x1E5AACF215683854), 0x5F5C4E53, 81},
  { /*   331 */ UINT64_C (0x3CB559E42AD070A8), 0xBEB89CA6, 81},
  { /*   332 */ UINT64_C (0x796AB3C855A0E151), 0x7D71394D, 81},
  { /*   333 */ UINT64_C (0xF2D56790AB41C2A2), 0xFAE27299, 81},
  { /*   334 */ UINT64_C (0x309114B688A6C086), 0xFEFA16EB, 82},
  { /*   335 */ UINT64_C (0x6122296D114D810D), 0xFDF42DD7, 82},
  { /*   336 */ UINT64_C (0xC24452DA229B021B), 0xFBE85BAE, 82},
  { /*   337 */ UINT64_C (0x26DA76F86D52339F), 0x3261ABF0, 83},
  { /*   338 */ UINT64_C (0x4DB4EDF0DAA4673E), 0x64C357DF, 83},
  { /*   339 */ UINT64_C (0x9B69DBE1B548CE7C), 0xC986AFBE, 83},
  { /*   340 */ UINT64_C (0x1F152BF9F10E8FB2), 0x8EB4898C, 84},
  { /*   341 */ UINT64_C (0x3E2A57F3E21D1F65), 0x1D691319, 84},
  { /*   342 */ UINT64_C (0x7C54AFE7C43A3ECA), 0x3AD22632, 84},
  { /*   343 */ UINT64_C (0xF8A95FCF88747D94), 0x75A44C64, 84},
  { /*   344 */ UINT64_C (0x31BB798FE8174C50), 0xE4540F47, 85},
  { /*   345 */ UINT64_C (0x6376F31FD02E98A1), 0xC8A81E8E, 85},
  { /*   346 */ UINT64_C (0xC6EDE63FA05D3143), 0x91503D1C, 85},
  { /*   347 */ UINT64_C (0x27C92E0CB9AC3D0D), 0x8376729F, 86},
  { /*   348 */ UINT64_C (0x4F925C1973587A1B), 0x06ECE53F, 86},
  { /*   349 */ UINT64_C (0x9F24B832E6B0F436), 0x0DD9CA7D, 86},
  { /*   350 */ UINT64_C (0x1FD424D6FAF030D7), 0x9C5EC219, 87},
  { /*   351 */ UINT64_C (0x3FA849ADF5E061AF), 0x38BD8432, 87},
  { /*   352 */ UINT64_C (0x7F50935BEBC0C35E), 0x717B0864, 87},
  { /*   353 */ UINT64_C (0xFEA126B7D78186BC), 0xE2F610C8, 87},
  { /*   354 */ UINT64_C (0x32ED07BE5E4D1AF2), 0x93CAD028, 88},
  { /*   355 */ UINT64_C (0x65DA0F7CBC9A35E5), 0x2795A050, 88},
  { /*   356 */ UINT64_C (0xCBB41EF979346BCA), 0x4F2B40A0, 88},
  { /*   357 */ UINT64_C (0x28BD9FCB7EA4158E), 0xDCA24020, 89},
  { /*   358 */ UINT64_C (0x517B3F96FD482B1D), 0xB9448040, 89},
  { /*   359 */ UINT64_C (0xA2F67F2DFA90563B), 0x72890080, 89},
  { /*   360 */ UINT64_C (0x2097B309321CDE0B), 0xE3B5001A, 90},
  { /*   361 */ UINT64_C (0x412F66126439BC17), 0xC76A0033, 90},
  { /*   362 */ UINT64_C (0x825ECC24C873782F), 0x8ED40067, 90},
  { /*   363 */ UINT64_C (0x1A12F5A0F4E3E4D6), 0x4FC40015, 91},
  { /*   364 */ UINT64_C (0x3425EB41E9C7C9AC), 0x9F880029, 91},
  { /*   365 */ UINT64_C (0x684BD683D38F9359), 0x3F100052, 91},
  { /*   366 */ UINT64_C (0xD097AD07A71F26B2), 0x7E2000A4, 91},
  { /*   367 */ UINT64_C (0x29B7EF67EE396E23), 0xB2D33354, 92},
  { /*   368 */ UINT64_C (0x536FDECFDC72DC47), 0x65A666A8, 92},
  { /*   369 */ UINT64_C (0xA6DFBD9FB8E5B88E), 0xCB4CCD50, 92},
  { /*   370 */ UINT64_C (0x215FF2B98B6124E9), 0x5BDC2910, 93},
  { /*   371 */ UINT64_C (0x42BFE57316C249D2), 0xB7B85220, 93},
  { /*   372 */ UINT64_C (0x857FCAE62D8493A5), 0x6F70A440, 93},
  { /*   373 */ UINT64_C (0x1AB328946F80EA54), 0x497CEDA6, 94},
  { /*   374 */ UINT64_C (0x35665128DF01D4A8), 0x92F9DB4D, 94},
  { /*   375 */ UINT64_C (0x6ACCA251BE03A951), 0x25F3B69A, 94},
  { /*   376 */ UINT64_C (0xD59944A37C0752A2), 0x4BE76D33, 94},
  { /*   377 */ UINT64_C (0x2AB840ED7F34AA20), 0x7594AF71, 95},
  { /*   378 */ UINT64_C (0x557081DAFE695440), 0xEB295EE1, 95},
  { /*   379 */ UINT64_C (0xAAE103B5FCD2A881), 0xD652BDC3, 95},
  { /*   380 */ UINT64_C (0x222D00BDFF5D54E6), 0xC476F2C1, 96},
  { /*   381 */ UINT64_C (0x445A017BFEBAA9CD), 0x88EDE581, 96},
  { /*   382 */ UINT64_C (0x88B402F7FD75539B), 0x11DBCB02, 96},
  { /*   383 */ UINT64_C (0x1B5733CB32B110B8), 0x9D2BF567, 97},
  { /*   384 */ UINT64_C (0x36AE679665622171), 0x3A57EACE, 97},
  { /*   385 */ UINT64_C (0x6D5CCF2CCAC442E2), 0x74AFD59B, 97},
  { /*   386 */ UINT64_C (0xDAB99E59958885C4), 0xE95FAB37, 97},
  { /*   387 */ UINT64_C (0x2BBEB9451DE81AC0), 0xFB7988A5, 98},
  { /*   388 */ UINT64_C (0x577D728A3BD03581), 0xF6F31149, 98},
  { /*   389 */ UINT64_C (0xAEFAE51477A06B03), 0xEDE62292, 98},
  { /*   390 */ UINT64_C (0x22FEFA9DB1867BCD), 0x95FAD3B7, 99},
  { /*   391 */ UINT64_C (0x45FDF53B630CF79B), 0x2BF5A76E, 99},
  { /*   392 */ UINT64_C (0x8BFBEA76C619EF36), 0x57EB4EDB, 99},
  { /*   393 */ UINT64_C (0x1BFF2EE48E052FD7), 0xAB2F0FC5, 100},
  { /*   394 */ UINT64_C (0x37FE5DC91C0A5FAF), 0x565E1F8B, 100},
  { /*   395 */ UINT64_C (0x6FFCBB923814BF5E), 0xACBC3F16, 100},
  { /*   396 */ UINT64_C (0xDFF9772470297EBD), 0x59787E2C, 100},
  { /*   397 */ UINT64_C (0x2CCB7E3A7CD51959), 0x11E4E609, 101},
  { /*   398 */ UINT64_C (0x5996FC74F9AA32B2), 0x23C9CC11, 101},
  { /*   399 */ UINT64_C (0xB32DF8E9F3546564), 0x47939823, 101},
  { /*   400 */ UINT64_C (0x23D5FE9530AA7AAD), 0xA7EA51A1, 102},
  { /*   401 */ UINT64_C (0x47ABFD2A6154F55B), 0x4FD4A341, 102},
  { /*   402 */ UINT64_C (0x8F57FA54C2A9EAB6), 0x9FA94682, 102},
  { /*   403 */ UINT64_C (0x1CAB3210F3BB9557), 0xB988414D, 103},
  { /*   404 */ UINT64_C (0x39566421E7772AAF), 0x7310829B, 103},
  { /*   405 */ UINT64_C (0x72ACC843CEEE555E), 0xE6210535, 103},
  { /*   406 */ UINT64_C (0xE55990879DDCAABD), 0xCC420A6A, 103},
  { /*   407 */ UINT64_C (0x2DDEB68185F8EEF2), 0xC2739BAF, 104},
  { /*   408 */ UINT64_C (0x5BBD6D030BF1DDE5), 0x84E7375E, 104},
  { /*   409 */ UINT64_C (0xB77ADA0617E3BBCB), 0x09CE6EBB, 104},
  { /*   410 */ UINT64_C (0x24B22B9AD193F25B), 0xCEC2E2F2, 105},
  { /*   411 */ UINT64_C (0x49645735A327E4B7), 0x9D85C5E5, 105},
  { /*   412 */ UINT64_C (0x92C8AE6B464FC96F), 0x3B0B8BC9, 105},
  { /*   413 */ UINT64_C (0x1D5B561574765B7C), 0xA568B58F, 106},
  { /*   414 */ UINT64_C (0x3AB6AC2AE8ECB6F9), 0x4AD16B1D, 106},
  { /*   415 */ UINT64_C (0x756D5855D1D96DF2), 0x95A2D63A, 106},
  { /*   416 */ UINT64_C (0xEADAB0ABA3B2DBE5), 0x2B45AC75, 106},
  { /*   417 */ UINT64_C (0x2EF889BBED8A2BFA), 0xA241227E, 107},
  { /*   418 */ UINT64_C (0x5DF11377DB1457F5), 0x448244FC, 107},
  { /*   419 */ UINT64_C (0xBBE226EFB628AFEA), 0x890489F7, 107},
  { /*   420 */ UINT64_C (0x2593A163246E8995), 0x4E9A81FE, 108},
  { /*   421 */ UINT64_C (0x4B2742C648DD132A), 0x9D3503FC, 108},
  { /*   422 */ UINT64_C (0x964E858C91BA2655), 0x3A6A07F9, 108},
  { /*   423 */ UINT64_C (0x1E0FB44F50586E11), 0x0BAECE65, 109},
  { /*   424 */ UINT64_C (0x3C1F689EA0B0DC22), 0x175D9CCA, 109},
  { /*   425 */ UINT64_C (0x783ED13D4161B844), 0x2EBB3994, 109},
  { /*   426 */ UINT64_C (0xF07DA27A82C37088), 0x5D767328, 109},
  { /*   427 */ UINT64_C (0x3019207EE6F3E34E), 0x7917B0A2, 110},
  { /*   428 */ UINT64_C (0x603240FDCDE7C69C), 0xF22F6143, 110},
  { /*   429 */ UINT64_C (0xC06481FB9BCF8D39), 0xE45EC286, 110},
  { /*   430 */ UINT64_C (0x267A8065858FE90B), 0x9412F3B4, 111},
  { /*   431 */ UINT64_C (0x4CF500CB0B1FD217), 0x2825E769, 111},
  { /*   432 */ UINT64_C (0x99EA0196163FA42E), 0x504BCED2, 111},
  { /*   433 */ UINT64_C (0x1EC866B79E0CBA6F), 0xA9A8C2F7, 112},
  { /*   434 */ UINT64_C (0x3D90CD6F3C1974DF), 0x535185ED, 112},
  { /*   435 */ UINT64_C (0x7B219ADE7832E9BE), 0xA6A30BDB, 112},
  { /*   436 */ UINT64_C (0xF64335BCF065D37D), 0x4D4617B6, 112},
  { /*   437 */ UINT64_C (0x3140A458FCE12A4C), 0x42A79E58, 113},
  { /*   438 */ UINT64_C (0x628148B1F9C25498), 0x854F3CAF, 113},
  { /*   439 */ UINT64_C (0xC5029163F384A931), 0x0A9E795E, 113},
  { /*   440 */ UINT64_C (0x2766E9E0CA4DBB70), 0x3552E513, 114},
  { /*   441 */ UINT64_C (0x4ECDD3C1949B76E0), 0x6AA5CA26, 114},
  { /*   442 */ UINT64_C (0x9D9BA7832936EDC0), 0xD54B944C, 114},
  { /*   443 */ UINT64_C (0x1F8587E7083E2F8C), 0xF775840F, 115},
  { /*   444 */ UINT64_C (0x3F0B0FCE107C5F19), 0xEEEB081E, 115},
  { /*   445 */ UINT64_C (0x7E161F9C20F8BE33), 0xDDD6103C, 115},
  { /*   446 */ UINT64_C (0xFC2C3F3841F17C67), 0xBBAC2079, 115},
  { /*   447 */ UINT64_C (0x326F3FD80D304C14), 0xBF226CE5, 116},
  { /*   448 */ UINT64_C (0x64DE7FB01A609829), 0x7E44D9CA, 116},
  { /*   449 */ UINT64_C (0xC9BCFF6034C13052), 0xFC89B394, 116},
  { /*   450 */ UINT64_C (0x2858FFE00A8D09AA), 0x3281F0B7, 117},
  { /*   451 */ UINT64_C (0x50B1FFC0151A1354), 0x6503E16E, 117},
  { /*   452 */ UINT64_C (0xA163FF802A3426A8), 0xCA07C2DD, 117},
  { /*   453 */ UINT64_C (0x20473319A20A6E21), 0xC2018D5F, 118},
  { /*   454 */ UINT64_C (0x408E66334414DC43), 0x84031ABF, 118},
  { /*   455 */ UINT64_C (0x811CCC668829B887), 0x0806357D, 118},
  { /*   456 */ UINT64_C (0x19D28F47B4D524E7), 0xCE67A44C, 119},
  { /*   457 */ UINT64_C (0x33A51E8F69AA49CF), 0x9CCF4899, 119},
  { /*   458 */ UINT64_C (0x674A3D1ED354939F), 0x399E9131, 119},
  { /*   459 */ UINT64_C (0xCE947A3DA6A9273E), 0x733D2262, 119},
  { /*   460 */ UINT64_C (0x2950E53F87BB6E3F), 0xB0A5D3AD, 120},
  { /*   461 */ UINT64_C (0x52A1CA7F0F76DC7F), 0x614BA75A, 120},
  { /*   462 */ UINT64_C (0xA54394FE1EEDB8FE), 0xC2974EB5, 120},
  { /*   463 */ UINT64_C (0x210D8432D2FC5832), 0xF3B7DC8B, 121},
  { /*   464 */ UINT64_C (0x421B0865A5F8B065), 0xE76FB915, 121},
  { /*   465 */ UINT64_C (0x843610CB4BF160CB), 0xCEDF722A, 121},
  { /*   466 */ UINT64_C (0x1A71368F0F30468F), 0x295FE3A2, 122},
  { /*   467 */ UINT64_C (0x34E26D1E1E608D1E), 0x52BFC744, 122},
  { /*   468 */ UINT64_C (0x69C4DA3C3CC11A3C), 0xA57F8E88, 122},
  { /*   469 */ UINT64_C (0xD389B47879823479), 0x4AFF1D11, 122},
  { /*   470 */ UINT64_C (0x2A4EBDB1B1E6D74B), 0x75663903, 123},
  { /*   471 */ UINT64_C (0x549D7B6363CDAE96), 0xEACC7207, 123},
  { /*   472 */ UINT64_C (0xA93AF6C6C79B5D2D), 0xD598E40D, 123},
  { /*   473 */ UINT64_C (0x21D897C15B1F12A2), 0xC451C736, 124},
  { /*   474 */ UINT64_C (0x43B12F82B63E2545), 0x88A38E6C, 124},
  { /*   475 */ UINT64_C (0x87625F056C7C4A8B), 0x11471CD7, 124},
  { /*   476 */ UINT64_C (0x1B13AC9AAF4C0EE8), 0x9D0E38F8, 125},
  { /*   477 */ UINT64_C (0x362759355E981DD1), 0x3A1C71F0, 125},
  { /*   478 */ UINT64_C (0x6C4EB26ABD303BA2), 0x7438E3E0, 125},
  { /*   479 */ UINT64_C (0xD89D64D57A607744), 0xE871C7BF, 125},
  { /*   480 */ UINT64_C (0x2B52ADC44BACE4A7), 0x61B05B26, 126},
  { /*   481 */ UINT64_C (0x56A55B889759C94E), 0xC360B64C, 126},
  { /*   482 */ UINT64_C (0xAD4AB7112EB3929D), 0x86C16C99, 126},
  { /*   483 */ UINT64_C (0x22A88B036FBD83B9), 0x1AF37C1F, 127},
  { /*   484 */ UINT64_C (0x45511606DF7B0772), 0x35E6F83D, 127},
  { /*   485 */ UINT64_C (0x8AA22C0DBEF60EE4), 0x6BCDF07A, 127},
  { /*   486 */ UINT64_C (0x1BBA08CF8C979C94), 0x158F967F, 128},
  { /*   487 */ UINT64_C (0x3774119F192F3928), 0x2B1F2CFE, 128},
  { /*   488 */ UINT64_C (0x6EE8233E325E7250), 0x563E59FB, 128},
  { /*   489 */ UINT64_C (0xDDD0467C64BCE4A0), 0xAC7CB3F7, 128},
  { /*   490 */ UINT64_C (0x2C5CDAE5ADBF60EC), 0xEF4C23FE, 129},
  { /*   491 */ UINT64_C (0x58B9B5CB5B7EC1D9), 0xDE9847FC, 129},
  { /*   492 */ UINT64_C (0xB1736B96B6FD83B3), 0xBD308FF9, 129},
  { /*   493 */ UINT64_C (0x237D7BEAF165E723), 0xF2A34FFF, 130},
  { /*   494 */ UINT64_C (0x46FAF7D5E2CBCE47), 0xE5469FFD, 130},
  { /*   495 */ UINT64_C (0x8DF5EFABC5979C8F), 0xCA8D3FFA, 130},
  { /*   496 */ UINT64_C (0x1C6463225AB7EC1C), 0xC21C3FFF, 131},
  { /*   497 */ UINT64_C (0x38C8C644B56FD839), 0x84387FFE, 131},
  { /*   498 */ UINT64_C (0x71918C896ADFB073), 0x0870FFFB, 131},
  { /*   499 */ UINT64_C (0xE3231912D5BF60E6), 0x10E1FFF7, 131},
  { /*   500 */ UINT64_C (0x2D6D6B6A2ABFE02E), 0x03606665, 132},
  { /*   501 */ UINT64_C (0x5ADAD6D4557FC05C), 0x06C0CCC9, 132},
  { /*   502 */ UINT64_C (0xB5B5ADA8AAFF80B8), 0x0D819992, 132},
  { /*   503 */ UINT64_C (0x24578921BBCCB358), 0x02B3851D, 133},
  { /*   504 */ UINT64_C (0x48AF1243779966B0), 0x05670A3A, 133},
  { /*   505 */ UINT64_C (0x915E2486EF32CD60), 0x0ACE1475, 133},
  { /*   506 */ UINT64_C (0x1D12D41AFCA3C2AC), 0xCEF60417, 134},
  { /*   507 */ UINT64_C (0x3A25A835F9478559), 0x9DEC082F, 134},
  { /*   508 */ UINT64_C (0x744B506BF28F0AB3), 0x3BD8105D, 134},
  { /*   509 */ UINT64_C (0xE896A0D7E51E1566), 0x77B020BB, 134},
  { /*   510 */ UINT64_C (0x2E8486919439377A), 0xE4BCD359, 135},
  { /*   511 */ UINT64_C (0x5D090D2328726EF5), 0xC979A6B1, 135},
  { /*   512 */ UINT64_C (0xBA121A4650E4DDEB), 0x92F34D62, 135},
  { /*   513 */ UINT64_C (0x2536D20E102DC5FB), 0xEA30A914, 136},
  { /*   514 */ UINT64_C (0x4A6DA41C205B8BF7), 0xD4615227, 136},
  { /*   515 */ UINT64_C (0x94DB483840B717EF), 0xA8C2A44F, 136},
  { /*   516 */ UINT64_C (0x1DC574D80CF16B2F), 0xEE8D5410, 137},
  { /*   517 */ UINT64_C (0x3B8AE9B019E2D65F), 0xDD1AA81F, 137},
  { /*   518 */ UINT64_C (0x7715D36033C5ACBF), 0xBA35503F, 137},
  { /*   519 */ UINT64_C (0xEE2BA6C0678B597F), 0x746AA07E, 137},
  { /*   520 */ UINT64_C (0x2FA2548CE1824519), 0x7DAEECE6, 138},
  { /*   521 */ UINT64_C (0x5F44A919C3048A32), 0xFB5DD9CC, 138},
  { /*   522 */ UINT64_C (0xBE89523386091465), 0xF6BBB398, 138},
  { /*   523 */ UINT64_C (0x261B76D71ACE9DAD), 0xFE258A52, 139},
  { /*   524 */ UINT64_C (0x4C36EDAE359D3B5B), 0xFC4B14A3, 139},
  { /*   525 */ UINT64_C (0x986DDB5C6B3A76B7), 0xF8962946, 139},
  { /*   526 */ UINT64_C (0x1E7C5F127BD87E24), 0xCB513B74, 140},
  { /*   527 */ UINT64_C (0x3CF8BE24F7B0FC49), 0x96A276E9, 140},
  { /*   528 */ UINT64_C (0x79F17C49EF61F893), 0x2D44EDD2, 140},
  { /*   529 */ UINT64_C (0xF3E2F893DEC3F126), 0x5A89DBA4, 140},
  { /*   530 */ UINT64_C (0x30C6FE83F95A636E), 0x121B9254, 141},
  { /*   531 */ UINT64_C (0x618DFD07F2B4C6DC), 0x243724A8, 141},
  { /*   532 */ UINT64_C (0xC31BFA0FE5698DB8), 0x486E4950, 141},
  { /*   533 */ UINT64_C (0x2705986994484F8B), 0x41AFA843, 142},
  { /*   534 */ UINT64_C (0x4E0B30D328909F16), 0x835F5086, 142},
  { /*   535 */ UINT64_C (0x9C1661A651213E2D), 0x06BEA10D, 142},
  { /*   536 */ UINT64_C (0x1F37AD21436D0C6F), 0x67BFB9CF, 143},
  { /*   537 */ UINT64_C (0x3E6F5A4286DA18DE), 0xCF7F739F, 143},
  { /*   538 */ UINT64_C (0x7CDEB4850DB431BD), 0x9EFEE73D, 143},
  { /*   539 */ UINT64_C (0xF9BD690A1B68637B), 0x3DFDCE7B, 143},
  { /*   540 */ UINT64_C (0x31F2AE9B9F14E0B2), 0x3F99294C, 144},
  { /*   541 */ UINT64_C (0x63E55D373E29C164), 0x7F325297, 144},
  { /*   542 */ UINT64_C (0xC7CABA6E7C5382C8), 0xFE64A52F, 144},
  { /*   543 */ UINT64_C (0x27F5587C7F43E6F4), 0xFFADBAA3, 145},
  { /*   544 */ UINT64_C (0x4FEAB0F8FE87CDE9), 0xFF5B7546, 145},
  { /*   545 */ UINT64_C (0x9FD561F1FD0F9BD3), 0xFEB6EA8C, 145},
  { /*   546 */ UINT64_C (0x1FF779FD329CB8C3), 0xFFBE2EE9, 146},
  { /*   547 */ UINT64_C (0x3FEEF3FA65397187), 0xFF7C5DD2, 146},
  { /*   548 */ UINT64_C (0x7FDDE7F4CA72E30F), 0xFEF8BBA3, 146},
  { /*   549 */ UINT64_C (0xFFBBCFE994E5C61F), 0xFDF17746, 146},
  { /*   550 */ UINT64_C (0x33258FFB842DF46C), 0xCC637E41, 147},
  { /*   551 */ UINT64_C (0x664B1FF7085BE8D9), 0x98C6FC83, 147},
  { /*   552 */ UINT64_C (0xCC963FEE10B7D1B3), 0x318DF905, 147},
  { /*   553 */ UINT64_C (0x28EAD9960357F6BD), 0x704F9834, 148},
  { /*   554 */ UINT64_C (0x51D5B32C06AFED7A), 0xE09F3068, 148},
  { /*   555 */ UINT64_C (0xA3AB66580D5FDAF5), 0xC13E60D1, 148},
  { /*   556 */ UINT64_C (0x20BBE144CF799231), 0x26A6135D, 149},
  { /*   557 */ UINT64_C (0x4177C2899EF32462), 0x4D4C26BA, 149},
  { /*   558 */ UINT64_C (0x82EF85133DE648C4), 0x9A984D74, 149},
  { /*   559 */ UINT64_C (0x1A2FE76A3F9474F4), 0x1EEB42B1, 150},
  { /*   560 */ UINT64_C (0x345FCED47F28E9E8), 0x3DD68562, 150},
  { /*   561 */ UINT64_C (0x68BF9DA8FE51D3D0), 0x7BAD0AC3, 150},
  { /*   562 */ UINT64_C (0xD17F3B51FCA3A7A0), 0xF75A1586, 150},
  { /*   563 */ UINT64_C (0x29E63F1065BA54B9), 0xCB12044E, 151},
  { /*   564 */ UINT64_C (0x53CC7E20CB74A973), 0x9624089C, 151},
  { /*   565 */ UINT64_C (0xA798FC4196E952E7), 0x2C481138, 151},
  { /*   566 */ UINT64_C (0x2184FF405161DD61), 0x6F419D0B, 152},
  { /*   567 */ UINT64_C (0x4309FE80A2C3BAC2), 0xDE833A16, 152},
  { /*   568 */ UINT64_C (0x8613FD0145877585), 0xBD06742D, 152},
  { /*   569 */ UINT64_C (0x1AD0CC33744E4AB4), 0x59014A6F, 153},
  { /*   570 */ UINT64_C (0x35A19866E89C9568), 0xB20294DF, 153},
  { /*   571 */ UINT64_C (0x6B4330CDD1392AD1), 0x640529BE, 153},
  { /*   572 */ UINT64_C (0xD686619BA27255A2), 0xC80A537B, 153},
  { /*   573 */ UINT64_C (0x2AE7AD1F207D4453), 0xC19BAA4C, 154},
  { /*   574 */ UINT64_C (0x55CF5A3E40FA88A7), 0x83375498, 154},
  { /*   575 */ UINT64_C (0xAB9EB47C81F5114F), 0x066EA92F, 154},
  { /*   576 */ UINT64_C (0x2252F0E5B39769DC), 0x9AE2EEA3, 155},
  { /*   577 */ UINT64_C (0x44A5E1CB672ED3B9), 0x35C5DD46, 155},
  { /*   578 */ UINT64_C (0x894BC396CE5DA772), 0x6B8BBA8C, 155},
  { /*   579 */ UINT64_C (0x1B758D848FAC54B0), 0x7BE8BEE9, 156},
  { /*   580 */ UINT64_C (0x36EB1B091F58A960), 0xF7D17DD2, 156},
  { /*   581 */ UINT64_C (0x6DD636123EB152C1), 0xEFA2FBA3, 156},
  { /*   582 */ UINT64_C (0xDBAC6C247D62A583), 0xDF45F747, 156},
  { /*   583 */ UINT64_C (0x2BEF48D41913BAB3), 0xF97464A8, 157},
  { /*   584 */ UINT64_C (0x57DE91A832277567), 0xF2E8C94F, 157},
  { /*   585 */ UINT64_C (0xAFBD2350644EEACF), 0xE5D1929F, 157},
  { /*   586 */ UINT64_C (0x2325D3DCE0DC955C), 0xC7905086, 158},
  { /*   587 */ UINT64_C (0x464BA7B9C1B92AB9), 0x8F20A10C, 158},
  { /*   588 */ UINT64_C (0x8C974F7383725573), 0x1E414219, 158},
  { /*   589 */ UINT64_C (0x1C1E43171A4A1117), 0x060D0D38, 159},
  { /*   590 */ UINT64_C (0x383C862E3494222E), 0x0C1A1A70, 159},
  { /*   591 */ UINT64_C (0x70790C5C6928445C), 0x183434E1, 159},
  { /*   592 */ UINT64_C (0xE0F218B8D25088B8), 0x306869C1, 159},
  { /*   593 */ UINT64_C (0x2CFD3824F6DCE824), 0xD67B485A, 160},
  { /*   594 */ UINT64_C (0x59FA7049EDB9D049), 0xACF690B4, 160},
  { /*   595 */ UINT64_C (0xB3F4E093DB73A093), 0x59ED2167, 160},
  { /*   596 */ UINT64_C (0x23FDC683F8B0B9B7), 0x11FC39E1, 161},
  { /*   597 */ UINT64_C (0x47FB8D07F161736E), 0x23F873C3, 161},
  { /*   598 */ UINT64_C (0x8FF71A0FE2C2E6DC), 0x47F0E786, 161},
  { /*   599 */ UINT64_C (0x1CCB0536608D615F), 0x419694B4, 162},
  { /*   600 */ UINT64_C (0x39960A6CC11AC2BE), 0x832D2969, 162},
  { /*   601 */ UINT64_C (0x732C14D98235857D), 0x065A52D2, 162},
  { /*   602 */ UINT64_C (0xE65829B3046B0AFA), 0x0CB4A5A3, 162},
  { /*   603 */ UINT64_C (0x2E11A1F09A7BCEFE), 0xCF575454, 163},
  { /*   604 */ UINT64_C (0x5C2343E134F79DFD), 0x9EAEA8A8, 163},
  { /*   605 */ UINT64_C (0xB84687C269EF3BFB), 0x3D5D514F, 163},
  { /*   606 */ UINT64_C (0x24DAE7F3AEC97265), 0x72AC4376, 164},
  { /*   607 */ UINT64_C (0x49B5CFE75D92E4CA), 0xE55886ED, 164},
  { /*   608 */ UINT64_C (0x936B9FCEBB25C995), 0xCAB10DD9, 164},
  { /*   609 */ UINT64_C (0x1D7BECC2F23AC1EA), 0xC223692B, 165},
  { /*   610 */ UINT64_C (0x3AF7D985E47583D5), 0x8446D257, 165},
  { /*   611 */ UINT64_C (0x75EFB30BC8EB07AB), 0x088DA4AE, 165},
  { /*   612 */ UINT64_C (0xEBDF661791D60F56), 0x111B495B, 165},
  { /*   613 */ UINT64_C (0x2F2CAE04B6C46977), 0x9D057512, 166},
  { /*   614 */ UINT64_C (0x5E595C096D88D2EF), 0x3A0AEA24, 166},
  { /*   615 */ UINT64_C (0xBCB2B812DB11A5DE), 0x7415D449, 166},
  { /*   616 */ UINT64_C (0x25BD5803C569EDF9), 0x4A6AC40F, 167},
  { /*   617 */ UINT64_C (0x4B7AB0078AD3DBF2), 0x94D5881D, 167},
  { /*   618 */ UINT64_C (0x96F5600F15A7B7E5), 0x29AB103A, 167},
  { /*   619 */ UINT64_C (0x1E3113363787F194), 0x3B889CD8, 168},
  { /*   620 */ UINT64_C (0x3C62266C6F0FE328), 0x771139B1, 168},
  { /*   621 */ UINT64_C (0x78C44CD8DE1FC650), 0xEE227362, 168},
  { /*   622 */ UINT64_C (0xF18899B1BC3F8CA1), 0xDC44E6C4, 168},
  { /*   623 */ UINT64_C (0x304E85238C0CB5B9), 0xF8DA948E, 169},
  { /*   624 */ UINT64_C (0x609D0A4718196B73), 0xF1B5291B, 169},
  { /*   625 */ UINT64_C (0xC13A148E3032D6E7), 0xE36A5236, 169},
  { /*   626 */ UINT64_C (0x26A5374FA33D5E2E), 0x60AEDD3E, 170},
  { /*   627 */ UINT64_C (0x4D4A6E9F467ABC5C), 0xC15DBA7C, 170},
  { /*   628 */ UINT64_C (0x9A94DD3E8CF578B9), 0x82BB74F8, 170},
  { /*   629 */ UINT64_C (0x1EEA92A61C311825), 0x1A257DCB, 171},
  { /*   630 */ UINT64_C (0x3DD5254C3862304A), 0x344AFB96, 171},
  { /*   631 */ UINT64_C (0x7BAA4A9870C46094), 0x6895F72D, 171},
  { /*   632 */ UINT64_C (0xF7549530E188C128), 0xD12BEE5A, 171},
  { /*   633 */ UINT64_C (0x31775109C6B4F36E), 0x903BFC78, 172},
  { /*   634 */ UINT64_C (0x62EEA2138D69E6DD), 0x2077F8F1, 172},
  { /*   635 */ UINT64_C (0xC5DD44271AD3CDBA), 0x40EFF1E2, 172},
  { /*   636 */ UINT64_C (0x2792A73B055D8F8B), 0xA6966394, 173},
  { /*   637 */ UINT64_C (0x4F254E760ABB1F17), 0x4D2CC727, 173},
  { /*   638 */ UINT64_C (0x9E4A9CEC15763E2E), 0x9A598E4E, 173},
  { /*   639 */ UINT64_C (0x1FA885C8D117A609), 0x5211E943, 174},
  { /*   640 */ UINT64_C (0x3F510B91A22F4C12), 0xA423D286, 174},
  { /*   641 */ UINT64_C (0x7EA21723445E9825), 0x4847A50B, 174},
  { /*   642 */ UINT64_C (0xFD442E4688BD304A), 0x908F4A16, 174},
  { /*   643 */ UINT64_C (0x32A73C7481BF700E), 0xE9B64204, 175},
  { /*   644 */ UINT64_C (0x654E78E9037EE01D), 0xD36C8409, 175},
  { /*   645 */ UINT64_C (0xCA9CF1D206FDC03B), 0xA6D90812, 175},
  { /*   646 */ UINT64_C (0x2885C9F6CE32C00B), 0xEE2B6804, 176},
  { /*   647 */ UINT64_C (0x510B93ED9C658017), 0xDC56D007, 176},
  { /*   648 */ UINT64_C (0xA21727DB38CB002F), 0xB8ADA00E, 176},
  { /*   649 */ UINT64_C (0x206B07F8A4F5666F), 0xF1BC5336, 177},
  { /*   650 */ UINT64_C (0x40D60FF149EACCDF), 0xE378A66C, 177},
  { /*   651 */ UINT64_C (0x81AC1FE293D599BF), 0xC6F14CD8, 177},
  { /*   652 */ UINT64_C (0x19EF3993B72AB859), 0x8E304292, 178},
  { /*   653 */ UINT64_C (0x33DE73276E5570B3), 0x1C608523, 178},
  { /*   654 */ UINT64_C (0x67BCE64EDCAAE166), 0x38C10A47, 178},
  { /*   655 */ UINT64_C (0xCF79CC9DB955C2CC), 0x7182148D, 178},
  { /*   656 */ UINT64_C (0x297EC285F1DDF3C2), 0x7D1A041C, 179},
  { /*   657 */ UINT64_C (0x52FD850BE3BBE784), 0xFA340839, 179},
  { /*   658 */ UINT64_C (0xA5FB0A17C777CF09), 0xF4681071, 179},
  { /*   659 */ UINT64_C (0x21323537F4B18FCE), 0xCA7B367D, 180},
  { /*   660 */ UINT64_C (0x42646A6FE9631F9D), 0x94F66CFA, 180},
  { /*   661 */ UINT64_C (0x84C8D4DFD2C63F3B), 0x29ECD9F4, 180},
  { /*   662 */ UINT64_C (0x1A8E90F9908E0CA5), 0x6EC8F864, 181},
  { /*   663 */ UINT64_C (0x351D21F3211C194A), 0xDD91F0C8, 181},
  { /*   664 */ UINT64_C (0x6A3A43E642383295), 0xBB23E190, 181},
  { /*   665 */ UINT64_C (0xD47487CC8470652B), 0x7647C320, 181},
  { /*   666 */ UINT64_C (0x2A7DB4C280E3476F), 0x17A7F3D3, 182},
  { /*   667 */ UINT64_C (0x54FB698501C68EDE), 0x2F4FE7A6, 182},
  { /*   668 */ UINT64_C (0xA9F6D30A038D1DBC), 0x5E9FCF4D, 182},
  { /*   669 */ UINT64_C (0x21FE2A3533E905F2), 0x79532976, 183},
  { /*   670 */ UINT64_C (0x43FC546A67D20BE4), 0xF2A652EC, 183},
  { /*   671 */ UINT64_C (0x87F8A8D4CFA417C9), 0xE54CA5D7, 183},
  { /*   672 */ UINT64_C (0x1B31BB5DC320D18E), 0xC775BAC5, 184},
  { /*   673 */ UINT64_C (0x366376BB8641A31D), 0x8EEB7589, 184},
  { /*   674 */ UINT64_C (0x6CC6ED770C83463B), 0x1DD6EB12, 184},
  { /*   675 */ UINT64_C (0xD98DDAEE19068C76), 0x3BADD625, 184},
  { /*   676 */ UINT64_C (0x2B82C562D1CE1C17), 0xA5892AD4, 185},
  { /*   677 */ UINT64_C (0x57058AC5A39C382F), 0x4B1255A8, 185},
  { /*   678 */ UINT64_C (0xAE0B158B4738705E), 0x9624AB51, 185},
  { /*   679 */ UINT64_C (0x22CF044F0E3E7CDF), 0xB7A0EF10, 186},
  { /*   680 */ UINT64_C (0x459E089E1C7CF9BF), 0x6F41DE20, 186},
  { /*   681 */ UINT64_C (0x8B3C113C38F9F37E), 0xDE83BC41, 186},
  { /*   682 */ UINT64_C (0x1BD8D03F3E9863E6), 0x2C80BF40, 187},
  { /*   683 */ UINT64_C (0x37B1A07E7D30C7CC), 0x59017E80, 187},
  { /*   684 */ UINT64_C (0x6F6340FCFA618F98), 0xB202FD00, 187},
  { /*   685 */ UINT64_C (0xDEC681F9F4C31F31), 0x6405FA01, 187},
  { /*   686 */ UINT64_C (0x2C8E19FECA8D6CA3), 0x7A679867, 188},
  { /*   687 */ UINT64_C (0x591C33FD951AD946), 0xF4CF30CD, 188},
  { /*   688 */ UINT64_C (0xB23867FB2A35B28D), 0xE99E619A, 188},
  { /*   689 */ UINT64_C (0x23A4E198A20ABD4F), 0x951FAD1F, 189},
  { /*   690 */ UINT64_C (0x4749C33144157A9F), 0x2A3F5A3E, 189},
  { /*   691 */ UINT64_C (0x8E938662882AF53E), 0x547EB47B, 189},
  { /*   692 */ UINT64_C (0x1C83E7AD4E6EFDD9), 0x4419574C, 190},
  { /*   693 */ UINT64_C (0x3907CF5A9CDDFBB2), 0x8832AE98, 190},
  { /*   694 */ UINT64_C (0x720F9EB539BBF765), 0x10655D30, 190},
  { /*   695 */ UINT64_C (0xE41F3D6A7377EECA), 0x20CABA5F, 190},
  { /*   696 */ UINT64_C (0x2D9FD9154A4B2FC2), 0x068EF213, 191},
  { /*   697 */ UINT64_C (0x5B3FB22A94965F84), 0x0D1DE426, 191},
  { /*   698 */ UINT64_C (0xB67F6455292CBF08), 0x1A3BC84C, 191},
  { /*   699 */ UINT64_C (0x247FE0DDD508F301), 0x9ED8C1A9, 192},
  { /*   700 */ UINT64_C (0x48FFC1BBAA11E603), 0x3DB18352, 192},
  { /*   701 */ UINT64_C (0x91FF83775423CC06), 0x7B6306A3, 192},
  { /*   702 */ UINT64_C (0x1D331A4B10D3F59A), 0xE57A3487, 193},
  { /*   703 */ UINT64_C (0x3A66349621A7EB35), 0xCAF4690E, 193},
  { /*   704 */ UINT64_C (0x74CC692C434FD66B), 0x95E8D21C, 193},
  { /*   705 */ UINT64_C (0xE998D258869FACD7), 0x2BD1A438, 193},
  { /*   706 */ UINT64_C (0x2EB82A11B48655C4), 0xA25D20D8, 194},
  { /*   707 */ UINT64_C (0x5D705423690CAB89), 0x44BA41B0, 194},
  { /*   708 */ UINT64_C (0xBAE0A846D2195712), 0x89748360, 194},
  { /*   709 */ UINT64_C (0x256021A7C39EAB03), 0xB5174D7A, 195},
  { /*   710 */ UINT64_C (0x4AC0434F873D5607), 0x6A2E9AF3, 195},
  { /*   711 */ UINT64_C (0x9580869F0E7AAC0E), 0xD45D35E7, 195},
  { /*   712 */ UINT64_C (0x1DE6815302E5559C), 0x90DF712E, 196},
  { /*   713 */ UINT64_C (0x3BCD02A605CAAB39), 0x21BEE25C, 196},
  { /*   714 */ UINT64_C (0x779A054C0B955672), 0x437DC4B9, 196},
  { /*   715 */ UINT64_C (0xEF340A98172AACE4), 0x86FB8971, 196},
  { /*   716 */ UINT64_C (0x2FD735519E3BBC2D), 0xB498B517, 197},
  { /*   717 */ UINT64_C (0x5FAE6AA33C77785B), 0x69316A2D, 197},
  { /*   718 */ UINT64_C (0xBF5CD54678EEF0B6), 0xD262D45A, 197},
  { /*   719 */ UINT64_C (0x2645C4414B62FCF1), 0x5D46F745, 198},
  { /*   720 */ UINT64_C (0x4C8B888296C5F9E2), 0xBA8DEE8B, 198},
  { /*   721 */ UINT64_C (0x991711052D8BF3C5), 0x751BDD15, 198},
  { /*   722 */ UINT64_C (0x1E9E369AA2B59727), 0x7DD25F6B, 199},
  { /*   723 */ UINT64_C (0x3D3C6D35456B2E4E), 0xFBA4BED5, 199},
  { /*   724 */ UINT64_C (0x7A78DA6A8AD65C9D), 0xF7497DAB, 199},
  { /*   725 */ UINT64_C (0xF4F1B4D515ACB93B), 0xEE92FB55, 199},
  { /*   726 */ UINT64_C (0x30FD242A9DEF583F), 0x2FB6FF11, 200},
  { /*   727 */ UINT64_C (0x61FA48553BDEB07E), 0x5F6DFE22, 200},
  { /*   728 */ UINT64_C (0xC3F490AA77BD60FC), 0xBEDBFC44, 200},
  { /*   729 */ UINT64_C (0x2730E9BBB18C4698), 0xF2F8CC0E, 201},
  { /*   730 */ UINT64_C (0x4E61D37763188D31), 0xE5F1981B, 201},
  { /*   731 */ UINT64_C (0x9CC3A6EEC6311A63), 0xCBE33036, 201},
  { /*   732 */ UINT64_C (0x1F5A549627A36BAD), 0x8F2D700B, 202},
  { /*   733 */ UINT64_C (0x3EB4A92C4F46D75B), 0x1E5AE016, 202},
  { /*   734 */ UINT64_C (0x7D6952589E8DAEB6), 0x3CB5C02C, 202},
  { /*   735 */ UINT64_C (0xFAD2A4B13D1B5D6C), 0x796B8057, 202},
  { /*   736 */ UINT64_C (0x322A20F03F6BDF7C), 0x1848B345, 203},
  { /*   737 */ UINT64_C (0x645441E07ED7BEF8), 0x30916689, 203},
  { /*   738 */ UINT64_C (0xC8A883C0FDAF7DF0), 0x6122CD13, 203},
  { /*   739 */ UINT64_C (0x2821B3F365EFE5FC), 0xE03A2904, 204},
  { /*   740 */ UINT64_C (0x504367E6CBDFCBF9), 0xC0745207, 204},
  { /*   741 */ UINT64_C (0xA086CFCD97BF97F3), 0x80E8A40F, 204},
  { /*   742 */ UINT64_C (0x201AF65C518CB7FD), 0x802E8736, 205},
  { /*   743 */ UINT64_C (0x4035ECB8A3196FFB), 0x005D0E6C, 205},
  { /*   744 */ UINT64_C (0x806BD9714632DFF6), 0x00BA1CD9, 205},
  { /*   745 */ UINT64_C (0x19AF2B7D0E0A2CCA), 0xCCF205C5, 206},
  { /*   746 */ UINT64_C (0x335E56FA1C145995), 0x99E40B8A, 206},
  { /*   747 */ UINT64_C (0x66BCADF43828B32B), 0x33C81714, 206},
  { /*   748 */ UINT64_C (0xCD795BE870516656), 0x67902E27, 206},
  { /*   749 */ UINT64_C (0x29184594E3437ADE), 0x14B66FA1, 207},
  { /*   750 */ UINT64_C (0x52308B29C686F5BC), 0x296CDF43, 207},
  { /*   751 */ UINT64_C (0xA46116538D0DEB78), 0x52D9BE86, 207},
  { /*   752 */ UINT64_C (0x20E037AA4F692F18), 0x1091F2E8, 208},
  { /*   753 */ UINT64_C (0x41C06F549ED25E30), 0x2123E5CF, 208},
  { /*   754 */ UINT64_C (0x8380DEA93DA4BC60), 0x4247CB9E, 208},
  { /*   755 */ UINT64_C (0x1A4CF9550C5425AC), 0xDA0E5BEC, 209},
  { /*   756 */ UINT64_C (0x3499F2AA18A84B59), 0xB41CB7D9, 209},
  { /*   757 */ UINT64_C (0x6933E554315096B3), 0x68396FB2, 209},
  { /*   758 */ UINT64_C (0xD267CAA862A12D66), 0xD072DF64, 209},
  { /*   759 */ UINT64_C (0x2A14C221AD536F7A), 0xF67D5FE1, 210},
  { /*   760 */ UINT64_C (0x542984435AA6DEF5), 0xECFABFC2, 210},
  { /*   761 */ UINT64_C (0xA8530886B54DBDEB), 0xD9F57F83, 210},
  { /*   762 */ UINT64_C (0x21AA34E7BDDC592F), 0x2B977FE7, 211},
  { /*   763 */ UINT64_C (0x435469CF7BB8B25E), 0x572EFFCE, 211},
  { /*   764 */ UINT64_C (0x86A8D39EF77164BC), 0xAE5DFF9C, 211},
  { /*   765 */ UINT64_C (0x1AEE90B964B04758), 0xEFAC6652, 212},
  { /*   766 */ UINT64_C (0x35DD2172C9608EB1), 0xDF58CCA5, 212},
  { /*   767 */ UINT64_C (0x6BBA42E592C11D63), 0xBEB1994A, 212},
  { /*   768 */ UINT64_C (0xD77485CB25823AC7), 0x7D633293, 212},
  { /*   769 */ UINT64_C (0x2B174DF56DE6D88E), 0x4C470A1D, 213},
  { /*   770 */ UINT64_C (0x562E9BEADBCDB11C), 0x988E143B, 213},
  { /*   771 */ UINT64_C (0xAC5D37D5B79B6239), 0x311C2876, 213},
  { /*   772 */ UINT64_C (0x22790B2ABE5246D8), 0x3D05A1B1, 214},
  { /*   773 */ UINT64_C (0x44F216557CA48DB0), 0x7A0B4362, 214},
  { /*   774 */ UINT64_C (0x89E42CAAF9491B60), 0xF41686C5, 214},
  { /*   775 */ UINT64_C (0x1B9408EEFEA838AC), 0xFD9E1AF4, 215},
  { /*   776 */ UINT64_C (0x372811DDFD507159), 0xFB3C35E8, 215},
  { /*   777 */ UINT64_C (0x6E5023BBFAA0E2B3), 0xF6786BD0, 215},
  { /*   778 */ UINT64_C (0xDCA04777F541C567), 0xECF0D7A1, 215},
  { /*   779 */ UINT64_C (0x2C200E4B310D277B), 0x2F635E53, 216},
  { /*   780 */ UINT64_C (0x58401C96621A4EF6), 0x5EC6BCA7, 216},
  { /*   781 */ UINT64_C (0xB080392CC4349DEC), 0xBD8D794E, 216},
  { /*   782 */ UINT64_C (0x234CD83C273DB92F), 0x591C4B76, 217},
  { /*   783 */ UINT64_C (0x4699B0784E7B725E), 0xB23896EC, 217},
  { /*   784 */ UINT64_C (0x8D3360F09CF6E4BD), 0x64712DD8, 217},
  { /*   785 */ UINT64_C (0x1C3D79C9B8FE2DBF), 0x7A7D092B, 218},
  { /*   786 */ UINT64_C (0x387AF39371FC5B7E), 0xF4FA1256, 218},
  { /*   787 */ UINT64_C (0x70F5E726E3F8B6FD), 0xE9F424AD, 218},
  { /*   788 */ UINT64_C (0xE1EBCE4DC7F16DFB), 0xD3E84959, 218},
  { /*   789 */ UINT64_C (0x2D2F2942C196AF98), 0xC3FB41DF, 219},
  { /*   790 */ UINT64_C (0x5A5E5285832D5F31), 0x87F683BD, 219},
  { /*   791 */ UINT64_C (0xB4BCA50B065ABE63), 0x0FED077A, 219},
  { /*   792 */ UINT64_C (0x2425BA9BCE122613), 0xCFFC34B2, 220},
  { /*   793 */ UINT64_C (0x484B75379C244C27), 0x9FF86964, 220},
  { /*   794 */ UINT64_C (0x9096EA6F3848984F), 0x3FF0D2C8, 220},
  { /*   795 */ UINT64_C (0x1CEAFBAFD80E84DC), 0xA6635D5B, 221},
  { /*   796 */ UINT64_C (0x39D5F75FB01D09B9), 0x4CC6BAB7, 221},
  { /*   797 */ UINT64_C (0x73ABEEBF603A1372), 0x998D756D, 221},
  { /*   798 */ UINT64_C (0xE757DD7EC07426E5), 0x331AEADA, 221},
  { /*   799 */ UINT64_C (0x2E44C5E6267DA161), 0x0A38955F, 222},
  { /*   800 */ UINT64_C (0x5C898BCC4CFB42C2), 0x14712ABE, 222},
  { /*   801 */ UINT64_C (0xB913179899F68584), 0x28E2557B, 222},
  { /*   802 */ UINT64_C (0x2503D184EB97B44D), 0xA1C6DDE5, 223},
  { /*   803 */ UINT64_C (0x4A07A309D72F689B), 0x438DBBCB, 223},
  { /*   804 */ UINT64_C (0x940F4613AE5ED136), 0x871B7796, 223},
  { /*   805 */ UINT64_C (0x1D9CA79D894629D7), 0xB49F17EB, 224},
  { /*   806 */ UINT64_C (0x3B394F3B128C53AF), 0x693E2FD6, 224},
  { /*   807 */ UINT64_C (0x76729E762518A75E), 0xD27C5FAB, 224},
  { /*   808 */ UINT64_C (0xECE53CEC4A314EBD), 0xA4F8BF56, 224},
  { /*   809 */ UINT64_C (0x2F610C2F4209DC8C), 0x5431BFDE, 225},
  { /*   810 */ UINT64_C (0x5EC2185E8413B918), 0xA8637FBC, 225},
  { /*   811 */ UINT64_C (0xBD8430BD08277231), 0x50C6FF78, 225},
  { /*   812 */ UINT64_C (0x25E73CF29B3B16D6), 0xA9C1664B, 226},
  { /*   813 */ UINT64_C (0x4BCE79E536762DAD), 0x5382CC96, 226},
  { /*   814 */ UINT64_C (0x979CF3CA6CEC5B5A), 0xA705992D, 226},
  { /*   815 */ UINT64_C (0x1E5297287C2F4578), 0x87CDEB6F, 227},
  { /*   816 */ UINT64_C (0x3CA52E50F85E8AF1), 0x0F9BD6DF, 227},
  { /*   817 */ UINT64_C (0x794A5CA1F0BD15E2), 0x1F37ADBE, 227},
  { /*   818 */ UINT64_C (0xF294B943E17A2BC4), 0x3E6F5B7B, 227},
  { /*   819 */ UINT64_C (0x3084250D937ED58D), 0xA616457F, 228},
  { /*   820 */ UINT64_C (0x61084A1B26FDAB1B), 0x4C2C8AFE, 228},
  { /*   821 */ UINT64_C (0xC21094364DFB5636), 0x985915FC, 228},
  { /*   822 */ UINT64_C (0x26D01DA475FF113E), 0x1E783799, 229},
  { /*   823 */ UINT64_C (0x4DA03B48EBFE227C), 0x3CF06F32, 229},
  { /*   824 */ UINT64_C (0x9B407691D7FC44F8), 0x79E0DE63, 229},
  { /*   825 */ UINT64_C (0x1F0CE4839198DA98), 0x18602C7A, 230},
  { /*   826 */ UINT64_C (0x3E19C9072331B530), 0x30C058F5, 230},
  { /*   827 */ UINT64_C (0x7C33920E46636A60), 0x6180B1E9, 230},
  { /*   828 */ UINT64_C (0xF867241C8CC6D4C0), 0xC30163D2, 230},
  { /*   829 */ UINT64_C (0x31AE3A6C1C27C426), 0x8D66AD90, 231},
  { /*   830 */ UINT64_C (0x635C74D8384F884D), 0x1ACD5B21, 231},
  { /*   831 */ UINT64_C (0xC6B8E9B0709F109A), 0x359AB642, 231},
  { /*   832 */ UINT64_C (0x27BE952349B969B8), 0x711EF140, 232},
  { /*   833 */ UINT64_C (0x4F7D2A469372D370), 0xE23DE281, 232},
  { /*   834 */ UINT64_C (0x9EFA548D26E5A6E1), 0xC47BC501, 232},
  { /*   835 */ UINT64_C (0x1FCBAA82A1612160), 0x5A7F2767, 233},
  { /*   836 */ UINT64_C (0x3F97550542C242C0), 0xB4FE4ECD, 233},
  { /*   837 */ UINT64_C (0x7F2EAA0A85848581), 0x69FC9D9B, 233},
  { /*   838 */ UINT64_C (0xFE5D54150B090B02), 0xD3F93B35, 233},
  { /*   839 */ UINT64_C (0x32DF7737689B689A), 0x2A650BD7, 234},
  { /*   840 */ UINT64_C (0x65BEEE6ED136D134), 0x54CA17AF, 234},
  { /*   841 */ UINT64_C (0xCB7DDCDDA26DA268), 0xA9942F5E, 234},
  { /*   842 */ UINT64_C (0x28B2C5C5ED49207B), 0x551DA313, 235},
  { /*   843 */ UINT64_C (0x51658B8BDA9240F6), 0xAA3B4626, 235},
  { /*   844 */ UINT64_C (0xA2CB1717B52481ED), 0x54768C4B, 235},
  { /*   845 */ UINT64_C (0x208F049E576DB395), 0xDDB14F42, 236},
  { /*   846 */ UINT64_C (0x411E093CAEDB672B), 0xBB629E84, 236},
  { /*   847 */ UINT64_C (0x823C12795DB6CE57), 0x76C53D09, 236},
  { /*   848 */ UINT64_C (0x1A0C03B1DF8AF611), 0x7E27729B, 237},
  { /*   849 */ UINT64_C (0x34180763BF15EC22), 0xFC4EE537, 237},
  { /*   850 */ UINT64_C (0x68300EC77E2BD845), 0xF89DCA6D, 237},
  { /*   851 */ UINT64_C (0xD0601D8EFC57B08B), 0xF13B94DB, 237},
  { /*   852 */ UINT64_C (0x29ACD2B63277F01B), 0xFD0BEA92, 238},
  { /*   853 */ UINT64_C (0x5359A56C64EFE037), 0xFA17D524, 238},
  { /*   854 */ UINT64_C (0xA6B34AD8C9DFC06F), 0xF42FAA49, 238},
  { /*   855 */ UINT64_C (0x21570EF8285FF349), 0x973CBBA8, 239},
  { /*   856 */ UINT64_C (0x42AE1DF050BFE693), 0x2E797750, 239},
  { /*   857 */ UINT64_C (0x855C3BE0A17FCD26), 0x5CF2EEA1, 239},
  { /*   858 */ UINT64_C (0x1AAC0BF9B9E65C3A), 0xDF63C953, 240},
  { /*   859 */ UINT64_C (0x355817F373CCB875), 0xBEC792A7, 240},
  { /*   860 */ UINT64_C (0x6AB02FE6E79970EB), 0x7D8F254D, 240},
  { /*   861 */ UINT64_C (0xD5605FCDCF32E1D6), 0xFB1E4A9B, 240},
  { /*   862 */ UINT64_C (0x2AACDFF5F63D605E), 0x3239421F, 241},
  { /*   863 */ UINT64_C (0x5559BFEBEC7AC0BC), 0x6472843E, 241},
  { /*   864 */ UINT64_C (0xAAB37FD7D8F58178), 0xC8E5087C, 241},
  { /*   865 */ UINT64_C (0x2223E65E5E97804B), 0x5B6101B2, 242},
  { /*   866 */ UINT64_C (0x4447CCBCBD2F0096), 0xB6C20365, 242},
  { /*   867 */ UINT64_C (0x888F99797A5E012D), 0x6D8406C9, 242},
  { /*   868 */ UINT64_C (0x1B4FEB7EB212CD09), 0x15E7348F, 243},
  { /*   869 */ UINT64_C (0x369FD6FD64259A12), 0x2BCE691D, 243},
  { /*   870 */ UINT64_C (0x6D3FADFAC84B3424), 0x579CD23B, 243},
  { /*   871 */ UINT64_C (0xDA7F5BF590966848), 0xAF39A475, 243},
  { /*   872 */ UINT64_C (0x2BB31264501E14DB), 0x563EBA7E, 244},
  { /*   873 */ UINT64_C (0x576624C8A03C29B6), 0xAC7D74FC, 244},
  { /*   874 */ UINT64_C (0xAECC49914078536D), 0x58FAE9F7, 244},
  { /*   875 */ UINT64_C (0x22F5A850401810AF), 0x78322ECB, 245},
  { /*   876 */ UINT64_C (0x45EB50A08030215E), 0xF0645D96, 245},
  { /*   877 */ UINT64_C (0x8BD6A141006042BD), 0xE0C8BB2C, 245},
  { /*   878 */ UINT64_C (0x1BF7B9D9CCE00D59), 0x2CF4F23C, 246},
  { /*   879 */ UINT64_C (0x37EF73B399C01AB2), 0x59E9E478, 246},
  { /*   880 */ UINT64_C (0x6FDEE76733803564), 0xB3D3C8F0, 246},
  { /*   881 */ UINT64_C (0xDFBDCECE67006AC9), 0x67A791E1, 246},
  { /*   882 */ UINT64_C (0x2CBF8FC2E1667BC1), 0xE187E9FA, 247},
  { /*   883 */ UINT64_C (0x597F1F85C2CCF783), 0xC30FD3F3, 247},
  { /*   884 */ UINT64_C (0xB2FE3F0B8599EF07), 0x861FA7E7, 247},
  { /*   885 */ UINT64_C (0x23CC73024DEB9634), 0xB46CBB2E, 248},
  { /*   886 */ UINT64_C (0x4798E6049BD72C69), 0x68D9765C, 248},
  { /*   887 */ UINT64_C (0x8F31CC0937AE58D2), 0xD1B2ECB9, 248},
  { /*   888 */ UINT64_C (0x1CA38F350B22DE90), 0x9056FC25, 249},
  { /*   889 */ UINT64_C (0x39471E6A1645BD21), 0x20ADF84A, 249},
  { /*   890 */ UINT64_C (0x728E3CD42C8B7A42), 0x415BF094, 249},
  { /*   891 */ UINT64_C (0xE51C79A85916F484), 0x82B7E128, 249},
  { /*   892 */ UINT64_C (0x2DD27EBB4504974D), 0xB3BE603B, 250},
  { /*   893 */ UINT64_C (0x5BA4FD768A092E9B), 0x677CC076, 250},
  { /*   894 */ UINT64_C (0xB749FAED14125D36), 0xCEF980EC, 250},
  { /*   895 */ UINT64_C (0x24A865629D9D45D7), 0xC2FEB362, 251},
  { /*   896 */ UINT64_C (0x4950CAC53B3A8BAF), 0x85FD66C5, 251},
  { /*   897 */ UINT64_C (0x92A1958A7675175F), 0x0BFACD8A, 251},
  { /*   898 */ UINT64_C (0x1D53844EE47DD179), 0x68CBC2B5, 252},
  { /*   899 */ UINT64_C (0x3AA7089DC8FBA2F2), 0xD197856A, 252},
  { /*   900 */ UINT64_C (0x754E113B91F745E5), 0xA32F0AD5, 252},
  { /*   901 */ UINT64_C (0xEA9C227723EE8BCB), 0x465E15A9, 252},
  { /*   902 */ UINT64_C (0x2EEC06E4A0C94F28), 0xA7AC6ABB, 253},
  { /*   903 */ UINT64_C (0x5DD80DC941929E51), 0x4F58D577, 253},
  { /*   904 */ UINT64_C (0xBBB01B9283253CA2), 0x9EB1AAEE, 253},
  { /*   905 */ UINT64_C (0x25899F1D4D6DD8ED), 0x52F05563, 254},
  { /*   906 */ UINT64_C (0x4B133E3A9ADBB1DA), 0xA5E0AAC6, 254},
  { /*   907 */ UINT64_C (0x96267C7535B763B5), 0x4BC1558B, 254},
  { /*   908 */ UINT64_C (0x1E07B27DD78B13F1), 0x0F26AAB5, 255},
  { /*   909 */ UINT64_C (0x3C0F64FBAF1627E2), 0x1E4D556B, 255},
  { /*   910 */ UINT64_C (0x781EC9F75E2C4FC4), 0x3C9AAAD6, 255},
  { /*   911 */ UINT64_C (0xF03D93EEBC589F88), 0x793555AB, 255},
  { /*   912 */ UINT64_C (0x300C50C958DE864E), 0x7EA44455, 256},
  { /*   913 */ UINT64_C (0x6018A192B1BD0C9C), 0xFD4888AB, 256},
  { /*   914 */ UINT64_C (0xC0314325637A1939), 0xFA911156, 256},
  { /*   915 */ UINT64_C (0x267040A113E5383E), 0xCBB69D11, 257},
  { /*   916 */ UINT64_C (0x4CE0814227CA707D), 0x976D3A22, 257},
  { /*   917 */ UINT64_C (0x99C102844F94E0FB), 0x2EDA7445, 257},
  { /*   918 */ UINT64_C (0x1EC033B40FEA9365), 0x6FC54A74, 258},
  { /*   919 */ UINT64_C (0x3D8067681FD526CA), 0xDF8A94E8, 258},
  { /*   920 */ UINT64_C (0x7B00CED03FAA4D95), 0xBF1529D1, 258},
  { /*   921 */ UINT64_C (0xF6019DA07F549B2B), 0x7E2A53A1, 258},
  { /*   922 */ UINT64_C (0x313385ECE6441F08), 0xB2D543ED, 259},
  { /*   923 */ UINT64_C (0x62670BD9CC883E11), 0x65AA87DA, 259},
  { /*   924 */ UINT64_C (0xC4CE17B399107C22), 0xCB550FB4, 259},
  { /*   925 */ UINT64_C (0x275C6B23EB69B26D), 0x5BDDCFF1, 260},
  { /*   926 */ UINT64_C (0x4EB8D647D6D364DA), 0xB7BB9FE2, 260},
  { /*   927 */ UINT64_C (0x9D71AC8FADA6C9B5), 0x6F773FC3, 260},
  { /*   928 */ UINT64_C (0x1F7D228322BAF524), 0x497E3FF4, 261},
  { /*   929 */ UINT64_C (0x3EFA45064575EA48), 0x92FC7FE8, 261},
  { /*   930 */ UINT64_C (0x7DF48A0C8AEBD491), 0x25F8FFD0, 261},
  { /*   931 */ UINT64_C (0xFBE9141915D7A922), 0x4BF1FF9F, 261},
  { /*   932 */ UINT64_C (0x3261D0D1D12B21D3), 0xA8C9FFED, 262},
  { /*   933 */ UINT64_C (0x64C3A1A3A25643A7), 0x5193FFD9, 262},
  { /*   934 */ UINT64_C (0xC987434744AC874E), 0xA327FFB2, 262},
  { /*   935 */ UINT64_C (0x284E40A7DA88E7DC), 0x8707FFF0, 263},
  { /*   936 */ UINT64_C (0x509C814FB511CFB9), 0x0E0FFFE1, 263},
  { /*   937 */ UINT64_C (0xA139029F6A239F72), 0x1C1FFFC2, 263},
  { /*   938 */ UINT64_C (0x203E9A1FE2071FE3), 0x9F39998D, 264},
  { /*   939 */ UINT64_C (0x407D343FC40E3FC7), 0x3E73331A, 264},
  { /*   940 */ UINT64_C (0x80FA687F881C7F8E), 0x7CE66635, 264},
  { /*   941 */ UINT64_C (0x19CBAE7FE805B31C), 0x7F6147A4, 265},
  { /*   942 */ UINT64_C (0x33975CFFD00B6638), 0xFEC28F48, 265},
  { /*   943 */ UINT64_C (0x672EB9FFA016CC71), 0xFD851E91, 265},
  { /*   944 */ UINT64_C (0xCE5D73FF402D98E3), 0xFB0A3D21, 265},
  { /*   945 */ UINT64_C (0x2945E3FFD9A2B82D), 0x989BA5D3, 266},
  { /*   946 */ UINT64_C (0x528BC7FFB345705B), 0x31374BA7, 266},
  { /*   947 */ UINT64_C (0xA5178FFF668AE0B6), 0x626E974E, 266},
  { /*   948 */ UINT64_C (0x2104B66647B56024), 0x7A161E43, 267},
  { /*   949 */ UINT64_C (0x42096CCC8F6AC048), 0xF42C3C85, 267},
  { /*   950 */ UINT64_C (0x8412D9991ED58091), 0xE858790B, 267},
  { /*   951 */ UINT64_C (0x1A6A2B85062AB350), 0x61AB4B69, 268},
  { /*   952 */ UINT64_C (0x34D4570A0C5566A0), 0xC35696D1, 268},
  { /*   953 */ UINT64_C (0x69A8AE1418AACD41), 0x86AD2DA2, 268},
  { /*   954 */ UINT64_C (0xD3515C2831559A83), 0x0D5A5B45, 268},
  { /*   955 */ UINT64_C (0x2A4378D4D6AAB880), 0x9C454574, 269},
  { /*   956 */ UINT64_C (0x5486F1A9AD557101), 0x388A8AE8, 269},
  { /*   957 */ UINT64_C (0xA90DE3535AAAE202), 0x711515D1, 269},
  { /*   958 */ UINT64_C (0x21CF93DD7888939A), 0x169DD12A, 270},
  { /*   959 */ UINT64_C (0x439F27BAF1112734), 0x2D3BA253, 270},
  { /*   960 */ UINT64_C (0x873E4F75E2224E68), 0x5A7744A7, 270},
  { /*   961 */ UINT64_C (0x1B0C764AC6D3A948), 0x1217DA88, 271},
  { /*   962 */ UINT64_C (0x3618EC958DA75290), 0x242FB510, 271},
  { /*   963 */ UINT64_C (0x6C31D92B1B4EA520), 0x485F6A1F, 271},
  { /*   964 */ UINT64_C (0xD863B256369D4A40), 0x90BED43E, 271},
  { /*   965 */ UINT64_C (0x2B4723AAD7B90ED9), 0xB68C90D9, 272},
  { /*   966 */ UINT64_C (0x568E4755AF721DB3), 0x6D1921B3, 272},
  { /*   967 */ UINT64_C (0xAD1C8EAB5EE43B66), 0xDA324365, 272},
  { /*   968 */ UINT64_C (0x229F4FBBDFC73F14), 0x920A0D7B, 273},
  { /*   969 */ UINT64_C (0x453E9F77BF8E7E29), 0x24141AF5, 273},
  { /*   970 */ UINT64_C (0x8A7D3EEF7F1CFC52), 0x482835EA, 273},
  { /*   971 */ UINT64_C (0x1BB2A62FE638FF43), 0xA8080AC8, 274},
  { /*   972 */ UINT64_C (0x37654C5FCC71FE87), 0x50101591, 274},
  { /*   973 */ UINT64_C (0x6ECA98BF98E3FD0E), 0xA0202B22, 274},
  { /*   974 */ UINT64_C (0xDD95317F31C7FA1D), 0x40405644, 274},
  { /*   975 */ UINT64_C (0x2C5109E63D27FED2), 0xA6734474, 275},
  { /*   976 */ UINT64_C (0x58A213CC7A4FFDA5), 0x4CE688E8, 275},
  { /*   977 */ UINT64_C (0xB1442798F49FFB4A), 0x99CD11D0, 275},
  { /*   978 */ UINT64_C (0x237407EB641FFF0E), 0xEB8F69F6, 276},
  { /*   979 */ UINT64_C (0x46E80FD6C83FFE1D), 0xD71ED3ED, 276},
  { /*   980 */ UINT64_C (0x8DD01FAD907FFC3B), 0xAE3DA7D9, 276},
  { /*   981 */ UINT64_C (0x1C5CD322B67FFF3F), 0x22D92192, 277},
  { /*   982 */ UINT64_C (0x38B9A6456CFFFE7E), 0x45B24324, 277},
  { /*   983 */ UINT64_C (0x71734C8AD9FFFCFC), 0x8B648648, 277},
  { /*   984 */ UINT64_C (0xE2E69915B3FFF9F9), 0x16C90C8F, 277},
  { /*   985 */ UINT64_C (0x2D6151D123FFFECB), 0x6AF50283, 278},
  { /*   986 */ UINT64_C (0x5AC2A3A247FFFD96), 0xD5EA0506, 278},
  { /*   987 */ UINT64_C (0xB58547448FFFFB2D), 0xABD40A0C, 278},
  { /*   988 */ UINT64_C (0x244DDB0DB666656F), 0x88C40202, 279},
  { /*   989 */ UINT64_C (0x489BB61B6CCCCADF), 0x11880405, 279},
  { /*   990 */ UINT64_C (0x91376C36D99995BE), 0x2310080A, 279},
  { /*   991 */ UINT64_C (0x1D0B15A491EB8459), 0x3A366802, 280},
  { /*   992 */ UINT64_C (0x3A162B4923D708B2), 0x746CD004, 280},
  { /*   993 */ UINT64_C (0x742C569247AE1164), 0xE8D9A008, 280},
  { /*   994 */ UINT64_C (0xE858AD248F5C22C9), 0xD1B34010, 280},
  { /*   995 */ UINT64_C (0x2E7822A0E978D3C1), 0xF6BD7336, 281},
  { /*   996 */ UINT64_C (0x5CF04541D2F1A783), 0xED7AE66D, 281},
  { /*   997 */ UINT64_C (0xB9E08A83A5E34F07), 0xDAF5CCD9, 281},
  { /*   998 */ UINT64_C (0x252CE880BAC70FCE), 0x5EFDF5C5, 282},
  { /*   999 */ UINT64_C (0x4A59D101758E1F9C), 0xBDFBEB8A, 282},
  { /*  1000 */ UINT64_C (0x94B3A202EB1C3F39), 0x7BF7D714, 282},
  { /*  1001 */ UINT64_C (0x1DBD86CD6238D971), 0xE597F7D1, 283},
  { /*  1002 */ UINT64_C (0x3B7B0D9AC471B2E3), 0xCB2FEFA2, 283},
  { /*  1003 */ UINT64_C (0x76F61B3588E365C7), 0x965FDF43, 283},
  { /*  1004 */ UINT64_C (0xEDEC366B11C6CB8F), 0x2CBFBE87, 283},
  { /*  1005 */ UINT64_C (0x2F95A47BD05AF583), 0x08F3261B, 284},
  { /*  1006 */ UINT64_C (0x5F2B48F7A0B5EB06), 0x11E64C36, 284},
  { /*  1007 */ UINT64_C (0xBE5691EF416BD60C), 0x23CC986C, 284},
  { /*  1008 */ UINT64_C (0x261150630D159135), 0xA0C284E2, 285},
  { /*  1009 */ UINT64_C (0x4C22A0C61A2B226B), 0x418509C5, 285},
  { /*  1010 */ UINT64_C (0x9845418C345644D6), 0x830A1389, 285},
  { /*  1011 */ UINT64_C (0x1E74404F3DAADA91), 0x4D686A4F, 286},
  { /*  1012 */ UINT64_C (0x3CE8809E7B55B522), 0x9AD0D49D, 286},
  { /*  1013 */ UINT64_C (0x79D1013CF6AB6A45), 0x35A1A93B, 286},
  { /*  1014 */ UINT64_C (0xF3A20279ED56D48A), 0x6B435275, 286},
  { /*  1015 */ UINT64_C (0x30BA007EC9115DB5), 0x48A7107E, 287},
  { /*  1016 */ UINT64_C (0x617400FD9222BB6A), 0x914E20FC, 287},
  { /*  1017 */ UINT64_C (0xC2E801FB244576D5), 0x229C41F8, 287},
  { /*  1018 */ UINT64_C (0x26FB3398A0DAB15D), 0xD3B8D9FE, 288},
  { /*  1019 */ UINT64_C (0x4DF6673141B562BB), 0xA771B3FD, 288},
  { /*  1020 */ UINT64_C (0x9BECCE62836AC577), 0x4EE367F9, 288},
  { /*  1021 */ UINT64_C (0x1F2F5C7A1A488DE4), 0xA960AE65, 289},
  { /*  1022 */ UINT64_C (0x3E5EB8F434911BC9), 0x52C15CCA, 289},
  { /*  1023 */ UINT64_C (0x7CBD71E869223792), 0xA582B994, 289},
  { /*  1024 */ UINT64_C (0xF97AE3D0D2446F25), 0x4B057328, 289},
  { /*  1025 */ UINT64_C (0x31E560C35D40E307), 0x75677D6E, 290},
};

static const struct
{
  uint64_t mul;
  int32_t exp;
} fpowers2[] = {
  { /*  -151 */ UINT64_C (0x5190F96B91344AE4), -57},
  { /*  -150 */ UINT64_C (0xA321F2D7226895C8), -57},
  { /*  -149 */ UINT64_C (0x20A063C4A07B5128), -56},
  { /*  -148 */ UINT64_C (0x4140C78940F6A250), -56},
  { /*  -147 */ UINT64_C (0x82818F1281ED44A0), -56},
  { /*  -146 */ UINT64_C (0x1A19E96A19FC40ED), -55},
  { /*  -145 */ UINT64_C (0x3433D2D433F881D9), -55},
  { /*  -144 */ UINT64_C (0x6867A5A867F103B3), -55},
  { /*  -143 */ UINT64_C (0xD0CF4B50CFE20766), -55},
  { /*  -142 */ UINT64_C (0x29C30F1029939B14), -54},
  { /*  -141 */ UINT64_C (0x53861E2053273629), -54},
  { /*  -140 */ UINT64_C (0xA70C3C40A64E6C52), -54},
  { /*  -139 */ UINT64_C (0x2168D8D9BADC7C10), -53},
  { /*  -138 */ UINT64_C (0x42D1B1B375B8F821), -53},
  { /*  -137 */ UINT64_C (0x85A36366EB71F041), -53},
  { /*  -136 */ UINT64_C (0x1ABA4714957D300D), -52},
  { /*  -135 */ UINT64_C (0x35748E292AFA601A), -52},
  { /*  -134 */ UINT64_C (0x6AE91C5255F4C034), -52},
  { /*  -133 */ UINT64_C (0xD5D238A4ABE98068), -52},
  { /*  -132 */ UINT64_C (0x2AC3A4EDBBFB8015), -51},
  { /*  -131 */ UINT64_C (0x558749DB77F7002A), -51},
  { /*  -130 */ UINT64_C (0xAB0E93B6EFEE0054), -51},
  { /*  -129 */ UINT64_C (0x22361D8AFCC93344), -50},
  { /*  -128 */ UINT64_C (0x446C3B15F9926688), -50},
  { /*  -127 */ UINT64_C (0x88D8762BF324CD10), -50},
  { /*  -126 */ UINT64_C (0x1B5E7E08CA3A8F6A), -49},
  { /*  -125 */ UINT64_C (0x36BCFC1194751ED3), -49},
  { /*  -124 */ UINT64_C (0x6D79F82328EA3DA6), -49},
  { /*  -123 */ UINT64_C (0xDAF3F04651D47B4C), -49},
  { /*  -122 */ UINT64_C (0x2BCA63414390E576), -48},
  { /*  -121 */ UINT64_C (0x5794C6828721CAEB), -48},
  { /*  -120 */ UINT64_C (0xAF298D050E4395D7), -48},
  { /*  -119 */ UINT64_C (0x23084F676940B791), -47},
  { /*  -118 */ UINT64_C (0x46109ECED2816F23), -47},
  { /*  -117 */ UINT64_C (0x8C213D9DA502DE45), -47},
  { /*  -116 */ UINT64_C (0x1C06A5EC5433C60E), -46},
  { /*  -115 */ UINT64_C (0x380D4BD8A8678C1C), -46},
  { /*  -114 */ UINT64_C (0x701A97B150CF1837), -46},
  { /*  -113 */ UINT64_C (0xE0352F62A19E306F), -46},
  { /*  -112 */ UINT64_C (0x2CD76FE086B93CE3), -45},
  { /*  -111 */ UINT64_C (0x59AEDFC10D7279C6), -45},
  { /*  -110 */ UINT64_C (0xB35DBF821AE4F38C), -45},
  { /*  -109 */ UINT64_C (0x23DF8CB39EFA971C), -44},
  { /*  -108 */ UINT64_C (0x47BF19673DF52E38), -44},
  { /*  -107 */ UINT64_C (0x8F7E32CE7BEA5C70), -44},
  { /*  -106 */ UINT64_C (0x1CB2D6F618C878E3), -43},
  { /*  -105 */ UINT64_C (0x3965ADEC3190F1C6), -43},
  { /*  -104 */ UINT64_C (0x72CB5BD86321E38D), -43},
  { /*  -103 */ UINT64_C (0xE596B7B0C643C719), -43},
  { /*  -102 */ UINT64_C (0x2DEAF189C140C16B), -42},
  { /*  -101 */ UINT64_C (0x5BD5E313828182D7), -42},
  { /*  -100 */ UINT64_C (0xB7ABC627050305AE), -42},
  { /*   -99 */ UINT64_C (0x24BBF46E3433CDF0), -41},
  { /*   -98 */ UINT64_C (0x4977E8DC68679BDF), -41},
  { /*   -97 */ UINT64_C (0x92EFD1B8D0CF37BE), -41},
  { /*   -96 */ UINT64_C (0x1D6329F1C35CA4C0), -40},
  { /*   -95 */ UINT64_C (0x3AC653E386B9497F), -40},
  { /*   -94 */ UINT64_C (0x758CA7C70D7292FF), -40},
  { /*   -93 */ UINT64_C (0x178287F49C4A1D66), -39},
  { /*   -92 */ UINT64_C (0x2F050FE938943ACC), -39},
  { /*   -91 */ UINT64_C (0x5E0A1FD271287599), -39},
  { /*   -90 */ UINT64_C (0xBC143FA4E250EB31), -39},
  { /*   -89 */ UINT64_C (0x259DA6542D43623D), -38},
  { /*   -88 */ UINT64_C (0x4B3B4CA85A86C47A), -38},
  { /*   -87 */ UINT64_C (0x96769950B50D88F4), -38},
  { /*   -86 */ UINT64_C (0x1E17B84357691B64), -37},
  { /*   -85 */ UINT64_C (0x3C2F7086AED236C8), -37},
  { /*   -84 */ UINT64_C (0x785EE10D5DA46D90), -37},
  { /*   -83 */ UINT64_C (0x1812F9CF7920E2B6), -36},
  { /*   -82 */ UINT64_C (0x3025F39EF241C56D), -36},
  { /*   -81 */ UINT64_C (0x604BE73DE4838ADA), -36},
  { /*   -80 */ UINT64_C (0xC097CE7BC90715B3), -36},
  { /*   -79 */ UINT64_C (0x2684C2E58E9B0457), -35},
  { /*   -78 */ UINT64_C (0x4D0985CB1D3608AE), -35},
  { /*   -77 */ UINT64_C (0x9A130B963A6C115C), -35},
  { /*   -76 */ UINT64_C (0x1ED09BEAD87C0379), -34},
  { /*   -75 */ UINT64_C (0x3DA137D5B0F806F2), -34},
  { /*   -74 */ UINT64_C (0x7B426FAB61F00DE3), -34},
  { /*   -73 */ UINT64_C (0x18A6E32246C99C61), -33},
  { /*   -72 */ UINT64_C (0x314DC6448D9338C1), -33},
  { /*   -71 */ UINT64_C (0x629B8C891B267183), -33},
  { /*   -70 */ UINT64_C (0xC5371912364CE305), -33},
  { /*   -69 */ UINT64_C (0x27716B6A0ADC2D67), -32},
  { /*   -68 */ UINT64_C (0x4EE2D6D415B85ACF), -32},
  { /*   -67 */ UINT64_C (0x9DC5ADA82B70B59E), -32},
  { /*   -66 */ UINT64_C (0x1F8DEF8808B02453), -31},
  { /*   -65 */ UINT64_C (0x3F1BDF10116048A6), -31},
  { /*   -64 */ UINT64_C (0x7E37BE2022C0914B), -31},
  { /*   -63 */ UINT64_C (0x193E5939A08CE9DC), -30},
  { /*   -62 */ UINT64_C (0x327CB2734119D3B8), -30},
  { /*   -61 */ UINT64_C (0x64F964E68233A76F), -30},
  { /*   -60 */ UINT64_C (0xC9F2C9CD04674EDF), -30},
  { /*   -59 */ UINT64_C (0x2863C1F5CDAE42F9), -29},
  { /*   -58 */ UINT64_C (0x50C783EB9B5C85F3), -29},
  { /*   -57 */ UINT64_C (0xA18F07D736B90BE5), -29},
  { /*   -56 */ UINT64_C (0x204FCE5E3E250261), -28},
  { /*   -55 */ UINT64_C (0x409F9CBC7C4A04C2), -28},
  { /*   -54 */ UINT64_C (0x813F3978F8940984), -28},
  { /*   -53 */ UINT64_C (0x19D971E4FE8401E7), -27},
  { /*   -52 */ UINT64_C (0x33B2E3C9FD0803CF), -27},
  { /*   -51 */ UINT64_C (0x6765C793FA10079D), -27},
  { /*   -50 */ UINT64_C (0xCECB8F27F4200F3A), -27},
  { /*   -49 */ UINT64_C (0x295BE96E64066972), -26},
  { /*   -48 */ UINT64_C (0x52B7D2DCC80CD2E4), -26},
  { /*   -47 */ UINT64_C (0xA56FA5B99019A5C8), -26},
  { /*   -46 */ UINT64_C (0x2116545850052128), -25},
  { /*   -45 */ UINT64_C (0x422CA8B0A00A4250), -25},
  { /*   -44 */ UINT64_C (0x84595161401484A0), -25},
  { /*   -43 */ UINT64_C (0x1A784379D99DB420), -24},
  { /*   -42 */ UINT64_C (0x34F086F3B33B6840), -24},
  { /*   -41 */ UINT64_C (0x69E10DE76676D080), -24},
  { /*   -40 */ UINT64_C (0xD3C21BCECCEDA100), -24},
  { /*   -39 */ UINT64_C (0x2A5A058FC295ED00), -23},
  { /*   -38 */ UINT64_C (0x54B40B1F852BDA00), -23},
  { /*   -37 */ UINT64_C (0xA968163F0A57B400), -23},
  { /*   -36 */ UINT64_C (0x21E19E0C9BAB2400), -22},
  { /*   -35 */ UINT64_C (0x43C33C1937564800), -22},
  { /*   -34 */ UINT64_C (0x878678326EAC9000), -22},
  { /*   -33 */ UINT64_C (0x1B1AE4D6E2EF5000), -21},
  { /*   -32 */ UINT64_C (0x3635C9ADC5DEA000), -21},
  { /*   -31 */ UINT64_C (0x6C6B935B8BBD4000), -21},
  { /*   -30 */ UINT64_C (0xD8D726B7177A8000), -21},
  { /*   -29 */ UINT64_C (0x2B5E3AF16B188000), -20},
  { /*   -28 */ UINT64_C (0x56BC75E2D6310000), -20},
  { /*   -27 */ UINT64_C (0xAD78EBC5AC620000), -20},
  { /*   -26 */ UINT64_C (0x22B1C8C1227A0000), -19},
  { /*   -25 */ UINT64_C (0x4563918244F40000), -19},
  { /*   -24 */ UINT64_C (0x8AC7230489E80000), -19},
  { /*   -23 */ UINT64_C (0x1BC16D674EC80000), -18},
  { /*   -22 */ UINT64_C (0x3782DACE9D900000), -18},
  { /*   -21 */ UINT64_C (0x6F05B59D3B200000), -18},
  { /*   -20 */ UINT64_C (0xDE0B6B3A76400000), -18},
  { /*   -19 */ UINT64_C (0x2C68AF0BB1400000), -17},
  { /*   -18 */ UINT64_C (0x58D15E1762800000), -17},
  { /*   -17 */ UINT64_C (0xB1A2BC2EC5000000), -17},
  { /*   -16 */ UINT64_C (0x2386F26FC1000000), -16},
  { /*   -15 */ UINT64_C (0x470DE4DF82000000), -16},
  { /*   -14 */ UINT64_C (0x8E1BC9BF04000000), -16},
  { /*   -13 */ UINT64_C (0x1C6BF52634000000), -15},
  { /*   -12 */ UINT64_C (0x38D7EA4C68000000), -15},
  { /*   -11 */ UINT64_C (0x71AFD498D0000000), -15},
  { /*   -10 */ UINT64_C (0xE35FA931A0000000), -15},
  { /*    -9 */ UINT64_C (0x2D79883D20000000), -14},
  { /*    -8 */ UINT64_C (0x5AF3107A40000000), -14},
  { /*    -7 */ UINT64_C (0xB5E620F480000000), -14},
  { /*    -6 */ UINT64_C (0x246139CA80000000), -13},
  { /*    -5 */ UINT64_C (0x48C2739500000000), -13},
  { /*    -4 */ UINT64_C (0x9184E72A00000000), -13},
  { /*    -3 */ UINT64_C (0x1D1A94A200000000), -12},
  { /*    -2 */ UINT64_C (0x3A35294400000000), -12},
  { /*    -1 */ UINT64_C (0x746A528800000000), -12},
  { /*     0 */ UINT64_C (0x174876E800000000), -11},
  { /*     1 */ UINT64_C (0x2E90EDD000000000), -11},
  { /*     2 */ UINT64_C (0x5D21DBA000000000), -11},
  { /*     3 */ UINT64_C (0xBA43B74000000000), -11},
  { /*     4 */ UINT64_C (0x2540BE4000000000), -10},
  { /*     5 */ UINT64_C (0x4A817C8000000000), -10},
  { /*     6 */ UINT64_C (0x9502F90000000000), -10},
  { /*     7 */ UINT64_C (0x1DCD650000000000), -9},
  { /*     8 */ UINT64_C (0x3B9ACA0000000000), -9},
  { /*     9 */ UINT64_C (0x7735940000000000), -9},
  { /*    10 */ UINT64_C (0x17D7840000000000), -8},
  { /*    11 */ UINT64_C (0x2FAF080000000000), -8},
  { /*    12 */ UINT64_C (0x5F5E100000000000), -8},
  { /*    13 */ UINT64_C (0xBEBC200000000000), -8},
  { /*    14 */ UINT64_C (0x2625A00000000000), -7},
  { /*    15 */ UINT64_C (0x4C4B400000000000), -7},
  { /*    16 */ UINT64_C (0x9896800000000000), -7},
  { /*    17 */ UINT64_C (0x1E84800000000000), -6},
  { /*    18 */ UINT64_C (0x3D09000000000000), -6},
  { /*    19 */ UINT64_C (0x7A12000000000000), -6},
  { /*    20 */ UINT64_C (0x186A000000000000), -5},
  { /*    21 */ UINT64_C (0x30D4000000000000), -5},
  { /*    22 */ UINT64_C (0x61A8000000000000), -5},
  { /*    23 */ UINT64_C (0xC350000000000000), -5},
  { /*    24 */ UINT64_C (0x2710000000000000), -4},
  { /*    25 */ UINT64_C (0x4E20000000000000), -4},
  { /*    26 */ UINT64_C (0x9C40000000000000), -4},
  { /*    27 */ UINT64_C (0x1F40000000000000), -3},
  { /*    28 */ UINT64_C (0x3E80000000000000), -3},
  { /*    29 */ UINT64_C (0x7D00000000000000), -3},
  { /*    30 */ UINT64_C (0x1900000000000000), -2},
  { /*    31 */ UINT64_C (0x3200000000000000), -2},
  { /*    32 */ UINT64_C (0x6400000000000000), -2},
  { /*    33 */ UINT64_C (0xC800000000000000), -2},
  { /*    34 */ UINT64_C (0x2800000000000000), -1},
  { /*    35 */ UINT64_C (0x5000000000000000), -1},
  { /*    36 */ UINT64_C (0xA000000000000000), -1},
  { /*    37 */ UINT64_C (0x2000000000000000), 0},
  { /*    38 */ UINT64_C (0x4000000000000000), 0},
  { /*    39 */ UINT64_C (0x8000000000000000), 0},
  { /*    40 */ UINT64_C (0x199999999999999A), 1},
  { /*    41 */ UINT64_C (0x3333333333333333), 1},
  { /*    42 */ UINT64_C (0x6666666666666666), 1},
  { /*    43 */ UINT64_C (0xCCCCCCCCCCCCCCCD), 1},
  { /*    44 */ UINT64_C (0x28F5C28F5C28F5C3), 2},
  { /*    45 */ UINT64_C (0x51EB851EB851EB85), 2},
  { /*    46 */ UINT64_C (0xA3D70A3D70A3D70A), 2},
  { /*    47 */ UINT64_C (0x20C49BA5E353F7CF), 3},
  { /*    48 */ UINT64_C (0x4189374BC6A7EF9E), 3},
  { /*    49 */ UINT64_C (0x83126E978D4FDF3B), 3},
  { /*    50 */ UINT64_C (0x1A36E2EB1C432CA5), 4},
  { /*    51 */ UINT64_C (0x346DC5D63886594B), 4},
  { /*    52 */ UINT64_C (0x68DB8BAC710CB296), 4},
  { /*    53 */ UINT64_C (0xD1B71758E219652C), 4},
  { /*    54 */ UINT64_C (0x29F16B11C6D1E109), 5},
  { /*    55 */ UINT64_C (0x53E2D6238DA3C212), 5},
  { /*    56 */ UINT64_C (0xA7C5AC471B478423), 5},
  { /*    57 */ UINT64_C (0x218DEF416BDB1A6D), 6},
  { /*    58 */ UINT64_C (0x431BDE82D7B634DB), 6},
  { /*    59 */ UINT64_C (0x8637BD05AF6C69B6), 6},
  { /*    60 */ UINT64_C (0x1AD7F29ABCAF4858), 7},
  { /*    61 */ UINT64_C (0x35AFE535795E90AF), 7},
  { /*    62 */ UINT64_C (0x6B5FCA6AF2BD215E), 7},
  { /*    63 */ UINT64_C (0xD6BF94D5E57A42BC), 7},
  { /*    64 */ UINT64_C (0x2AF31DC4611873BF), 8},
  { /*    65 */ UINT64_C (0x55E63B88C230E77E), 8},
  { /*    66 */ UINT64_C (0xABCC77118461CEFD), 8},
  { /*    67 */ UINT64_C (0x225C17D04DAD2966), 9},
  { /*    68 */ UINT64_C (0x44B82FA09B5A52CC), 9},
  { /*    69 */ UINT64_C (0x89705F4136B4A597), 9},
  { /*    70 */ UINT64_C (0x1B7CDFD9D7BDBAB8), 10},
  { /*    71 */ UINT64_C (0x36F9BFB3AF7B7570), 10},
  { /*    72 */ UINT64_C (0x6DF37F675EF6EADF), 10},
  { /*    73 */ UINT64_C (0xDBE6FECEBDEDD5BF), 10},
  { /*    74 */ UINT64_C (0x2BFAFFC2F2C92AC0), 11},
  { /*    75 */ UINT64_C (0x57F5FF85E592557F), 11},
  { /*    76 */ UINT64_C (0xAFEBFF0BCB24AAFF), 11},
  { /*    77 */ UINT64_C (0x232F33025BD42233), 12},
  { /*    78 */ UINT64_C (0x465E6604B7A84466), 12},
  { /*    79 */ UINT64_C (0x8CBCCC096F5088CC), 12},
  { /*    80 */ UINT64_C (0x1C25C268497681C2), 13},
  { /*    81 */ UINT64_C (0x384B84D092ED0385), 13},
  { /*    82 */ UINT64_C (0x709709A125DA070A), 13},
  { /*    83 */ UINT64_C (0xE12E13424BB40E13), 13},
  { /*    84 */ UINT64_C (0x2D09370D42573604), 14},
  { /*    85 */ UINT64_C (0x5A126E1A84AE6C08), 14},
  { /*    86 */ UINT64_C (0xB424DC35095CD80F), 14},
  { /*    87 */ UINT64_C (0x24075F3DCEAC2B36), 15},
  { /*    88 */ UINT64_C (0x480EBE7B9D58566D), 15},
  { /*    89 */ UINT64_C (0x901D7CF73AB0ACD9), 15},
  { /*    90 */ UINT64_C (0x1CD2B297D889BC2B), 16},
  { /*    91 */ UINT64_C (0x39A5652FB1137857), 16},
  { /*    92 */ UINT64_C (0x734ACA5F6226F0AE), 16},
  { /*    93 */ UINT64_C (0xE69594BEC44DE15B), 16},
  { /*    94 */ UINT64_C (0x2E1DEA8C8DA92D12), 17},
  { /*    95 */ UINT64_C (0x5C3BD5191B525A25), 17},
  { /*    96 */ UINT64_C (0xB877AA3236A4B449), 17},
  { /*    97 */ UINT64_C (0x24E4BBA3A4875742), 18},
  { /*    98 */ UINT64_C (0x49C97747490EAE84), 18},
  { /*    99 */ UINT64_C (0x9392EE8E921D5D07), 18},
  { /*   100 */ UINT64_C (0x1D83C94FB6D2AC35), 19},
  { /*   101 */ UINT64_C (0x3B07929F6DA55869), 19},
  { /*   102 */ UINT64_C (0x760F253EDB4AB0D3), 19},
  { /*   103 */ UINT64_C (0x179CA10C9242235D), 20},
  { /*   104 */ UINT64_C (0x2F394219248446BB), 20},
  { /*   105 */ UINT64_C (0x5E72843249088D75), 20},
  { /*   106 */ UINT64_C (0xBCE5086492111AEB), 20},
  { /*   107 */ UINT64_C (0x25C768141D369EFC), 21},
  { /*   108 */ UINT64_C (0x4B8ED0283A6D3DF7), 21},
  { /*   109 */ UINT64_C (0x971DA05074DA7BEF), 21},
  { /*   110 */ UINT64_C (0x1E392010175EE596), 22},
  { /*   111 */ UINT64_C (0x3C7240202EBDCB2C), 22},
  { /*   112 */ UINT64_C (0x78E480405D7B9659), 22},
  { /*   113 */ UINT64_C (0x182DB34012B25145), 23},
  { /*   114 */ UINT64_C (0x305B66802564A28A), 23},
  { /*   115 */ UINT64_C (0x60B6CD004AC94514), 23},
  { /*   116 */ UINT64_C (0xC16D9A0095928A27), 23},
  { /*   117 */ UINT64_C (0x26AF8533511D4ED5), 24},
  { /*   118 */ UINT64_C (0x4D5F0A66A23A9DA9), 24},
  { /*   119 */ UINT64_C (0x9ABE14CD44753B53), 24},
  { /*   120 */ UINT64_C (0x1EF2D0F5DA7DD8AA), 25},
  { /*   121 */ UINT64_C (0x3DE5A1EBB4FBB154), 25},
  { /*   122 */ UINT64_C (0x7BCB43D769F762A9), 25},
  { /*   123 */ UINT64_C (0x18C240C4AECB13BB), 26},
  { /*   124 */ UINT64_C (0x318481895D962777), 26},
  { /*   125 */ UINT64_C (0x63090312BB2C4EED), 26},
  { /*   126 */ UINT64_C (0xC612062576589DDB), 26},
  { /*   127 */ UINT64_C (0x279D346DE4781F92), 27},
  { /*   128 */ UINT64_C (0x4F3A68DBC8F03F24), 27},
  { /*   129 */ UINT64_C (0x9E74D1B791E07E48), 27},
};

#if 0
/* gcc -g -O3 -Wall b.c -o b -lmpfr -lgmp */
#include <stdio.h>
#include <mpfr.h>

#define	PRINT_VALUE	0

int
main (void)
{
  int i;
  int cmp;
  unsigned long value;
  unsigned long value1;
  unsigned long value2;
  long exp;
  mpfr_t p64;
  mpfr_t p32;
  mpfr_t t1;
  mpfr_t t2;
  mpfr_t t3;
  mpfr_t t4;
  unsigned long long v;

  mpfr_init2 (p64, 128);
  mpfr_init2 (p32, 128);
  mpfr_init2 (t1, 128);
  mpfr_init2 (t2, 128);
  mpfr_init2 (t3, 128);
  mpfr_init2 (t4, 128);
  mpfr_ui_pow_ui (p64, 2, 64, MPFR_RNDN);
  mpfr_ui_pow_ui (p32, 2, 32, MPFR_RNDN);

  printf ("static const struct\n");
  printf ("{\n");
  printf ("  uint64_t mul1;\n");
  printf ("  uint32_t mul2;\n");
  printf ("  int32_t exp;\n");
  printf ("} dpowers10[] = {\n");
  for (i = -362; i <= 309; i++) {
    /* t1 = powl (10.0, i); */
    mpfr_set_ui (t3, 10, MPFR_RNDN);
    mpfr_pow_si (t1, t3, i, MPFR_RNDN);
    if (i < 0) {
      /* t2 = log10l (1.0 / t1) / log10l (2.0); */
      mpfr_set_ui (t3, 1, MPFR_RNDN);
      mpfr_div (t3, t3, t1, MPFR_RNDN);
      mpfr_log10 (t2, t3, MPFR_RNDN);
      mpfr_set_ui (t3, 2, MPFR_RNDN);
      mpfr_log10 (t3, t3, MPFR_RNDN);
      mpfr_div (t2, t2, t3, MPFR_RNDN);
      /* t2 = floorl (t2 + 66); */
      mpfr_add_ui (t3, t2, 66, MPFR_RNDN);
      mpfr_rint_floor (t2, t3, MPFR_RNDN);
      /* while (powl (2.0, t2) * t1 >= p64)  */
      mpfr_set_ui (t3, 2, MPFR_RNDN);
      mpfr_pow (t3, t3, t2, MPFR_RNDN);
      mpfr_mul (t3, t3, t1, MPFR_RNDN);
      cmp = mpfr_cmp (t3, p64);
      while (cmp >= 0) {
	/* t2 -= 1.0; */
	mpfr_sub_ui (t2, t2, 1, MPFR_RNDN);
	mpfr_set_ui (t3, 2, MPFR_RNDN);
	mpfr_pow (t3, t3, t2, MPFR_RNDN);
	mpfr_mul (t3, t3, t1, MPFR_RNDN);
	cmp = mpfr_cmp (t3, p64);
      }
      mpfr_mul (t3, t3, p32, MPFR_RNDN);
      mpfr_round (t3, t3);
      mpfr_div (t4, t3, p32, MPFR_RNDN);
      mpfr_floor (t4, t4);
      value1 = mpfr_get_ui (t4, MPFR_RNDN);
      mpfr_mul (t4, t4, p32, MPFR_RNDN);
      mpfr_sub (t4, t3, t4, MPFR_RNDN);
      value2 = mpfr_get_ui (t4, MPFR_RNDN);
      exp = mpfr_get_si (t2, MPFR_RNDN);
#if PRINT_VALUE
      mpfr_printf (" /* %40.16Rf */ ", t3);
#endif
      printf ("  { /* %5d */ UINT64_C (0x%016lX), 0x%08lX, %4ld },\n", i,
	      value1, value2, -exp);
    }
    else {
      /* t2 = log10l (t1) / log10l (2.0); */
      mpfr_log10 (t2, t1, MPFR_RNDN);
      mpfr_set_ui (t3, 2, MPFR_RNDN);
      mpfr_log10 (t3, t3, MPFR_RNDN);
      mpfr_div (t2, t2, t3, MPFR_RNDN);
      /* t2 = floorl (t2 - 66); */
      mpfr_sub_ui (t3, t2, 66, MPFR_RNDN);
      mpfr_rint_floor (t2, t3, MPFR_RNDN);
      /* while (t1 / powl (2.0, t2) >= p64) */
      mpfr_set_ui (t3, 2, MPFR_RNDN);
      mpfr_pow (t3, t3, t2, MPFR_RNDN);
      mpfr_div (t3, t1, t3, MPFR_RNDN);
      cmp = mpfr_cmp (t3, p64);
      while (cmp >= 0) {
	/* t2 += 1.0; */
	mpfr_add_ui (t2, t2, 1, MPFR_RNDN);
	mpfr_set_ui (t3, 2, MPFR_RNDN);
	mpfr_pow (t3, t3, t2, MPFR_RNDN);
	mpfr_div (t3, t1, t3, MPFR_RNDN);
	cmp = mpfr_cmp (t3, p64);
      }
      mpfr_mul (t3, t3, p32, MPFR_RNDN);
      mpfr_round (t3, t3);
      mpfr_div (t4, t3, p32, MPFR_RNDN);
      mpfr_floor (t4, t4);
      value1 = mpfr_get_ui (t4, MPFR_RNDN);
      mpfr_mul (t4, t4, p32, MPFR_RNDN);
      mpfr_sub (t4, t3, t4, MPFR_RNDN);
      value2 = mpfr_get_ui (t4, MPFR_RNDN);
      exp = mpfr_get_si (t2, MPFR_RNDN);
#if PRINT_VALUE
      mpfr_printf (" /* %40.16Rf */ ", t3);
#endif
      printf ("  { /* %5d */ UINT64_C (0x%016lX), 0x%08lX, %4ld },\n", i,
	      value1, value2, exp);
    }
  }
  printf ("};\n\n");
  printf ("static const struct\n");
  printf ("{\n");
  printf ("  uint64_t mul;\n");
  printf ("  int32_t exp;\n");
  printf ("} fpowers10[] = {\n");
  for (i = -64; i <= 39; i++) {
    /* t1 = powl (10.0, i); */
    mpfr_set_ui (t3, 10, MPFR_RNDN);
    mpfr_pow_si (t1, t3, i, MPFR_RNDN);
    if (i < 0) {
      /* t2 = log10l (1.0 / t1) / log10l (2.0); */
      mpfr_set_ui (t3, 1, MPFR_RNDN);
      mpfr_div (t3, t3, t1, MPFR_RNDN);
      mpfr_log10 (t2, t3, MPFR_RNDN);
      mpfr_set_ui (t3, 2, MPFR_RNDN);
      mpfr_log10 (t3, t3, MPFR_RNDN);
      mpfr_div (t2, t2, t3, MPFR_RNDN);
      /* t2 = floorl (t2 + 34); */
      mpfr_add_ui (t3, t2, 34, MPFR_RNDN);
      mpfr_rint_floor (t2, t3, MPFR_RNDN);
      /* while (powl (2.0, t2) * t1 >= p32)  */
      mpfr_set_ui (t3, 2, MPFR_RNDN);
      mpfr_pow (t3, t3, t2, MPFR_RNDN);
      mpfr_mul (t3, t3, t1, MPFR_RNDN);
      cmp = mpfr_cmp (t3, p32);
      while (cmp >= 0) {
	/* t2 -= 1.0; */
	mpfr_sub_ui (t2, t2, 1, MPFR_RNDN);
	mpfr_set_ui (t3, 2, MPFR_RNDN);
	mpfr_pow (t3, t3, t2, MPFR_RNDN);
	mpfr_mul (t3, t3, t1, MPFR_RNDN);
	cmp = mpfr_cmp (t3, p32);
      }
      mpfr_mul (t3, t3, p32, MPFR_RNDN);
      mpfr_round (t3, t3);
      value = mpfr_get_ui (t3, MPFR_RNDN);
      exp = mpfr_get_si (t2, MPFR_RNDN);
#if PRINT_VALUE
      mpfr_printf (" /* %40.16Rf */ ", t3);
#endif
      printf ("  { /* %5d */ UINT64_C (0x%016lX), %4ld },\n", i, value, -exp);
    }
    else {
      /* t2 = log10l (t1) / log10l (2.0); */
      mpfr_log10 (t2, t1, MPFR_RNDN);
      mpfr_set_ui (t3, 2, MPFR_RNDN);
      mpfr_log10 (t3, t3, MPFR_RNDN);
      mpfr_div (t2, t2, t3, MPFR_RNDN);
      /* t2 = floorl (t2 - 34); */
      mpfr_sub_ui (t3, t2, 34, MPFR_RNDN);
      mpfr_rint_floor (t2, t3, MPFR_RNDN);
      /* while (t1 / powl (2.0, t2) >= p32) */
      mpfr_set_ui (t3, 2, MPFR_RNDN);
      mpfr_pow (t3, t3, t2, MPFR_RNDN);
      mpfr_div (t3, t1, t3, MPFR_RNDN);
      cmp = mpfr_cmp (t3, p32);
      while (cmp >= 0) {
	/* t2 += 1.0; */
	mpfr_add_ui (t2, t2, 1, MPFR_RNDN);
	mpfr_set_ui (t3, 2, MPFR_RNDN);
	mpfr_pow (t3, t3, t2, MPFR_RNDN);
	mpfr_div (t3, t1, t3, MPFR_RNDN);
	cmp = mpfr_cmp (t3, p32);
      }
      mpfr_mul (t3, t3, p32, MPFR_RNDN);
      mpfr_round (t3, t3);
      value = mpfr_get_ui (t3, MPFR_RNDN);
      exp = mpfr_get_si (t2, MPFR_RNDN);
#if PRINT_VALUE
      mpfr_printf (" /* %40.16Rf */ ", t3);
#endif
      printf ("  { /* %5d */ UINT64_C (0x%016lX), %4ld },\n", i, value, exp);
    }
  }
  printf ("};\n\n");
  mpfr_clear (p64);
  mpfr_clear (p32);
  mpfr_clear (t1);
  mpfr_clear (t2);
  mpfr_clear (t3);
  mpfr_clear (t4);
  mpfr_free_cache ();
  printf ("static const uint64_t ipowers64[] = {\n");
  v = 1;
  for (i = 0; i < 20; i++) {
    printf ("  UINT64_C (%llu),\n", v);
    v *= 10;
  }
  printf ("};\n\n");
  printf ("static const uint32_t ipowers32[] = {\n");
  v = 1;
  for (i = 0; i < 10; i++) {
    printf ("  %llu,\n", v);
    v *= 10.0;
  }
  printf ("};\n");
  return 0;
}
#endif

static const struct
{
  uint64_t mul1;
  uint32_t mul2;
  int32_t exp;
} dpowers10[] = {
  { /*  -362 */ UINT64_C (0xB05135F614BB847D), 0x8DA339E8, -1266},
  { /*  -361 */ UINT64_C (0xDC65837399EA659C), 0xF10C0861, -1263},
  { /*  -360 */ UINT64_C (0x89BF722840327F82), 0x16A7853D, -1259},
  { /*  -359 */ UINT64_C (0xAC2F4EB2503F1F62), 0x9C51668C, -1256},
  { /*  -358 */ UINT64_C (0xD73B225EE44EE73B), 0x4365C02F, -1253},
  { /*  -357 */ UINT64_C (0x8684F57B4EB15085), 0x0A1F981D, -1249},
  { /*  -356 */ UINT64_C (0xA82632DA225DA4A6), 0x4CA77E25, -1246},
  { /*  -355 */ UINT64_C (0xD22FBF90AAF50DCF), 0xDFD15DAE, -1243},
  { /*  -354 */ UINT64_C (0x835DD7BA6AD928A1), 0xEBE2DA8D, -1239},
  { /*  -353 */ UINT64_C (0xA4354DA9058F72CA), 0x66DB9130, -1236},
  { /*  -352 */ UINT64_C (0xCD42A11346F34F7D), 0x0092757C, -1233},
  { /*  -351 */ UINT64_C (0x8049A4AC0C5811AE), 0x205B896D, -1229},
  { /*  -350 */ UINT64_C (0xA05C0DD70F6E1619), 0xA8726BC9, -1226},
  { /*  -349 */ UINT64_C (0xC873114CD3499BA0), 0x128F06BB, -1223},
  { /*  -348 */ UINT64_C (0xFA8FD5A0081C0288), 0x1732C86A, -1220},
  { /*  -347 */ UINT64_C (0x9C99E58405118195), 0x0E7FBD42, -1216},
  { /*  -346 */ UINT64_C (0xC3C05EE50655E1FA), 0x521FAC93, -1213},
  { /*  -345 */ UINT64_C (0xF4B0769E47EB5A78), 0xE6A797B7, -1210},
  { /*  -344 */ UINT64_C (0x98EE4A22ECF3188B), 0x9028BED3, -1206},
  { /*  -343 */ UINT64_C (0xBF29DCABA82FDEAE), 0x7432EE87, -1203},
  { /*  -342 */ UINT64_C (0xEEF453D6923BD65A), 0x113FAA29, -1200},
  { /*  -341 */ UINT64_C (0x9558B4661B6565F8), 0x4AC7CA5A, -1196},
  { /*  -340 */ UINT64_C (0xBAAEE17FA23EBF76), 0x5D79BCF0, -1193},
  { /*  -339 */ UINT64_C (0xE95A99DF8ACE6F53), 0xF4D82C2C, -1190},
  { /*  -338 */ UINT64_C (0x91D8A02BB6C10594), 0x79071B9C, -1186},
  { /*  -337 */ UINT64_C (0xB64EC836A47146F9), 0x9748E282, -1183},
  { /*  -336 */ UINT64_C (0xE3E27A444D8D98B7), 0xFD1B1B23, -1180},
  { /*  -335 */ UINT64_C (0x8E6D8C6AB0787F72), 0xFE30F0F6, -1176},
  { /*  -334 */ UINT64_C (0xB208EF855C969F4F), 0xBDBD2D33, -1173},
  { /*  -333 */ UINT64_C (0xDE8B2B66B3BC4723), 0xAD2C7880, -1170},
  { /*  -332 */ UINT64_C (0x8B16FB203055AC76), 0x4C3BCB50, -1166},
  { /*  -331 */ UINT64_C (0xADDCB9E83C6B1793), 0xDF4ABE24, -1163},
  { /*  -330 */ UINT64_C (0xD953E8624B85DD78), 0xD71D6DAD, -1160},
  { /*  -329 */ UINT64_C (0x87D4713D6F33AA6B), 0x8672648C, -1156},
  { /*  -328 */ UINT64_C (0xA9C98D8CCB009506), 0x680EFDAF, -1153},
  { /*  -327 */ UINT64_C (0xD43BF0EFFDC0BA48), 0x0212BD1B, -1150},
  { /*  -326 */ UINT64_C (0x84A57695FE98746D), 0x014BB631, -1146},
  { /*  -325 */ UINT64_C (0xA5CED43B7E3E9188), 0x419EA3BD, -1143},
  { /*  -324 */ UINT64_C (0xCF42894A5DCE35EA), 0x52064CAD, -1140},
  { /*  -323 */ UINT64_C (0x818995CE7AA0E1B2), 0x7343EFEC, -1136},
  { /*  -322 */ UINT64_C (0xA1EBFB4219491A1F), 0x1014EBE7, -1133},
  { /*  -321 */ UINT64_C (0xCA66FA129F9B60A6), 0xD41A26E0, -1130},
  { /*  -320 */ UINT64_C (0xFD00B897478238D0), 0x8920B099, -1127},
  { /*  -319 */ UINT64_C (0x9E20735E8CB16382), 0x55B46E5F, -1123},
  { /*  -318 */ UINT64_C (0xC5A890362FDDBC62), 0xEB2189F7, -1120},
  { /*  -317 */ UINT64_C (0xF712B443BBD52B7B), 0xA5E9EC75, -1117},
  { /*  -316 */ UINT64_C (0x9A6BB0AA55653B2D), 0x47B233C9, -1113},
  { /*  -315 */ UINT64_C (0xC1069CD4EABE89F8), 0x999EC0BB, -1110},
  { /*  -314 */ UINT64_C (0xF148440A256E2C76), 0xC00670EA, -1107},
  { /*  -313 */ UINT64_C (0x96CD2A865764DBCA), 0x38040692, -1103},
  { /*  -312 */ UINT64_C (0xBC807527ED3E12BC), 0xC6050837, -1100},
  { /*  -311 */ UINT64_C (0xEBA09271E88D976B), 0xF7864A45, -1097},
  { /*  -310 */ UINT64_C (0x93445B8731587EA3), 0x7AB3EE6B, -1093},
  { /*  -309 */ UINT64_C (0xB8157268FDAE9E4C), 0x5960EA06, -1090},
  { /*  -308 */ UINT64_C (0xE61ACF033D1A45DF), 0x6FB92487, -1087},
  { /*  -307 */ UINT64_C (0x8FD0C16206306BAB), 0xA5D3B6D4, -1083},
  { /*  -306 */ UINT64_C (0xB3C4F1BA87BC8696), 0x8F48A48A, -1080},
  { /*  -305 */ UINT64_C (0xE0B62E2929ABA83C), 0x331ACDAC, -1077},
  { /*  -304 */ UINT64_C (0x8C71DCD9BA0B4925), 0x9FF0C08B, -1073},
  { /*  -303 */ UINT64_C (0xAF8E5410288E1B6F), 0x07ECF0AE, -1070},
  { /*  -302 */ UINT64_C (0xDB71E91432B1A24A), 0xC9E82CDA, -1067},
  { /*  -301 */ UINT64_C (0x892731AC9FAF056E), 0xBE311C08, -1063},
  { /*  -300 */ UINT64_C (0xAB70FE17C79AC6CA), 0x6DBD630A, -1060},
  { /*  -299 */ UINT64_C (0xD64D3D9DB981787D), 0x092CBBCD, -1057},
  { /*  -298 */ UINT64_C (0x85F0468293F0EB4E), 0x25BBF560, -1053},
  { /*  -297 */ UINT64_C (0xA76C582338ED2621), 0xAF2AF2B8, -1050},
  { /*  -296 */ UINT64_C (0xD1476E2C07286FAA), 0x1AF5AF66, -1047},
  { /*  -295 */ UINT64_C (0x82CCA4DB847945CA), 0x50D98DA0, -1043},
  { /*  -294 */ UINT64_C (0xA37FCE126597973C), 0xE50FF108, -1040},
  { /*  -293 */ UINT64_C (0xCC5FC196FEFD7D0C), 0x1E53ED4A, -1037},
  { /*  -292 */ UINT64_C (0xFF77B1FCBEBCDC4F), 0x25E8E89C, -1034},
  { /*  -291 */ UINT64_C (0x9FAACF3DF73609B1), 0x77B19162, -1030},
  { /*  -290 */ UINT64_C (0xC795830D75038C1D), 0xD59DF5BA, -1027},
  { /*  -289 */ UINT64_C (0xF97AE3D0D2446F25), 0x4B057328, -1024},
  { /*  -288 */ UINT64_C (0x9BECCE62836AC577), 0x4EE367F9, -1020},
  { /*  -287 */ UINT64_C (0xC2E801FB244576D5), 0x229C41F8, -1017},
  { /*  -286 */ UINT64_C (0xF3A20279ED56D48A), 0x6B435275, -1014},
  { /*  -285 */ UINT64_C (0x9845418C345644D6), 0x830A1389, -1010},
  { /*  -284 */ UINT64_C (0xBE5691EF416BD60C), 0x23CC986C, -1007},
  { /*  -283 */ UINT64_C (0xEDEC366B11C6CB8F), 0x2CBFBE87, -1004},
  { /*  -282 */ UINT64_C (0x94B3A202EB1C3F39), 0x7BF7D714, -1000},
  { /*  -281 */ UINT64_C (0xB9E08A83A5E34F07), 0xDAF5CCD9, -997},
  { /*  -280 */ UINT64_C (0xE858AD248F5C22C9), 0xD1B34010, -994},
  { /*  -279 */ UINT64_C (0x91376C36D99995BE), 0x2310080A, -990},
  { /*  -278 */ UINT64_C (0xB58547448FFFFB2D), 0xABD40A0C, -987},
  { /*  -277 */ UINT64_C (0xE2E69915B3FFF9F9), 0x16C90C8F, -984},
  { /*  -276 */ UINT64_C (0x8DD01FAD907FFC3B), 0xAE3DA7D9, -980},
  { /*  -275 */ UINT64_C (0xB1442798F49FFB4A), 0x99CD11D0, -977},
  { /*  -274 */ UINT64_C (0xDD95317F31C7FA1D), 0x40405644, -974},
  { /*  -273 */ UINT64_C (0x8A7D3EEF7F1CFC52), 0x482835EA, -970},
  { /*  -272 */ UINT64_C (0xAD1C8EAB5EE43B66), 0xDA324365, -967},
  { /*  -271 */ UINT64_C (0xD863B256369D4A40), 0x90BED43E, -964},
  { /*  -270 */ UINT64_C (0x873E4F75E2224E68), 0x5A7744A7, -960},
  { /*  -269 */ UINT64_C (0xA90DE3535AAAE202), 0x711515D1, -957},
  { /*  -268 */ UINT64_C (0xD3515C2831559A83), 0x0D5A5B45, -954},
  { /*  -267 */ UINT64_C (0x8412D9991ED58091), 0xE858790B, -950},
  { /*  -266 */ UINT64_C (0xA5178FFF668AE0B6), 0x626E974E, -947},
  { /*  -265 */ UINT64_C (0xCE5D73FF402D98E3), 0xFB0A3D21, -944},
  { /*  -264 */ UINT64_C (0x80FA687F881C7F8E), 0x7CE66635, -940},
  { /*  -263 */ UINT64_C (0xA139029F6A239F72), 0x1C1FFFC2, -937},
  { /*  -262 */ UINT64_C (0xC987434744AC874E), 0xA327FFB2, -934},
  { /*  -261 */ UINT64_C (0xFBE9141915D7A922), 0x4BF1FF9F, -931},
  { /*  -260 */ UINT64_C (0x9D71AC8FADA6C9B5), 0x6F773FC3, -927},
  { /*  -259 */ UINT64_C (0xC4CE17B399107C22), 0xCB550FB4, -924},
  { /*  -258 */ UINT64_C (0xF6019DA07F549B2B), 0x7E2A53A1, -921},
  { /*  -257 */ UINT64_C (0x99C102844F94E0FB), 0x2EDA7445, -917},
  { /*  -256 */ UINT64_C (0xC0314325637A1939), 0xFA911156, -914},
  { /*  -255 */ UINT64_C (0xF03D93EEBC589F88), 0x793555AB, -911},
  { /*  -254 */ UINT64_C (0x96267C7535B763B5), 0x4BC1558B, -907},
  { /*  -253 */ UINT64_C (0xBBB01B9283253CA2), 0x9EB1AAEE, -904},
  { /*  -252 */ UINT64_C (0xEA9C227723EE8BCB), 0x465E15A9, -901},
  { /*  -251 */ UINT64_C (0x92A1958A7675175F), 0x0BFACD8A, -897},
  { /*  -250 */ UINT64_C (0xB749FAED14125D36), 0xCEF980EC, -894},
  { /*  -249 */ UINT64_C (0xE51C79A85916F484), 0x82B7E128, -891},
  { /*  -248 */ UINT64_C (0x8F31CC0937AE58D2), 0xD1B2ECB9, -887},
  { /*  -247 */ UINT64_C (0xB2FE3F0B8599EF07), 0x861FA7E7, -884},
  { /*  -246 */ UINT64_C (0xDFBDCECE67006AC9), 0x67A791E1, -881},
  { /*  -245 */ UINT64_C (0x8BD6A141006042BD), 0xE0C8BB2C, -877},
  { /*  -244 */ UINT64_C (0xAECC49914078536D), 0x58FAE9F7, -874},
  { /*  -243 */ UINT64_C (0xDA7F5BF590966848), 0xAF39A475, -871},
  { /*  -242 */ UINT64_C (0x888F99797A5E012D), 0x6D8406C9, -867},
  { /*  -241 */ UINT64_C (0xAAB37FD7D8F58178), 0xC8E5087C, -864},
  { /*  -240 */ UINT64_C (0xD5605FCDCF32E1D6), 0xFB1E4A9B, -861},
  { /*  -239 */ UINT64_C (0x855C3BE0A17FCD26), 0x5CF2EEA1, -857},
  { /*  -238 */ UINT64_C (0xA6B34AD8C9DFC06F), 0xF42FAA49, -854},
  { /*  -237 */ UINT64_C (0xD0601D8EFC57B08B), 0xF13B94DB, -851},
  { /*  -236 */ UINT64_C (0x823C12795DB6CE57), 0x76C53D09, -847},
  { /*  -235 */ UINT64_C (0xA2CB1717B52481ED), 0x54768C4B, -844},
  { /*  -234 */ UINT64_C (0xCB7DDCDDA26DA268), 0xA9942F5E, -841},
  { /*  -233 */ UINT64_C (0xFE5D54150B090B02), 0xD3F93B35, -838},
  { /*  -232 */ UINT64_C (0x9EFA548D26E5A6E1), 0xC47BC501, -834},
  { /*  -231 */ UINT64_C (0xC6B8E9B0709F109A), 0x359AB642, -831},
  { /*  -230 */ UINT64_C (0xF867241C8CC6D4C0), 0xC30163D2, -828},
  { /*  -229 */ UINT64_C (0x9B407691D7FC44F8), 0x79E0DE63, -824},
  { /*  -228 */ UINT64_C (0xC21094364DFB5636), 0x985915FC, -821},
  { /*  -227 */ UINT64_C (0xF294B943E17A2BC4), 0x3E6F5B7B, -818},
  { /*  -226 */ UINT64_C (0x979CF3CA6CEC5B5A), 0xA705992D, -814},
  { /*  -225 */ UINT64_C (0xBD8430BD08277231), 0x50C6FF78, -811},
  { /*  -224 */ UINT64_C (0xECE53CEC4A314EBD), 0xA4F8BF56, -808},
  { /*  -223 */ UINT64_C (0x940F4613AE5ED136), 0x871B7796, -804},
  { /*  -222 */ UINT64_C (0xB913179899F68584), 0x28E2557B, -801},
  { /*  -221 */ UINT64_C (0xE757DD7EC07426E5), 0x331AEADA, -798},
  { /*  -220 */ UINT64_C (0x9096EA6F3848984F), 0x3FF0D2C8, -794},
  { /*  -219 */ UINT64_C (0xB4BCA50B065ABE63), 0x0FED077A, -791},
  { /*  -218 */ UINT64_C (0xE1EBCE4DC7F16DFB), 0xD3E84959, -788},
  { /*  -217 */ UINT64_C (0x8D3360F09CF6E4BD), 0x64712DD8, -784},
  { /*  -216 */ UINT64_C (0xB080392CC4349DEC), 0xBD8D794E, -781},
  { /*  -215 */ UINT64_C (0xDCA04777F541C567), 0xECF0D7A1, -778},
  { /*  -214 */ UINT64_C (0x89E42CAAF9491B60), 0xF41686C5, -774},
  { /*  -213 */ UINT64_C (0xAC5D37D5B79B6239), 0x311C2876, -771},
  { /*  -212 */ UINT64_C (0xD77485CB25823AC7), 0x7D633293, -768},
  { /*  -211 */ UINT64_C (0x86A8D39EF77164BC), 0xAE5DFF9C, -764},
  { /*  -210 */ UINT64_C (0xA8530886B54DBDEB), 0xD9F57F83, -761},
  { /*  -209 */ UINT64_C (0xD267CAA862A12D66), 0xD072DF64, -758},
  { /*  -208 */ UINT64_C (0x8380DEA93DA4BC60), 0x4247CB9E, -754},
  { /*  -207 */ UINT64_C (0xA46116538D0DEB78), 0x52D9BE86, -751},
  { /*  -206 */ UINT64_C (0xCD795BE870516656), 0x67902E27, -748},
  { /*  -205 */ UINT64_C (0x806BD9714632DFF6), 0x00BA1CD9, -744},
  { /*  -204 */ UINT64_C (0xA086CFCD97BF97F3), 0x80E8A40F, -741},
  { /*  -203 */ UINT64_C (0xC8A883C0FDAF7DF0), 0x6122CD13, -738},
  { /*  -202 */ UINT64_C (0xFAD2A4B13D1B5D6C), 0x796B8057, -735},
  { /*  -201 */ UINT64_C (0x9CC3A6EEC6311A63), 0xCBE33036, -731},
  { /*  -200 */ UINT64_C (0xC3F490AA77BD60FC), 0xBEDBFC44, -728},
  { /*  -199 */ UINT64_C (0xF4F1B4D515ACB93B), 0xEE92FB55, -725},
  { /*  -198 */ UINT64_C (0x991711052D8BF3C5), 0x751BDD15, -721},
  { /*  -197 */ UINT64_C (0xBF5CD54678EEF0B6), 0xD262D45A, -718},
  { /*  -196 */ UINT64_C (0xEF340A98172AACE4), 0x86FB8971, -715},
  { /*  -195 */ UINT64_C (0x9580869F0E7AAC0E), 0xD45D35E7, -711},
  { /*  -194 */ UINT64_C (0xBAE0A846D2195712), 0x89748360, -708},
  { /*  -193 */ UINT64_C (0xE998D258869FACD7), 0x2BD1A438, -705},
  { /*  -192 */ UINT64_C (0x91FF83775423CC06), 0x7B6306A3, -701},
  { /*  -191 */ UINT64_C (0xB67F6455292CBF08), 0x1A3BC84C, -698},
  { /*  -190 */ UINT64_C (0xE41F3D6A7377EECA), 0x20CABA5F, -695},
  { /*  -189 */ UINT64_C (0x8E938662882AF53E), 0x547EB47B, -691},
  { /*  -188 */ UINT64_C (0xB23867FB2A35B28D), 0xE99E619A, -688},
  { /*  -187 */ UINT64_C (0xDEC681F9F4C31F31), 0x6405FA01, -685},
  { /*  -186 */ UINT64_C (0x8B3C113C38F9F37E), 0xDE83BC41, -681},
  { /*  -185 */ UINT64_C (0xAE0B158B4738705E), 0x9624AB51, -678},
  { /*  -184 */ UINT64_C (0xD98DDAEE19068C76), 0x3BADD625, -675},
  { /*  -183 */ UINT64_C (0x87F8A8D4CFA417C9), 0xE54CA5D7, -671},
  { /*  -182 */ UINT64_C (0xA9F6D30A038D1DBC), 0x5E9FCF4D, -668},
  { /*  -181 */ UINT64_C (0xD47487CC8470652B), 0x7647C320, -665},
  { /*  -180 */ UINT64_C (0x84C8D4DFD2C63F3B), 0x29ECD9F4, -661},
  { /*  -179 */ UINT64_C (0xA5FB0A17C777CF09), 0xF4681071, -658},
  { /*  -178 */ UINT64_C (0xCF79CC9DB955C2CC), 0x7182148D, -655},
  { /*  -177 */ UINT64_C (0x81AC1FE293D599BF), 0xC6F14CD8, -651},
  { /*  -176 */ UINT64_C (0xA21727DB38CB002F), 0xB8ADA00E, -648},
  { /*  -175 */ UINT64_C (0xCA9CF1D206FDC03B), 0xA6D90812, -645},
  { /*  -174 */ UINT64_C (0xFD442E4688BD304A), 0x908F4A16, -642},
  { /*  -173 */ UINT64_C (0x9E4A9CEC15763E2E), 0x9A598E4E, -638},
  { /*  -172 */ UINT64_C (0xC5DD44271AD3CDBA), 0x40EFF1E2, -635},
  { /*  -171 */ UINT64_C (0xF7549530E188C128), 0xD12BEE5A, -632},
  { /*  -170 */ UINT64_C (0x9A94DD3E8CF578B9), 0x82BB74F8, -628},
  { /*  -169 */ UINT64_C (0xC13A148E3032D6E7), 0xE36A5236, -625},
  { /*  -168 */ UINT64_C (0xF18899B1BC3F8CA1), 0xDC44E6C4, -622},
  { /*  -167 */ UINT64_C (0x96F5600F15A7B7E5), 0x29AB103A, -618},
  { /*  -166 */ UINT64_C (0xBCB2B812DB11A5DE), 0x7415D449, -615},
  { /*  -165 */ UINT64_C (0xEBDF661791D60F56), 0x111B495B, -612},
  { /*  -164 */ UINT64_C (0x936B9FCEBB25C995), 0xCAB10DD9, -608},
  { /*  -163 */ UINT64_C (0xB84687C269EF3BFB), 0x3D5D514F, -605},
  { /*  -162 */ UINT64_C (0xE65829B3046B0AFA), 0x0CB4A5A3, -602},
  { /*  -161 */ UINT64_C (0x8FF71A0FE2C2E6DC), 0x47F0E786, -598},
  { /*  -160 */ UINT64_C (0xB3F4E093DB73A093), 0x59ED2167, -595},
  { /*  -159 */ UINT64_C (0xE0F218B8D25088B8), 0x306869C1, -592},
  { /*  -158 */ UINT64_C (0x8C974F7383725573), 0x1E414219, -588},
  { /*  -157 */ UINT64_C (0xAFBD2350644EEACF), 0xE5D1929F, -585},
  { /*  -156 */ UINT64_C (0xDBAC6C247D62A583), 0xDF45F747, -582},
  { /*  -155 */ UINT64_C (0x894BC396CE5DA772), 0x6B8BBA8C, -578},
  { /*  -154 */ UINT64_C (0xAB9EB47C81F5114F), 0x066EA92F, -575},
  { /*  -153 */ UINT64_C (0xD686619BA27255A2), 0xC80A537B, -572},
  { /*  -152 */ UINT64_C (0x8613FD0145877585), 0xBD06742D, -568},
  { /*  -151 */ UINT64_C (0xA798FC4196E952E7), 0x2C481138, -565},
  { /*  -150 */ UINT64_C (0xD17F3B51FCA3A7A0), 0xF75A1586, -562},
  { /*  -149 */ UINT64_C (0x82EF85133DE648C4), 0x9A984D74, -558},
  { /*  -148 */ UINT64_C (0xA3AB66580D5FDAF5), 0xC13E60D1, -555},
  { /*  -147 */ UINT64_C (0xCC963FEE10B7D1B3), 0x318DF905, -552},
  { /*  -146 */ UINT64_C (0xFFBBCFE994E5C61F), 0xFDF17746, -549},
  { /*  -145 */ UINT64_C (0x9FD561F1FD0F9BD3), 0xFEB6EA8C, -545},
  { /*  -144 */ UINT64_C (0xC7CABA6E7C5382C8), 0xFE64A52F, -542},
  { /*  -143 */ UINT64_C (0xF9BD690A1B68637B), 0x3DFDCE7B, -539},
  { /*  -142 */ UINT64_C (0x9C1661A651213E2D), 0x06BEA10D, -535},
  { /*  -141 */ UINT64_C (0xC31BFA0FE5698DB8), 0x486E4950, -532},
  { /*  -140 */ UINT64_C (0xF3E2F893DEC3F126), 0x5A89DBA4, -529},
  { /*  -139 */ UINT64_C (0x986DDB5C6B3A76B7), 0xF8962946, -525},
  { /*  -138 */ UINT64_C (0xBE89523386091465), 0xF6BBB398, -522},
  { /*  -137 */ UINT64_C (0xEE2BA6C0678B597F), 0x746AA07E, -519},
  { /*  -136 */ UINT64_C (0x94DB483840B717EF), 0xA8C2A44F, -515},
  { /*  -135 */ UINT64_C (0xBA121A4650E4DDEB), 0x92F34D62, -512},
  { /*  -134 */ UINT64_C (0xE896A0D7E51E1566), 0x77B020BB, -509},
  { /*  -133 */ UINT64_C (0x915E2486EF32CD60), 0x0ACE1475, -505},
  { /*  -132 */ UINT64_C (0xB5B5ADA8AAFF80B8), 0x0D819992, -502},
  { /*  -131 */ UINT64_C (0xE3231912D5BF60E6), 0x10E1FFF7, -499},
  { /*  -130 */ UINT64_C (0x8DF5EFABC5979C8F), 0xCA8D3FFA, -495},
  { /*  -129 */ UINT64_C (0xB1736B96B6FD83B3), 0xBD308FF9, -492},
  { /*  -128 */ UINT64_C (0xDDD0467C64BCE4A0), 0xAC7CB3F7, -489},
  { /*  -127 */ UINT64_C (0x8AA22C0DBEF60EE4), 0x6BCDF07A, -485},
  { /*  -126 */ UINT64_C (0xAD4AB7112EB3929D), 0x86C16C99, -482},
  { /*  -125 */ UINT64_C (0xD89D64D57A607744), 0xE871C7BF, -479},
  { /*  -124 */ UINT64_C (0x87625F056C7C4A8B), 0x11471CD7, -475},
  { /*  -123 */ UINT64_C (0xA93AF6C6C79B5D2D), 0xD598E40D, -472},
  { /*  -122 */ UINT64_C (0xD389B47879823479), 0x4AFF1D11, -469},
  { /*  -121 */ UINT64_C (0x843610CB4BF160CB), 0xCEDF722A, -465},
  { /*  -120 */ UINT64_C (0xA54394FE1EEDB8FE), 0xC2974EB5, -462},
  { /*  -119 */ UINT64_C (0xCE947A3DA6A9273E), 0x733D2262, -459},
  { /*  -118 */ UINT64_C (0x811CCC668829B887), 0x0806357D, -455},
  { /*  -117 */ UINT64_C (0xA163FF802A3426A8), 0xCA07C2DD, -452},
  { /*  -116 */ UINT64_C (0xC9BCFF6034C13052), 0xFC89B394, -449},
  { /*  -115 */ UINT64_C (0xFC2C3F3841F17C67), 0xBBAC2079, -446},
  { /*  -114 */ UINT64_C (0x9D9BA7832936EDC0), 0xD54B944C, -442},
  { /*  -113 */ UINT64_C (0xC5029163F384A931), 0x0A9E795E, -439},
  { /*  -112 */ UINT64_C (0xF64335BCF065D37D), 0x4D4617B6, -436},
  { /*  -111 */ UINT64_C (0x99EA0196163FA42E), 0x504BCED2, -432},
  { /*  -110 */ UINT64_C (0xC06481FB9BCF8D39), 0xE45EC286, -429},
  { /*  -109 */ UINT64_C (0xF07DA27A82C37088), 0x5D767328, -426},
  { /*  -108 */ UINT64_C (0x964E858C91BA2655), 0x3A6A07F9, -422},
  { /*  -107 */ UINT64_C (0xBBE226EFB628AFEA), 0x890489F7, -419},
  { /*  -106 */ UINT64_C (0xEADAB0ABA3B2DBE5), 0x2B45AC75, -416},
  { /*  -105 */ UINT64_C (0x92C8AE6B464FC96F), 0x3B0B8BC9, -412},
  { /*  -104 */ UINT64_C (0xB77ADA0617E3BBCB), 0x09CE6EBB, -409},
  { /*  -103 */ UINT64_C (0xE55990879DDCAABD), 0xCC420A6A, -406},
  { /*  -102 */ UINT64_C (0x8F57FA54C2A9EAB6), 0x9FA94682, -402},
  { /*  -101 */ UINT64_C (0xB32DF8E9F3546564), 0x47939823, -399},
  { /*  -100 */ UINT64_C (0xDFF9772470297EBD), 0x59787E2C, -396},
  { /*   -99 */ UINT64_C (0x8BFBEA76C619EF36), 0x57EB4EDB, -392},
  { /*   -98 */ UINT64_C (0xAEFAE51477A06B03), 0xEDE62292, -389},
  { /*   -97 */ UINT64_C (0xDAB99E59958885C4), 0xE95FAB37, -386},
  { /*   -96 */ UINT64_C (0x88B402F7FD75539B), 0x11DBCB02, -382},
  { /*   -95 */ UINT64_C (0xAAE103B5FCD2A881), 0xD652BDC3, -379},
  { /*   -94 */ UINT64_C (0xD59944A37C0752A2), 0x4BE76D33, -376},
  { /*   -93 */ UINT64_C (0x857FCAE62D8493A5), 0x6F70A440, -372},
  { /*   -92 */ UINT64_C (0xA6DFBD9FB8E5B88E), 0xCB4CCD50, -369},
  { /*   -91 */ UINT64_C (0xD097AD07A71F26B2), 0x7E2000A4, -366},
  { /*   -90 */ UINT64_C (0x825ECC24C873782F), 0x8ED40067, -362},
  { /*   -89 */ UINT64_C (0xA2F67F2DFA90563B), 0x72890080, -359},
  { /*   -88 */ UINT64_C (0xCBB41EF979346BCA), 0x4F2B40A0, -356},
  { /*   -87 */ UINT64_C (0xFEA126B7D78186BC), 0xE2F610C8, -353},
  { /*   -86 */ UINT64_C (0x9F24B832E6B0F436), 0x0DD9CA7D, -349},
  { /*   -85 */ UINT64_C (0xC6EDE63FA05D3143), 0x91503D1C, -346},
  { /*   -84 */ UINT64_C (0xF8A95FCF88747D94), 0x75A44C64, -343},
  { /*   -83 */ UINT64_C (0x9B69DBE1B548CE7C), 0xC986AFBE, -339},
  { /*   -82 */ UINT64_C (0xC24452DA229B021B), 0xFBE85BAE, -336},
  { /*   -81 */ UINT64_C (0xF2D56790AB41C2A2), 0xFAE27299, -333},
  { /*   -80 */ UINT64_C (0x97C560BA6B0919A5), 0xDCCD87A0, -329},
  { /*   -79 */ UINT64_C (0xBDB6B8E905CB600F), 0x5400E988, -326},
  { /*   -78 */ UINT64_C (0xED246723473E3813), 0x290123EA, -323},
  { /*   -77 */ UINT64_C (0x9436C0760C86E30B), 0xF9A0B672, -319},
  { /*   -76 */ UINT64_C (0xB94470938FA89BCE), 0xF808E40F, -316},
  { /*   -75 */ UINT64_C (0xE7958CB87392C2C2), 0xB60B1D12, -313},
  { /*   -74 */ UINT64_C (0x90BD77F3483BB9B9), 0xB1C6F22B, -309},
  { /*   -73 */ UINT64_C (0xB4ECD5F01A4AA828), 0x1E38AEB6, -306},
  { /*   -72 */ UINT64_C (0xE2280B6C20DD5232), 0x25C6DA64, -303},
  { /*   -71 */ UINT64_C (0x8D590723948A535F), 0x579C487E, -299},
  { /*   -70 */ UINT64_C (0xB0AF48EC79ACE837), 0x2D835A9E, -296},
  { /*   -69 */ UINT64_C (0xDCDB1B2798182244), 0xF8E43145, -293},
  { /*   -68 */ UINT64_C (0x8A08F0F8BF0F156B), 0x1B8E9ECB, -289},
  { /*   -67 */ UINT64_C (0xAC8B2D36EED2DAC5), 0xE272467E, -286},
  { /*   -66 */ UINT64_C (0xD7ADF884AA879177), 0x5B0ED81E, -283},
  { /*   -65 */ UINT64_C (0x86CCBB52EA94BAEA), 0x98E94713, -279},
  { /*   -64 */ UINT64_C (0xA87FEA27A539E9A5), 0x3F2398D7, -276},
  { /*   -63 */ UINT64_C (0xD29FE4B18E88640E), 0x8EEC7F0D, -273},
  { /*   -62 */ UINT64_C (0x83A3EEEEF9153E89), 0x1953CF68, -269},
  { /*   -61 */ UINT64_C (0xA48CEAAAB75A8E2B), 0x5FA8C342, -266},
  { /*   -60 */ UINT64_C (0xCDB02555653131B6), 0x3792F413, -263},
  { /*   -59 */ UINT64_C (0x808E17555F3EBF11), 0xE2BBD88C, -259},
  { /*   -58 */ UINT64_C (0xA0B19D2AB70E6ED6), 0x5B6ACEAF, -256},
  { /*   -57 */ UINT64_C (0xC8DE047564D20A8B), 0xF245825A, -253},
  { /*   -56 */ UINT64_C (0xFB158592BE068D2E), 0xEED6E2F1, -250},
  { /*   -55 */ UINT64_C (0x9CED737BB6C4183D), 0x55464DD7, -246},
  { /*   -54 */ UINT64_C (0xC428D05AA4751E4C), 0xAA97E14C, -243},
  { /*   -53 */ UINT64_C (0xF53304714D9265DF), 0xD53DD99F, -240},
  { /*   -52 */ UINT64_C (0x993FE2C6D07B7FAB), 0xE546A804, -236},
  { /*   -51 */ UINT64_C (0xBF8FDB78849A5F96), 0xDE985204, -233},
  { /*   -50 */ UINT64_C (0xEF73D256A5C0F77C), 0x963E6686, -230},
  { /*   -49 */ UINT64_C (0x95A8637627989AAD), 0xDDE70013, -226},
  { /*   -48 */ UINT64_C (0xBB127C53B17EC159), 0x5560C018, -223},
  { /*   -47 */ UINT64_C (0xE9D71B689DDE71AF), 0xAAB8F01E, -220},
  { /*   -46 */ UINT64_C (0x9226712162AB070D), 0xCAB39613, -216},
  { /*   -45 */ UINT64_C (0xB6B00D69BB55C8D1), 0x3D607B98, -213},
  { /*   -44 */ UINT64_C (0xE45C10C42A2B3B05), 0x8CB89A7E, -210},
  { /*   -43 */ UINT64_C (0x8EB98A7A9A5B04E3), 0x77F3608F, -206},
  { /*   -42 */ UINT64_C (0xB267ED1940F1C61C), 0x55F038B2, -203},
  { /*   -41 */ UINT64_C (0xDF01E85F912E37A3), 0x6B6C46DF, -200},
  { /*   -40 */ UINT64_C (0x8B61313BBABCE2C6), 0x2323AC4B, -196},
  { /*   -39 */ UINT64_C (0xAE397D8AA96C1B77), 0xABEC975E, -193},
  { /*   -38 */ UINT64_C (0xD9C7DCED53C72255), 0x96E7BD36, -190},
  { /*   -37 */ UINT64_C (0x881CEA14545C7575), 0x7E50D641, -186},
  { /*   -36 */ UINT64_C (0xAA242499697392D2), 0xDDE50BD2, -183},
  { /*   -35 */ UINT64_C (0xD4AD2DBFC3D07787), 0x955E4EC6, -180},
  { /*   -34 */ UINT64_C (0x84EC3C97DA624AB4), 0xBD5AF13C, -176},
  { /*   -33 */ UINT64_C (0xA6274BBDD0FADD61), 0xECB1AD8B, -173},
  { /*   -32 */ UINT64_C (0xCFB11EAD453994BA), 0x67DE18EE, -170},
  { /*   -31 */ UINT64_C (0x81CEB32C4B43FCF4), 0x80EACF95, -166},
  { /*   -30 */ UINT64_C (0xA2425FF75E14FC31), 0xA125837A, -163},
  { /*   -29 */ UINT64_C (0xCAD2F7F5359A3B3E), 0x096EE458, -160},
  { /*   -28 */ UINT64_C (0xFD87B5F28300CA0D), 0x8BCA9D6E, -157},
  { /*   -27 */ UINT64_C (0x9E74D1B791E07E48), 0x775EA265, -153},
  { /*   -26 */ UINT64_C (0xC612062576589DDA), 0x95364AFE, -150},
  { /*   -25 */ UINT64_C (0xF79687AED3EEC551), 0x3A83DDBE, -147},
  { /*   -24 */ UINT64_C (0x9ABE14CD44753B52), 0xC4926A96, -143},
  { /*   -23 */ UINT64_C (0xC16D9A0095928A27), 0x75B7053C, -140},
  { /*   -22 */ UINT64_C (0xF1C90080BAF72CB1), 0x5324C68B, -137},
  { /*   -21 */ UINT64_C (0x971DA05074DA7BEE), 0xD3F6FC17, -133},
  { /*   -20 */ UINT64_C (0xBCE5086492111AEA), 0x88F4BB1D, -130},
  { /*   -19 */ UINT64_C (0xEC1E4A7DB69561A5), 0x2B31E9E4, -127},
  { /*   -18 */ UINT64_C (0x9392EE8E921D5D07), 0x3AFF322E, -123},
  { /*   -17 */ UINT64_C (0xB877AA3236A4B449), 0x09BEFEBA, -120},
  { /*   -16 */ UINT64_C (0xE69594BEC44DE15B), 0x4C2EBE68, -117},
  { /*   -15 */ UINT64_C (0x901D7CF73AB0ACD9), 0x0F9D3701, -113},
  { /*   -14 */ UINT64_C (0xB424DC35095CD80F), 0x538484C2, -110},
  { /*   -13 */ UINT64_C (0xE12E13424BB40E13), 0x2865A5F2, -107},
  { /*   -12 */ UINT64_C (0x8CBCCC096F5088CB), 0xF93F87B7, -103},
  { /*   -11 */ UINT64_C (0xAFEBFF0BCB24AAFE), 0xF78F69A5, -100},
  { /*   -10 */ UINT64_C (0xDBE6FECEBDEDD5BE), 0xB573440E, -97},
  { /*    -9 */ UINT64_C (0x89705F4136B4A597), 0x31680A89, -93},
  { /*    -8 */ UINT64_C (0xABCC77118461CEFC), 0xFDC20D2B, -90},
  { /*    -7 */ UINT64_C (0xD6BF94D5E57A42BC), 0x3D329076, -87},
  { /*    -6 */ UINT64_C (0x8637BD05AF6C69B5), 0xA63F9A4A, -83},
  { /*    -5 */ UINT64_C (0xA7C5AC471B478423), 0x0FCF80DC, -80},
  { /*    -4 */ UINT64_C (0xD1B71758E219652B), 0xD3C36113, -77},
  { /*    -3 */ UINT64_C (0x83126E978D4FDF3B), 0x645A1CAC, -73},
  { /*    -2 */ UINT64_C (0xA3D70A3D70A3D70A), 0x3D70A3D7, -70},
  { /*    -1 */ UINT64_C (0xCCCCCCCCCCCCCCCC), 0xCCCCCCCD, -67},
  { /*     0 */ UINT64_C (0x8000000000000000), 0x00000000, -63},
  { /*     1 */ UINT64_C (0xA000000000000000), 0x00000000, -60},
  { /*     2 */ UINT64_C (0xC800000000000000), 0x00000000, -57},
  { /*     3 */ UINT64_C (0xFA00000000000000), 0x00000000, -54},
  { /*     4 */ UINT64_C (0x9C40000000000000), 0x00000000, -50},
  { /*     5 */ UINT64_C (0xC350000000000000), 0x00000000, -47},
  { /*     6 */ UINT64_C (0xF424000000000000), 0x00000000, -44},
  { /*     7 */ UINT64_C (0x9896800000000000), 0x00000000, -40},
  { /*     8 */ UINT64_C (0xBEBC200000000000), 0x00000000, -37},
  { /*     9 */ UINT64_C (0xEE6B280000000000), 0x00000000, -34},
  { /*    10 */ UINT64_C (0x9502F90000000000), 0x00000000, -30},
  { /*    11 */ UINT64_C (0xBA43B74000000000), 0x00000000, -27},
  { /*    12 */ UINT64_C (0xE8D4A51000000000), 0x00000000, -24},
  { /*    13 */ UINT64_C (0x9184E72A00000000), 0x00000000, -20},
  { /*    14 */ UINT64_C (0xB5E620F480000000), 0x00000000, -17},
  { /*    15 */ UINT64_C (0xE35FA931A0000000), 0x00000000, -14},
  { /*    16 */ UINT64_C (0x8E1BC9BF04000000), 0x00000000, -10},
  { /*    17 */ UINT64_C (0xB1A2BC2EC5000000), 0x00000000, -7},
  { /*    18 */ UINT64_C (0xDE0B6B3A76400000), 0x00000000, -4},
  { /*    19 */ UINT64_C (0x8AC7230489E80000), 0x00000000, 0},
  { /*    20 */ UINT64_C (0xAD78EBC5AC620000), 0x00000000, 3},
  { /*    21 */ UINT64_C (0xD8D726B7177A8000), 0x00000000, 6},
  { /*    22 */ UINT64_C (0x878678326EAC9000), 0x00000000, 10},
  { /*    23 */ UINT64_C (0xA968163F0A57B400), 0x00000000, 13},
  { /*    24 */ UINT64_C (0xD3C21BCECCEDA100), 0x00000000, 16},
  { /*    25 */ UINT64_C (0x84595161401484A0), 0x00000000, 20},
  { /*    26 */ UINT64_C (0xA56FA5B99019A5C8), 0x00000000, 23},
  { /*    27 */ UINT64_C (0xCECB8F27F4200F3A), 0x00000000, 26},
  { /*    28 */ UINT64_C (0x813F3978F8940984), 0x40000000, 30},
  { /*    29 */ UINT64_C (0xA18F07D736B90BE5), 0x50000000, 33},
  { /*    30 */ UINT64_C (0xC9F2C9CD04674EDE), 0xA4000000, 36},
  { /*    31 */ UINT64_C (0xFC6F7C4045812296), 0x4D000000, 39},
  { /*    32 */ UINT64_C (0x9DC5ADA82B70B59D), 0xF0200000, 43},
  { /*    33 */ UINT64_C (0xC5371912364CE305), 0x6C280000, 46},
  { /*    34 */ UINT64_C (0xF684DF56C3E01BC6), 0xC7320000, 49},
  { /*    35 */ UINT64_C (0x9A130B963A6C115C), 0x3C7F4000, 53},
  { /*    36 */ UINT64_C (0xC097CE7BC90715B3), 0x4B9F1000, 56},
  { /*    37 */ UINT64_C (0xF0BDC21ABB48DB20), 0x1E86D400, 59},
  { /*    38 */ UINT64_C (0x96769950B50D88F4), 0x13144480, 63},
  { /*    39 */ UINT64_C (0xBC143FA4E250EB31), 0x17D955A0, 66},
  { /*    40 */ UINT64_C (0xEB194F8E1AE525FD), 0x5DCFAB08, 69},
  { /*    41 */ UINT64_C (0x92EFD1B8D0CF37BE), 0x5AA1CAE5, 73},
  { /*    42 */ UINT64_C (0xB7ABC627050305AD), 0xF14A3D9E, 76},
  { /*    43 */ UINT64_C (0xE596B7B0C643C719), 0x6D9CCD06, 79},
  { /*    44 */ UINT64_C (0x8F7E32CE7BEA5C6F), 0xE4820024, 83},
  { /*    45 */ UINT64_C (0xB35DBF821AE4F38B), 0xDDA2802D, 86},
  { /*    46 */ UINT64_C (0xE0352F62A19E306E), 0xD50B2038, 89},
  { /*    47 */ UINT64_C (0x8C213D9DA502DE45), 0x4526F423, 93},
  { /*    48 */ UINT64_C (0xAF298D050E4395D6), 0x9670B12B, 96},
  { /*    49 */ UINT64_C (0xDAF3F04651D47B4C), 0x3C0CDD76, 99},
  { /*    50 */ UINT64_C (0x88D8762BF324CD0F), 0xA5880A6A, 103},
  { /*    51 */ UINT64_C (0xAB0E93B6EFEE0053), 0x8EEA0D04, 106},
  { /*    52 */ UINT64_C (0xD5D238A4ABE98068), 0x72A49046, 109},
  { /*    53 */ UINT64_C (0x85A36366EB71F041), 0x47A6DA2B, 113},
  { /*    54 */ UINT64_C (0xA70C3C40A64E6C51), 0x999090B6, 116},
  { /*    55 */ UINT64_C (0xD0CF4B50CFE20765), 0xFFF4B4E4, 119},
  { /*    56 */ UINT64_C (0x82818F1281ED449F), 0xBFF8F10E, 123},
  { /*    57 */ UINT64_C (0xA321F2D7226895C7), 0xAFF72D52, 126},
  { /*    58 */ UINT64_C (0xCBEA6F8CEB02BB39), 0x9BF4F8A7, 129},
  { /*    59 */ UINT64_C (0xFEE50B7025C36A08), 0x02F236D0, 132},
  { /*    60 */ UINT64_C (0x9F4F2726179A2245), 0x01D76242, 136},
  { /*    61 */ UINT64_C (0xC722F0EF9D80AAD6), 0x424D3AD3, 139},
  { /*    62 */ UINT64_C (0xF8EBAD2B84E0D58B), 0xD2E08987, 142},
  { /*    63 */ UINT64_C (0x9B934C3B330C8577), 0x63CC55F5, 146},
  { /*    64 */ UINT64_C (0xC2781F49FFCFA6D5), 0x3CBF6B72, 149},
  { /*    65 */ UINT64_C (0xF316271C7FC3908A), 0x8BEF464E, 152},
  { /*    66 */ UINT64_C (0x97EDD871CFDA3A56), 0x97758BF1, 156},
  { /*    67 */ UINT64_C (0xBDE94E8E43D0C8EC), 0x3D52EEED, 159},
  { /*    68 */ UINT64_C (0xED63A231D4C4FB27), 0x4CA7AAA8, 162},
  { /*    69 */ UINT64_C (0x945E455F24FB1CF8), 0x8FE8CAA9, 166},
  { /*    70 */ UINT64_C (0xB975D6B6EE39E436), 0xB3E2FD54, 169},
  { /*    71 */ UINT64_C (0xE7D34C64A9C85D44), 0x60DBBCA8, 172},
  { /*    72 */ UINT64_C (0x90E40FBEEA1D3A4A), 0xBC8955E9, 176},
  { /*    73 */ UINT64_C (0xB51D13AEA4A488DD), 0x6BABAB64, 179},
  { /*    74 */ UINT64_C (0xE264589A4DCDAB14), 0xC696963C, 182},
  { /*    75 */ UINT64_C (0x8D7EB76070A08AEC), 0xFC1E1DE6, 186},
  { /*    76 */ UINT64_C (0xB0DE65388CC8ADA8), 0x3B25A55F, 189},
  { /*    77 */ UINT64_C (0xDD15FE86AFFAD912), 0x49EF0EB7, 192},
  { /*    78 */ UINT64_C (0x8A2DBF142DFCC7AB), 0x6E356932, 196},
  { /*    79 */ UINT64_C (0xACB92ED9397BF996), 0x49C2C37F, 199},
  { /*    80 */ UINT64_C (0xD7E77A8F87DAF7FB), 0xDC33745F, 202},
  { /*    81 */ UINT64_C (0x86F0AC99B4E8DAFD), 0x69A028BB, 206},
  { /*    82 */ UINT64_C (0xA8ACD7C0222311BC), 0xC40832EA, 209},
  { /*    83 */ UINT64_C (0xD2D80DB02AABD62B), 0xF50A3FA5, 212},
  { /*    84 */ UINT64_C (0x83C7088E1AAB65DB), 0x792667C7, 216},
  { /*    85 */ UINT64_C (0xA4B8CAB1A1563F52), 0x577001B9, 219},
  { /*    86 */ UINT64_C (0xCDE6FD5E09ABCF26), 0xED4C0227, 222},
  { /*    87 */ UINT64_C (0x80B05E5AC60B6178), 0x544F8158, 226},
  { /*    88 */ UINT64_C (0xA0DC75F1778E39D6), 0x696361AE, 229},
  { /*    89 */ UINT64_C (0xC913936DD571C84C), 0x03BC3A1A, 232},
  { /*    90 */ UINT64_C (0xFB5878494ACE3A5F), 0x04AB48A0, 235},
  { /*    91 */ UINT64_C (0x9D174B2DCEC0E47B), 0x62EB0D64, 239},
  { /*    92 */ UINT64_C (0xC45D1DF942711D9A), 0x3BA5D0BD, 242},
  { /*    93 */ UINT64_C (0xF5746577930D6500), 0xCA8F44EC, 245},
  { /*    94 */ UINT64_C (0x9968BF6ABBE85F20), 0x7E998B14, 249},
  { /*    95 */ UINT64_C (0xBFC2EF456AE276E8), 0x9E3FEDD9, 252},
  { /*    96 */ UINT64_C (0xEFB3AB16C59B14A2), 0xC5CFE94F, 255},
  { /*    97 */ UINT64_C (0x95D04AEE3B80ECE5), 0xBBA1F1D1, 259},
  { /*    98 */ UINT64_C (0xBB445DA9CA61281F), 0x2A8A6E46, 262},
  { /*    99 */ UINT64_C (0xEA1575143CF97226), 0xF52D09D7, 265},
  { /*   100 */ UINT64_C (0x924D692CA61BE758), 0x593C2626, 269},
  { /*   101 */ UINT64_C (0xB6E0C377CFA2E12E), 0x6F8B2FB0, 272},
  { /*   102 */ UINT64_C (0xE498F455C38B997A), 0x0B6DFB9C, 275},
  { /*   103 */ UINT64_C (0x8EDF98B59A373FEC), 0x4724BD42, 279},
  { /*   104 */ UINT64_C (0xB2977EE300C50FE7), 0x58EDEC92, 282},
  { /*   105 */ UINT64_C (0xDF3D5E9BC0F653E1), 0x2F2967B6, 285},
  { /*   106 */ UINT64_C (0x8B865B215899F46C), 0xBD79E0D2, 289},
  { /*   107 */ UINT64_C (0xAE67F1E9AEC07187), 0xECD85907, 292},
  { /*   108 */ UINT64_C (0xDA01EE641A708DE9), 0xE80E6F48, 295},
  { /*   109 */ UINT64_C (0x884134FE908658B2), 0x3109058D, 299},
  { /*   110 */ UINT64_C (0xAA51823E34A7EEDE), 0xBD4B46F0, 302},
  { /*   111 */ UINT64_C (0xD4E5E2CDC1D1EA96), 0x6C9E18AC, 305},
  { /*   112 */ UINT64_C (0x850FADC09923329E), 0x03E2CF6C, 309},
  { /*   113 */ UINT64_C (0xA6539930BF6BFF45), 0x84DB8347, 312},
  { /*   114 */ UINT64_C (0xCFE87F7CEF46FF16), 0xE6126418, 315},
  { /*   115 */ UINT64_C (0x81F14FAE158C5F6E), 0x4FCB7E8F, 319},
  { /*   116 */ UINT64_C (0xA26DA3999AEF7749), 0xE3BE5E33, 322},
  { /*   117 */ UINT64_C (0xCB090C8001AB551C), 0x5CADF5C0, 325},
  { /*   118 */ UINT64_C (0xFDCB4FA002162A63), 0x73D97330, 328},
  { /*   119 */ UINT64_C (0x9E9F11C4014DDA7E), 0x2867E7FE, 332},
  { /*   120 */ UINT64_C (0xC646D63501A1511D), 0xB281E1FD, 335},
  { /*   121 */ UINT64_C (0xF7D88BC24209A565), 0x1F225A7D, 338},
  { /*   122 */ UINT64_C (0x9AE757596946075F), 0x3375788E, 342},
  { /*   123 */ UINT64_C (0xC1A12D2FC3978937), 0x0052D6B1, 345},
  { /*   124 */ UINT64_C (0xF209787BB47D6B84), 0xC0678C5E, 348},
  { /*   125 */ UINT64_C (0x9745EB4D50CE6332), 0xF840B7BB, 352},
  { /*   126 */ UINT64_C (0xBD176620A501FBFF), 0xB650E5A9, 355},
  { /*   127 */ UINT64_C (0xEC5D3FA8CE427AFF), 0xA3E51F14, 358},
  { /*   128 */ UINT64_C (0x93BA47C980E98CDF), 0xC66F336C, 362},
  { /*   129 */ UINT64_C (0xB8A8D9BBE123F017), 0xB80B0047, 365},
  { /*   130 */ UINT64_C (0xE6D3102AD96CEC1D), 0xA60DC059, 368},
  { /*   131 */ UINT64_C (0x9043EA1AC7E41392), 0x87C89838, 372},
  { /*   132 */ UINT64_C (0xB454E4A179DD1877), 0x29BABE46, 375},
  { /*   133 */ UINT64_C (0xE16A1DC9D8545E94), 0xF4296DD7, 378},
  { /*   134 */ UINT64_C (0x8CE2529E2734BB1D), 0x1899E4A6, 382},
  { /*   135 */ UINT64_C (0xB01AE745B101E9E4), 0x5EC05DD0, 385},
  { /*   136 */ UINT64_C (0xDC21A1171D42645D), 0x76707544, 388},
  { /*   137 */ UINT64_C (0x899504AE72497EBA), 0x6A06494A, 392},
  { /*   138 */ UINT64_C (0xABFA45DA0EDBDE69), 0x0487DB9D, 395},
  { /*   139 */ UINT64_C (0xD6F8D7509292D603), 0x45A9D284, 398},
  { /*   140 */ UINT64_C (0x865B86925B9BC5C2), 0x0B8A2393, 402},
  { /*   141 */ UINT64_C (0xA7F26836F282B732), 0x8E6CAC77, 405},
  { /*   142 */ UINT64_C (0xD1EF0244AF2364FF), 0x3207D795, 408},
  { /*   143 */ UINT64_C (0x8335616AED761F1F), 0x7F44E6BD, 412},
  { /*   144 */ UINT64_C (0xA402B9C5A8D3A6E7), 0x5F16206D, 415},
  { /*   145 */ UINT64_C (0xCD036837130890A1), 0x36DBA888, 418},
  { /*   146 */ UINT64_C (0x802221226BE55A64), 0xC2494955, 422},
  { /*   147 */ UINT64_C (0xA02AA96B06DEB0FD), 0xF2DB9BAA, 425},
  { /*   148 */ UINT64_C (0xC83553C5C8965D3D), 0x6F928295, 428},
  { /*   149 */ UINT64_C (0xFA42A8B73ABBF48C), 0xCB77233A, 431},
  { /*   150 */ UINT64_C (0x9C69A97284B578D7), 0xFF2A7604, 435},
  { /*   151 */ UINT64_C (0xC38413CF25E2D70D), 0xFEF51385, 438},
  { /*   152 */ UINT64_C (0xF46518C2EF5B8CD1), 0x7EB25866, 441},
  { /*   153 */ UINT64_C (0x98BF2F79D5993802), 0xEF2F7740, 445},
  { /*   154 */ UINT64_C (0xBEEEFB584AFF8603), 0xAAFB5510, 448},
  { /*   155 */ UINT64_C (0xEEAABA2E5DBF6784), 0x95BA2A54, 451},
  { /*   156 */ UINT64_C (0x952AB45CFA97A0B2), 0xDD945A74, 455},
  { /*   157 */ UINT64_C (0xBA756174393D88DF), 0x94F97112, 458},
  { /*   158 */ UINT64_C (0xE912B9D1478CEB17), 0x7A37CD56, 461},
  { /*   159 */ UINT64_C (0x91ABB422CCB812EE), 0xAC62E056, 465},
  { /*   160 */ UINT64_C (0xB616A12B7FE617AA), 0x577B986B, 468},
  { /*   161 */ UINT64_C (0xE39C49765FDF9D94), 0xED5A7E86, 471},
  { /*   162 */ UINT64_C (0x8E41ADE9FBEBC27D), 0x14588F14, 475},
  { /*   163 */ UINT64_C (0xB1D219647AE6B31C), 0x596EB2D9, 478},
  { /*   164 */ UINT64_C (0xDE469FBD99A05FE3), 0x6FCA5F8F, 481},
  { /*   165 */ UINT64_C (0x8AEC23D680043BEE), 0x25DE7BB9, 485},
  { /*   166 */ UINT64_C (0xADA72CCC20054AE9), 0xAF561AA8, 488},
  { /*   167 */ UINT64_C (0xD910F7FF28069DA4), 0x1B2BA152, 491},
  { /*   168 */ UINT64_C (0x87AA9AFF79042286), 0x90FB44D3, 495},
  { /*   169 */ UINT64_C (0xA99541BF57452B28), 0x353A1608, 498},
  { /*   170 */ UINT64_C (0xD3FA922F2D1675F2), 0x42889B8A, 501},
  { /*   171 */ UINT64_C (0x847C9B5D7C2E09B7), 0x69956136, 505},
  { /*   172 */ UINT64_C (0xA59BC234DB398C25), 0x43FAB983, 508},
  { /*   173 */ UINT64_C (0xCF02B2C21207EF2E), 0x94F967E4, 511},
  { /*   174 */ UINT64_C (0x8161AFB94B44F57D), 0x1D1BE0EF, 515},
  { /*   175 */ UINT64_C (0xA1BA1BA79E1632DC), 0x6462D92A, 518},
  { /*   176 */ UINT64_C (0xCA28A291859BBF93), 0x7D7B8F75, 521},
  { /*   177 */ UINT64_C (0xFCB2CB35E702AF78), 0x5CDA7352, 524},
  { /*   178 */ UINT64_C (0x9DEFBF01B061ADAB), 0x3A088813, 528},
  { /*   179 */ UINT64_C (0xC56BAEC21C7A1916), 0x088AAA18, 531},
  { /*   180 */ UINT64_C (0xF6C69A72A3989F5B), 0x8AAD549E, 534},
  { /*   181 */ UINT64_C (0x9A3C2087A63F6399), 0x36AC54E3, 538},
  { /*   182 */ UINT64_C (0xC0CB28A98FCF3C7F), 0x84576A1C, 541},
  { /*   183 */ UINT64_C (0xF0FDF2D3F3C30B9F), 0x656D44A3, 544},
  { /*   184 */ UINT64_C (0x969EB7C47859E743), 0x9F644AE6, 548},
  { /*   185 */ UINT64_C (0xBC4665B596706114), 0x873D5D9F, 551},
  { /*   186 */ UINT64_C (0xEB57FF22FC0C7959), 0xA90CB507, 554},
  { /*   187 */ UINT64_C (0x9316FF75DD87CBD8), 0x09A7F124, 558},
  { /*   188 */ UINT64_C (0xB7DCBF5354E9BECE), 0x0C11ED6D, 561},
  { /*   189 */ UINT64_C (0xE5D3EF282A242E81), 0x8F1668C9, 564},
  { /*   190 */ UINT64_C (0x8FA475791A569D10), 0xF96E017D, 568},
  { /*   191 */ UINT64_C (0xB38D92D760EC4455), 0x37C981DD, 571},
  { /*   192 */ UINT64_C (0xE070F78D3927556A), 0x85BBE254, 574},
  { /*   193 */ UINT64_C (0x8C469AB843B89562), 0x93956D74, 578},
  { /*   194 */ UINT64_C (0xAF58416654A6BABB), 0x387AC8D2, 581},
  { /*   195 */ UINT64_C (0xDB2E51BFE9D0696A), 0x06997B06, 584},
  { /*   196 */ UINT64_C (0x88FCF317F22241E2), 0x441FECE4, 588},
  { /*   197 */ UINT64_C (0xAB3C2FDDEEAAD25A), 0xD527E81D, 591},
  { /*   198 */ UINT64_C (0xD60B3BD56A5586F1), 0x8A71E224, 594},
  { /*   199 */ UINT64_C (0x85C7056562757456), 0xF6872D56, 598},
  { /*   200 */ UINT64_C (0xA738C6BEBB12D16C), 0xB428F8AC, 601},
  { /*   201 */ UINT64_C (0xD106F86E69D785C7), 0xE13336D7, 604},
  { /*   202 */ UINT64_C (0x82A45B450226B39C), 0xECC00246, 608},
  { /*   203 */ UINT64_C (0xA34D721642B06084), 0x27F002D8, 611},
  { /*   204 */ UINT64_C (0xCC20CE9BD35C78A5), 0x31EC038E, 614},
  { /*   205 */ UINT64_C (0xFF290242C83396CE), 0x7E670471, 617},
  { /*   206 */ UINT64_C (0x9F79A169BD203E41), 0x0F0062C7, 621},
  { /*   207 */ UINT64_C (0xC75809C42C684DD1), 0x52C07B79, 624},
  { /*   208 */ UINT64_C (0xF92E0C3537826145), 0xA7709A57, 627},
  { /*   209 */ UINT64_C (0x9BBCC7A142B17CCB), 0x88A66076, 631},
  { /*   210 */ UINT64_C (0xC2ABF989935DDBFE), 0x6ACFF894, 634},
  { /*   211 */ UINT64_C (0xF356F7EBF83552FE), 0x0583F6B9, 637},
  { /*   212 */ UINT64_C (0x98165AF37B2153DE), 0xC3727A33, 641},
  { /*   213 */ UINT64_C (0xBE1BF1B059E9A8D6), 0x744F18C0, 644},
  { /*   214 */ UINT64_C (0xEDA2EE1C7064130C), 0x1162DEF0, 647},
  { /*   215 */ UINT64_C (0x9485D4D1C63E8BE7), 0x8ADDCB56, 651},
  { /*   216 */ UINT64_C (0xB9A74A0637CE2EE1), 0x6D953E2C, 654},
  { /*   217 */ UINT64_C (0xE8111C87C5C1BA99), 0xC8FA8DB7, 657},
  { /*   218 */ UINT64_C (0x910AB1D4DB9914A0), 0x1D9C9892, 661},
  { /*   219 */ UINT64_C (0xB54D5E4A127F59C8), 0x2503BEB7, 664},
  { /*   220 */ UINT64_C (0xE2A0B5DC971F303A), 0x2E44AE65, 667},
  { /*   221 */ UINT64_C (0x8DA471A9DE737E24), 0x5CEAECFF, 671},
  { /*   222 */ UINT64_C (0xB10D8E1456105DAD), 0x7425A83F, 674},
  { /*   223 */ UINT64_C (0xDD50F1996B947518), 0xD12F124E, 677},
  { /*   224 */ UINT64_C (0x8A5296FFE33CC92F), 0x82BD6B71, 681},
  { /*   225 */ UINT64_C (0xACE73CBFDC0BFB7B), 0x636CC64D, 684},
  { /*   226 */ UINT64_C (0xD8210BEFD30EFA5A), 0x3C47F7E0, 687},
  { /*   227 */ UINT64_C (0x8714A775E3E95C78), 0x65ACFAEC, 691},
  { /*   228 */ UINT64_C (0xA8D9D1535CE3B396), 0x7F1839A7, 694},
  { /*   229 */ UINT64_C (0xD31045A8341CA07C), 0x1EDE4811, 697},
  { /*   230 */ UINT64_C (0x83EA2B892091E44D), 0x934AED0B, 701},
  { /*   231 */ UINT64_C (0xA4E4B66B68B65D60), 0xF81DA84D, 704},
  { /*   232 */ UINT64_C (0xCE1DE40642E3F4B9), 0x36251261, 707},
  { /*   233 */ UINT64_C (0x80D2AE83E9CE78F3), 0xC1D72B7C, 711},
  { /*   234 */ UINT64_C (0xA1075A24E4421730), 0xB24CF65C, 714},
  { /*   235 */ UINT64_C (0xC94930AE1D529CFC), 0xDEE033F2, 717},
  { /*   236 */ UINT64_C (0xFB9B7CD9A4A7443C), 0x169840EF, 720},
  { /*   237 */ UINT64_C (0x9D412E0806E88AA5), 0x8E1F2895, 724},
  { /*   238 */ UINT64_C (0xC491798A08A2AD4E), 0xF1A6F2BB, 727},
  { /*   239 */ UINT64_C (0xF5B5D7EC8ACB58A2), 0xAE10AF69, 730},
  { /*   240 */ UINT64_C (0x9991A6F3D6BF1765), 0xACCA6DA2, 734},
  { /*   241 */ UINT64_C (0xBFF610B0CC6EDD3F), 0x17FD090A, 737},
  { /*   242 */ UINT64_C (0xEFF394DCFF8A948E), 0xDDFC4B4D, 740},
  { /*   243 */ UINT64_C (0x95F83D0A1FB69CD9), 0x4ABDAF10, 744},
  { /*   244 */ UINT64_C (0xBB764C4CA7A4440F), 0x9D6D1AD4, 747},
  { /*   245 */ UINT64_C (0xEA53DF5FD18D5513), 0x84C86189, 750},
  { /*   246 */ UINT64_C (0x92746B9BE2F8552C), 0x32FD3CF6, 754},
  { /*   247 */ UINT64_C (0xB7118682DBB66A77), 0x3FBC8C33, 757},
  { /*   248 */ UINT64_C (0xE4D5E82392A40515), 0x0FABAF40, 760},
  { /*   249 */ UINT64_C (0x8F05B1163BA6832D), 0x29CB4D88, 764},
  { /*   250 */ UINT64_C (0xB2C71D5BCA9023F8), 0x743E20EA, 767},
  { /*   251 */ UINT64_C (0xDF78E4B2BD342CF6), 0x914DA924, 770},
  { /*   252 */ UINT64_C (0x8BAB8EEFB6409C1A), 0x1AD089B7, 774},
  { /*   253 */ UINT64_C (0xAE9672ABA3D0C320), 0xA184AC24, 777},
  { /*   254 */ UINT64_C (0xDA3C0F568CC4F3E8), 0xC9E5D72E, 780},
  { /*   255 */ UINT64_C (0x8865899617FB1871), 0x7E2FA67C, 784},
  { /*   256 */ UINT64_C (0xAA7EEBFB9DF9DE8D), 0xDDBB901C, 787},
  { /*   257 */ UINT64_C (0xD51EA6FA85785631), 0x552A7422, 790},
  { /*   258 */ UINT64_C (0x8533285C936B35DE), 0xD53A8896, 794},
  { /*   259 */ UINT64_C (0xA67FF273B8460356), 0x8A892ABB, 797},
  { /*   260 */ UINT64_C (0xD01FEF10A657842C), 0x2D2B756A, 800},
  { /*   261 */ UINT64_C (0x8213F56A67F6B29B), 0x9C3B2962, 804},
  { /*   262 */ UINT64_C (0xA298F2C501F45F42), 0x8349F3BB, 807},
  { /*   263 */ UINT64_C (0xCB3F2F7642717713), 0x241C70A9, 810},
  { /*   264 */ UINT64_C (0xFE0EFB53D30DD4D7), 0xED238CD4, 813},
  { /*   265 */ UINT64_C (0x9EC95D1463E8A506), 0xF4363804, 817},
  { /*   266 */ UINT64_C (0xC67BB4597CE2CE48), 0xB143C605, 820},
  { /*   267 */ UINT64_C (0xF81AA16FDC1B81DA), 0xDD94B787, 823},
  { /*   268 */ UINT64_C (0x9B10A4E5E9913128), 0xCA7CF2B4, 827},
  { /*   269 */ UINT64_C (0xC1D4CE1F63F57D72), 0xFD1C2F61, 830},
  { /*   270 */ UINT64_C (0xF24A01A73CF2DCCF), 0xBC633B39, 833},
  { /*   271 */ UINT64_C (0x976E41088617CA01), 0xD5BE0504, 837},
  { /*   272 */ UINT64_C (0xBD49D14AA79DBC82), 0x4B2D8645, 840},
  { /*   273 */ UINT64_C (0xEC9C459D51852BA2), 0xDDF8E7D6, 843},
  { /*   274 */ UINT64_C (0x93E1AB8252F33B45), 0xCABB90E6, 847},
  { /*   275 */ UINT64_C (0xB8DA1662E7B00A17), 0x3D6A751F, 850},
  { /*   276 */ UINT64_C (0xE7109BFBA19C0C9D), 0x0CC51267, 853},
  { /*   277 */ UINT64_C (0x906A617D450187E2), 0x27FB2B80, 857},
  { /*   278 */ UINT64_C (0xB484F9DC9641E9DA), 0xB1F9F661, 860},
  { /*   279 */ UINT64_C (0xE1A63853BBD26451), 0x5E7873F9, 863},
  { /*   280 */ UINT64_C (0x8D07E33455637EB2), 0xDB0B487B, 867},
  { /*   281 */ UINT64_C (0xB049DC016ABC5E5F), 0x91CE1A9A, 870},
  { /*   282 */ UINT64_C (0xDC5C5301C56B75F7), 0x7641A141, 873},
  { /*   283 */ UINT64_C (0x89B9B3E11B6329BA), 0xA9E904C8, 877},
  { /*   284 */ UINT64_C (0xAC2820D9623BF429), 0x546345FB, 880},
  { /*   285 */ UINT64_C (0xD732290FBACAF133), 0xA97C1779, 883},
  { /*   286 */ UINT64_C (0x867F59A9D4BED6C0), 0x49ED8EAC, 887},
  { /*   287 */ UINT64_C (0xA81F301449EE8C70), 0x5C68F257, 890},
  { /*   288 */ UINT64_C (0xD226FC195C6A2F8C), 0x73832EEC, 893},
  { /*   289 */ UINT64_C (0x83585D8FD9C25DB7), 0xC831FD54, 897},
  { /*   290 */ UINT64_C (0xA42E74F3D032F525), 0xBA3E7CA9, 900},
  { /*   291 */ UINT64_C (0xCD3A1230C43FB26F), 0x28CE1BD3, 903},
  { /*   292 */ UINT64_C (0x80444B5E7AA7CF85), 0x7980D164, 907},
  { /*   293 */ UINT64_C (0xA0555E361951C366), 0xD7E105BD, 910},
  { /*   294 */ UINT64_C (0xC86AB5C39FA63440), 0x8DD9472C, 913},
  { /*   295 */ UINT64_C (0xFA856334878FC150), 0xB14F98F7, 916},
  { /*   296 */ UINT64_C (0x9C935E00D4B9D8D2), 0x6ED1BF9A, 920},
  { /*   297 */ UINT64_C (0xC3B8358109E84F07), 0x0A862F81, 923},
  { /*   298 */ UINT64_C (0xF4A642E14C6262C8), 0xCD27BB61, 926},
  { /*   299 */ UINT64_C (0x98E7E9CCCFBD7DBD), 0x8038D51D, 930},
  { /*   300 */ UINT64_C (0xBF21E44003ACDD2C), 0xE0470A64, 933},
  { /*   301 */ UINT64_C (0xEEEA5D5004981478), 0x1858CCFD, 936},
  { /*   302 */ UINT64_C (0x95527A5202DF0CCB), 0x0F37801E, 940},
  { /*   303 */ UINT64_C (0xBAA718E68396CFFD), 0xD3056026, 943},
  { /*   304 */ UINT64_C (0xE950DF20247C83FD), 0x47C6B82F, 946},
  { /*   305 */ UINT64_C (0x91D28B7416CDD27E), 0x4CDC331D, 950},
  { /*   306 */ UINT64_C (0xB6472E511C81471D), 0xE0133FE5, 953},
  { /*   307 */ UINT64_C (0xE3D8F9E563A198E5), 0x58180FDE, 956},
  { /*   308 */ UINT64_C (0x8E679C2F5E44FF8F), 0x570F09EB, 960},
  { /*   309 */ UINT64_C (0xB201833B35D63F73), 0x2CD2CC65, 963},
};

static const struct
{
  uint64_t mul;
  int32_t exp;
} fpowers10[] = {
  { /*   -64 */ UINT64_C (0xA87FEA27A539E9A5), -244},
  { /*   -63 */ UINT64_C (0xD29FE4B18E88640F), -241},
  { /*   -62 */ UINT64_C (0x83A3EEEEF9153E89), -237},
  { /*   -61 */ UINT64_C (0xA48CEAAAB75A8E2B), -234},
  { /*   -60 */ UINT64_C (0xCDB02555653131B6), -231},
  { /*   -59 */ UINT64_C (0x808E17555F3EBF12), -227},
  { /*   -58 */ UINT64_C (0xA0B19D2AB70E6ED6), -224},
  { /*   -57 */ UINT64_C (0xC8DE047564D20A8C), -221},
  { /*   -56 */ UINT64_C (0xFB158592BE068D2F), -218},
  { /*   -55 */ UINT64_C (0x9CED737BB6C4183D), -214},
  { /*   -54 */ UINT64_C (0xC428D05AA4751E4D), -211},
  { /*   -53 */ UINT64_C (0xF53304714D9265E0), -208},
  { /*   -52 */ UINT64_C (0x993FE2C6D07B7FAC), -204},
  { /*   -51 */ UINT64_C (0xBF8FDB78849A5F97), -201},
  { /*   -50 */ UINT64_C (0xEF73D256A5C0F77D), -198},
  { /*   -49 */ UINT64_C (0x95A8637627989AAE), -194},
  { /*   -48 */ UINT64_C (0xBB127C53B17EC159), -191},
  { /*   -47 */ UINT64_C (0xE9D71B689DDE71B0), -188},
  { /*   -46 */ UINT64_C (0x9226712162AB070E), -184},
  { /*   -45 */ UINT64_C (0xB6B00D69BB55C8D1), -181},
  { /*   -44 */ UINT64_C (0xE45C10C42A2B3B06), -178},
  { /*   -43 */ UINT64_C (0x8EB98A7A9A5B04E3), -174},
  { /*   -42 */ UINT64_C (0xB267ED1940F1C61C), -171},
  { /*   -41 */ UINT64_C (0xDF01E85F912E37A3), -168},
  { /*   -40 */ UINT64_C (0x8B61313BBABCE2C6), -164},
  { /*   -39 */ UINT64_C (0xAE397D8AA96C1B78), -161},
  { /*   -38 */ UINT64_C (0xD9C7DCED53C72256), -158},
  { /*   -37 */ UINT64_C (0x881CEA14545C7575), -154},
  { /*   -36 */ UINT64_C (0xAA242499697392D3), -151},
  { /*   -35 */ UINT64_C (0xD4AD2DBFC3D07788), -148},
  { /*   -34 */ UINT64_C (0x84EC3C97DA624AB5), -144},
  { /*   -33 */ UINT64_C (0xA6274BBDD0FADD62), -141},
  { /*   -32 */ UINT64_C (0xCFB11EAD453994BA), -138},
  { /*   -31 */ UINT64_C (0x81CEB32C4B43FCF5), -134},
  { /*   -30 */ UINT64_C (0xA2425FF75E14FC32), -131},
  { /*   -29 */ UINT64_C (0xCAD2F7F5359A3B3E), -128},
  { /*   -28 */ UINT64_C (0xFD87B5F28300CA0E), -125},
  { /*   -27 */ UINT64_C (0x9E74D1B791E07E48), -121},
  { /*   -26 */ UINT64_C (0xC612062576589DDB), -118},
  { /*   -25 */ UINT64_C (0xF79687AED3EEC551), -115},
  { /*   -24 */ UINT64_C (0x9ABE14CD44753B53), -111},
  { /*   -23 */ UINT64_C (0xC16D9A0095928A27), -108},
  { /*   -22 */ UINT64_C (0xF1C90080BAF72CB1), -105},
  { /*   -21 */ UINT64_C (0x971DA05074DA7BEF), -101},
  { /*   -20 */ UINT64_C (0xBCE5086492111AEB), -98},
  { /*   -19 */ UINT64_C (0xEC1E4A7DB69561A5), -95},
  { /*   -18 */ UINT64_C (0x9392EE8E921D5D07), -91},
  { /*   -17 */ UINT64_C (0xB877AA3236A4B449), -88},
  { /*   -16 */ UINT64_C (0xE69594BEC44DE15B), -85},
  { /*   -15 */ UINT64_C (0x901D7CF73AB0ACD9), -81},
  { /*   -14 */ UINT64_C (0xB424DC35095CD80F), -78},
  { /*   -13 */ UINT64_C (0xE12E13424BB40E13), -75},
  { /*   -12 */ UINT64_C (0x8CBCCC096F5088CC), -71},
  { /*   -11 */ UINT64_C (0xAFEBFF0BCB24AAFF), -68},
  { /*   -10 */ UINT64_C (0xDBE6FECEBDEDD5BF), -65},
  { /*    -9 */ UINT64_C (0x89705F4136B4A597), -61},
  { /*    -8 */ UINT64_C (0xABCC77118461CEFD), -58},
  { /*    -7 */ UINT64_C (0xD6BF94D5E57A42BC), -55},
  { /*    -6 */ UINT64_C (0x8637BD05AF6C69B6), -51},
  { /*    -5 */ UINT64_C (0xA7C5AC471B478423), -48},
  { /*    -4 */ UINT64_C (0xD1B71758E219652C), -45},
  { /*    -3 */ UINT64_C (0x83126E978D4FDF3B), -41},
  { /*    -2 */ UINT64_C (0xA3D70A3D70A3D70A), -38},
  { /*    -1 */ UINT64_C (0xCCCCCCCCCCCCCCCD), -35},
  { /*     0 */ UINT64_C (0x8000000000000000), -31},
  { /*     1 */ UINT64_C (0xA000000000000000), -28},
  { /*     2 */ UINT64_C (0xC800000000000000), -25},
  { /*     3 */ UINT64_C (0xFA00000000000000), -22},
  { /*     4 */ UINT64_C (0x9C40000000000000), -18},
  { /*     5 */ UINT64_C (0xC350000000000000), -15},
  { /*     6 */ UINT64_C (0xF424000000000000), -12},
  { /*     7 */ UINT64_C (0x9896800000000000), -8},
  { /*     8 */ UINT64_C (0xBEBC200000000000), -5},
  { /*     9 */ UINT64_C (0xEE6B280000000000), -2},
  { /*    10 */ UINT64_C (0x9502F90000000000), 2},
  { /*    11 */ UINT64_C (0xBA43B74000000000), 5},
  { /*    12 */ UINT64_C (0xE8D4A51000000000), 8},
  { /*    13 */ UINT64_C (0x9184E72A00000000), 12},
  { /*    14 */ UINT64_C (0xB5E620F480000000), 15},
  { /*    15 */ UINT64_C (0xE35FA931A0000000), 18},
  { /*    16 */ UINT64_C (0x8E1BC9BF04000000), 22},
  { /*    17 */ UINT64_C (0xB1A2BC2EC5000000), 25},
  { /*    18 */ UINT64_C (0xDE0B6B3A76400000), 28},
  { /*    19 */ UINT64_C (0x8AC7230489E80000), 32},
  { /*    20 */ UINT64_C (0xAD78EBC5AC620000), 35},
  { /*    21 */ UINT64_C (0xD8D726B7177A8000), 38},
  { /*    22 */ UINT64_C (0x878678326EAC9000), 42},
  { /*    23 */ UINT64_C (0xA968163F0A57B400), 45},
  { /*    24 */ UINT64_C (0xD3C21BCECCEDA100), 48},
  { /*    25 */ UINT64_C (0x84595161401484A0), 52},
  { /*    26 */ UINT64_C (0xA56FA5B99019A5C8), 55},
  { /*    27 */ UINT64_C (0xCECB8F27F4200F3A), 58},
  { /*    28 */ UINT64_C (0x813F3978F8940984), 62},
  { /*    29 */ UINT64_C (0xA18F07D736B90BE5), 65},
  { /*    30 */ UINT64_C (0xC9F2C9CD04674EDF), 68},
  { /*    31 */ UINT64_C (0xFC6F7C4045812296), 71},
  { /*    32 */ UINT64_C (0x9DC5ADA82B70B59E), 75},
  { /*    33 */ UINT64_C (0xC5371912364CE305), 78},
  { /*    34 */ UINT64_C (0xF684DF56C3E01BC7), 81},
  { /*    35 */ UINT64_C (0x9A130B963A6C115C), 85},
  { /*    36 */ UINT64_C (0xC097CE7BC90715B3), 88},
  { /*    37 */ UINT64_C (0xF0BDC21ABB48DB20), 91},
  { /*    38 */ UINT64_C (0x96769950B50D88F4), 95},
  { /*    39 */ UINT64_C (0xBC143FA4E250EB31), 98},
};

static const uint64_t ipowers64[] = {
  UINT64_C (1),
  UINT64_C (10),
  UINT64_C (100),
  UINT64_C (1000),
  UINT64_C (10000),
  UINT64_C (100000),
  UINT64_C (1000000),
  UINT64_C (10000000),
  UINT64_C (100000000),
  UINT64_C (1000000000),
  UINT64_C (10000000000),
  UINT64_C (100000000000),
  UINT64_C (1000000000000),
  UINT64_C (10000000000000),
  UINT64_C (100000000000000),
  UINT64_C (1000000000000000),
  UINT64_C (10000000000000000),
  UINT64_C (100000000000000000),
  UINT64_C (1000000000000000000),
  UINT64_C (10000000000000000000),
};

static const uint32_t ipowers32[] = {
  1,
  10,
  100,
  1000,
  10000,
  100000,
  1000000,
  10000000,
  100000000,
  1000000000,
};

#if 0
/* gcc -g -O3 -Wall c.c -o c */
#include <stdio.h>

int
main (void)
{
  int t1;
  int t2;
  int t3;

  printf ("static const char num3[3000] = {\n");
  for (t1 = 0; t1 < 10; t1++) {
    for (t2 = 0; t2 < 10; t2++) {
      printf ("  ");
      for (t3 = 0; t3 < 5; t3++) {
	printf ("'%d', '%d', '%d', ", t1, t2, t3);
      }
      printf ("\n  ");
      for (; t3 < 10; t3++) {
	printf ("'%d', '%d', '%d', ", t1, t2, t3);
      }
      printf ("\n");
    }
  }
  printf ("};\n");
  return 0;
}
#endif

static const char num3[3000] = {
  '0', '0', '0', '0', '0', '1', '0', '0', '2', '0', '0', '3', '0', '0', '4',
  '0', '0', '5', '0', '0', '6', '0', '0', '7', '0', '0', '8', '0', '0', '9',
  '0', '1', '0', '0', '1', '1', '0', '1', '2', '0', '1', '3', '0', '1', '4',
  '0', '1', '5', '0', '1', '6', '0', '1', '7', '0', '1', '8', '0', '1', '9',
  '0', '2', '0', '0', '2', '1', '0', '2', '2', '0', '2', '3', '0', '2', '4',
  '0', '2', '5', '0', '2', '6', '0', '2', '7', '0', '2', '8', '0', '2', '9',
  '0', '3', '0', '0', '3', '1', '0', '3', '2', '0', '3', '3', '0', '3', '4',
  '0', '3', '5', '0', '3', '6', '0', '3', '7', '0', '3', '8', '0', '3', '9',
  '0', '4', '0', '0', '4', '1', '0', '4', '2', '0', '4', '3', '0', '4', '4',
  '0', '4', '5', '0', '4', '6', '0', '4', '7', '0', '4', '8', '0', '4', '9',
  '0', '5', '0', '0', '5', '1', '0', '5', '2', '0', '5', '3', '0', '5', '4',
  '0', '5', '5', '0', '5', '6', '0', '5', '7', '0', '5', '8', '0', '5', '9',
  '0', '6', '0', '0', '6', '1', '0', '6', '2', '0', '6', '3', '0', '6', '4',
  '0', '6', '5', '0', '6', '6', '0', '6', '7', '0', '6', '8', '0', '6', '9',
  '0', '7', '0', '0', '7', '1', '0', '7', '2', '0', '7', '3', '0', '7', '4',
  '0', '7', '5', '0', '7', '6', '0', '7', '7', '0', '7', '8', '0', '7', '9',
  '0', '8', '0', '0', '8', '1', '0', '8', '2', '0', '8', '3', '0', '8', '4',
  '0', '8', '5', '0', '8', '6', '0', '8', '7', '0', '8', '8', '0', '8', '9',
  '0', '9', '0', '0', '9', '1', '0', '9', '2', '0', '9', '3', '0', '9', '4',
  '0', '9', '5', '0', '9', '6', '0', '9', '7', '0', '9', '8', '0', '9', '9',
  '1', '0', '0', '1', '0', '1', '1', '0', '2', '1', '0', '3', '1', '0', '4',
  '1', '0', '5', '1', '0', '6', '1', '0', '7', '1', '0', '8', '1', '0', '9',
  '1', '1', '0', '1', '1', '1', '1', '1', '2', '1', '1', '3', '1', '1', '4',
  '1', '1', '5', '1', '1', '6', '1', '1', '7', '1', '1', '8', '1', '1', '9',
  '1', '2', '0', '1', '2', '1', '1', '2', '2', '1', '2', '3', '1', '2', '4',
  '1', '2', '5', '1', '2', '6', '1', '2', '7', '1', '2', '8', '1', '2', '9',
  '1', '3', '0', '1', '3', '1', '1', '3', '2', '1', '3', '3', '1', '3', '4',
  '1', '3', '5', '1', '3', '6', '1', '3', '7', '1', '3', '8', '1', '3', '9',
  '1', '4', '0', '1', '4', '1', '1', '4', '2', '1', '4', '3', '1', '4', '4',
  '1', '4', '5', '1', '4', '6', '1', '4', '7', '1', '4', '8', '1', '4', '9',
  '1', '5', '0', '1', '5', '1', '1', '5', '2', '1', '5', '3', '1', '5', '4',
  '1', '5', '5', '1', '5', '6', '1', '5', '7', '1', '5', '8', '1', '5', '9',
  '1', '6', '0', '1', '6', '1', '1', '6', '2', '1', '6', '3', '1', '6', '4',
  '1', '6', '5', '1', '6', '6', '1', '6', '7', '1', '6', '8', '1', '6', '9',
  '1', '7', '0', '1', '7', '1', '1', '7', '2', '1', '7', '3', '1', '7', '4',
  '1', '7', '5', '1', '7', '6', '1', '7', '7', '1', '7', '8', '1', '7', '9',
  '1', '8', '0', '1', '8', '1', '1', '8', '2', '1', '8', '3', '1', '8', '4',
  '1', '8', '5', '1', '8', '6', '1', '8', '7', '1', '8', '8', '1', '8', '9',
  '1', '9', '0', '1', '9', '1', '1', '9', '2', '1', '9', '3', '1', '9', '4',
  '1', '9', '5', '1', '9', '6', '1', '9', '7', '1', '9', '8', '1', '9', '9',
  '2', '0', '0', '2', '0', '1', '2', '0', '2', '2', '0', '3', '2', '0', '4',
  '2', '0', '5', '2', '0', '6', '2', '0', '7', '2', '0', '8', '2', '0', '9',
  '2', '1', '0', '2', '1', '1', '2', '1', '2', '2', '1', '3', '2', '1', '4',
  '2', '1', '5', '2', '1', '6', '2', '1', '7', '2', '1', '8', '2', '1', '9',
  '2', '2', '0', '2', '2', '1', '2', '2', '2', '2', '2', '3', '2', '2', '4',
  '2', '2', '5', '2', '2', '6', '2', '2', '7', '2', '2', '8', '2', '2', '9',
  '2', '3', '0', '2', '3', '1', '2', '3', '2', '2', '3', '3', '2', '3', '4',
  '2', '3', '5', '2', '3', '6', '2', '3', '7', '2', '3', '8', '2', '3', '9',
  '2', '4', '0', '2', '4', '1', '2', '4', '2', '2', '4', '3', '2', '4', '4',
  '2', '4', '5', '2', '4', '6', '2', '4', '7', '2', '4', '8', '2', '4', '9',
  '2', '5', '0', '2', '5', '1', '2', '5', '2', '2', '5', '3', '2', '5', '4',
  '2', '5', '5', '2', '5', '6', '2', '5', '7', '2', '5', '8', '2', '5', '9',
  '2', '6', '0', '2', '6', '1', '2', '6', '2', '2', '6', '3', '2', '6', '4',
  '2', '6', '5', '2', '6', '6', '2', '6', '7', '2', '6', '8', '2', '6', '9',
  '2', '7', '0', '2', '7', '1', '2', '7', '2', '2', '7', '3', '2', '7', '4',
  '2', '7', '5', '2', '7', '6', '2', '7', '7', '2', '7', '8', '2', '7', '9',
  '2', '8', '0', '2', '8', '1', '2', '8', '2', '2', '8', '3', '2', '8', '4',
  '2', '8', '5', '2', '8', '6', '2', '8', '7', '2', '8', '8', '2', '8', '9',
  '2', '9', '0', '2', '9', '1', '2', '9', '2', '2', '9', '3', '2', '9', '4',
  '2', '9', '5', '2', '9', '6', '2', '9', '7', '2', '9', '8', '2', '9', '9',
  '3', '0', '0', '3', '0', '1', '3', '0', '2', '3', '0', '3', '3', '0', '4',
  '3', '0', '5', '3', '0', '6', '3', '0', '7', '3', '0', '8', '3', '0', '9',
  '3', '1', '0', '3', '1', '1', '3', '1', '2', '3', '1', '3', '3', '1', '4',
  '3', '1', '5', '3', '1', '6', '3', '1', '7', '3', '1', '8', '3', '1', '9',
  '3', '2', '0', '3', '2', '1', '3', '2', '2', '3', '2', '3', '3', '2', '4',
  '3', '2', '5', '3', '2', '6', '3', '2', '7', '3', '2', '8', '3', '2', '9',
  '3', '3', '0', '3', '3', '1', '3', '3', '2', '3', '3', '3', '3', '3', '4',
  '3', '3', '5', '3', '3', '6', '3', '3', '7', '3', '3', '8', '3', '3', '9',
  '3', '4', '0', '3', '4', '1', '3', '4', '2', '3', '4', '3', '3', '4', '4',
  '3', '4', '5', '3', '4', '6', '3', '4', '7', '3', '4', '8', '3', '4', '9',
  '3', '5', '0', '3', '5', '1', '3', '5', '2', '3', '5', '3', '3', '5', '4',
  '3', '5', '5', '3', '5', '6', '3', '5', '7', '3', '5', '8', '3', '5', '9',
  '3', '6', '0', '3', '6', '1', '3', '6', '2', '3', '6', '3', '3', '6', '4',
  '3', '6', '5', '3', '6', '6', '3', '6', '7', '3', '6', '8', '3', '6', '9',
  '3', '7', '0', '3', '7', '1', '3', '7', '2', '3', '7', '3', '3', '7', '4',
  '3', '7', '5', '3', '7', '6', '3', '7', '7', '3', '7', '8', '3', '7', '9',
  '3', '8', '0', '3', '8', '1', '3', '8', '2', '3', '8', '3', '3', '8', '4',
  '3', '8', '5', '3', '8', '6', '3', '8', '7', '3', '8', '8', '3', '8', '9',
  '3', '9', '0', '3', '9', '1', '3', '9', '2', '3', '9', '3', '3', '9', '4',
  '3', '9', '5', '3', '9', '6', '3', '9', '7', '3', '9', '8', '3', '9', '9',
  '4', '0', '0', '4', '0', '1', '4', '0', '2', '4', '0', '3', '4', '0', '4',
  '4', '0', '5', '4', '0', '6', '4', '0', '7', '4', '0', '8', '4', '0', '9',
  '4', '1', '0', '4', '1', '1', '4', '1', '2', '4', '1', '3', '4', '1', '4',
  '4', '1', '5', '4', '1', '6', '4', '1', '7', '4', '1', '8', '4', '1', '9',
  '4', '2', '0', '4', '2', '1', '4', '2', '2', '4', '2', '3', '4', '2', '4',
  '4', '2', '5', '4', '2', '6', '4', '2', '7', '4', '2', '8', '4', '2', '9',
  '4', '3', '0', '4', '3', '1', '4', '3', '2', '4', '3', '3', '4', '3', '4',
  '4', '3', '5', '4', '3', '6', '4', '3', '7', '4', '3', '8', '4', '3', '9',
  '4', '4', '0', '4', '4', '1', '4', '4', '2', '4', '4', '3', '4', '4', '4',
  '4', '4', '5', '4', '4', '6', '4', '4', '7', '4', '4', '8', '4', '4', '9',
  '4', '5', '0', '4', '5', '1', '4', '5', '2', '4', '5', '3', '4', '5', '4',
  '4', '5', '5', '4', '5', '6', '4', '5', '7', '4', '5', '8', '4', '5', '9',
  '4', '6', '0', '4', '6', '1', '4', '6', '2', '4', '6', '3', '4', '6', '4',
  '4', '6', '5', '4', '6', '6', '4', '6', '7', '4', '6', '8', '4', '6', '9',
  '4', '7', '0', '4', '7', '1', '4', '7', '2', '4', '7', '3', '4', '7', '4',
  '4', '7', '5', '4', '7', '6', '4', '7', '7', '4', '7', '8', '4', '7', '9',
  '4', '8', '0', '4', '8', '1', '4', '8', '2', '4', '8', '3', '4', '8', '4',
  '4', '8', '5', '4', '8', '6', '4', '8', '7', '4', '8', '8', '4', '8', '9',
  '4', '9', '0', '4', '9', '1', '4', '9', '2', '4', '9', '3', '4', '9', '4',
  '4', '9', '5', '4', '9', '6', '4', '9', '7', '4', '9', '8', '4', '9', '9',
  '5', '0', '0', '5', '0', '1', '5', '0', '2', '5', '0', '3', '5', '0', '4',
  '5', '0', '5', '5', '0', '6', '5', '0', '7', '5', '0', '8', '5', '0', '9',
  '5', '1', '0', '5', '1', '1', '5', '1', '2', '5', '1', '3', '5', '1', '4',
  '5', '1', '5', '5', '1', '6', '5', '1', '7', '5', '1', '8', '5', '1', '9',
  '5', '2', '0', '5', '2', '1', '5', '2', '2', '5', '2', '3', '5', '2', '4',
  '5', '2', '5', '5', '2', '6', '5', '2', '7', '5', '2', '8', '5', '2', '9',
  '5', '3', '0', '5', '3', '1', '5', '3', '2', '5', '3', '3', '5', '3', '4',
  '5', '3', '5', '5', '3', '6', '5', '3', '7', '5', '3', '8', '5', '3', '9',
  '5', '4', '0', '5', '4', '1', '5', '4', '2', '5', '4', '3', '5', '4', '4',
  '5', '4', '5', '5', '4', '6', '5', '4', '7', '5', '4', '8', '5', '4', '9',
  '5', '5', '0', '5', '5', '1', '5', '5', '2', '5', '5', '3', '5', '5', '4',
  '5', '5', '5', '5', '5', '6', '5', '5', '7', '5', '5', '8', '5', '5', '9',
  '5', '6', '0', '5', '6', '1', '5', '6', '2', '5', '6', '3', '5', '6', '4',
  '5', '6', '5', '5', '6', '6', '5', '6', '7', '5', '6', '8', '5', '6', '9',
  '5', '7', '0', '5', '7', '1', '5', '7', '2', '5', '7', '3', '5', '7', '4',
  '5', '7', '5', '5', '7', '6', '5', '7', '7', '5', '7', '8', '5', '7', '9',
  '5', '8', '0', '5', '8', '1', '5', '8', '2', '5', '8', '3', '5', '8', '4',
  '5', '8', '5', '5', '8', '6', '5', '8', '7', '5', '8', '8', '5', '8', '9',
  '5', '9', '0', '5', '9', '1', '5', '9', '2', '5', '9', '3', '5', '9', '4',
  '5', '9', '5', '5', '9', '6', '5', '9', '7', '5', '9', '8', '5', '9', '9',
  '6', '0', '0', '6', '0', '1', '6', '0', '2', '6', '0', '3', '6', '0', '4',
  '6', '0', '5', '6', '0', '6', '6', '0', '7', '6', '0', '8', '6', '0', '9',
  '6', '1', '0', '6', '1', '1', '6', '1', '2', '6', '1', '3', '6', '1', '4',
  '6', '1', '5', '6', '1', '6', '6', '1', '7', '6', '1', '8', '6', '1', '9',
  '6', '2', '0', '6', '2', '1', '6', '2', '2', '6', '2', '3', '6', '2', '4',
  '6', '2', '5', '6', '2', '6', '6', '2', '7', '6', '2', '8', '6', '2', '9',
  '6', '3', '0', '6', '3', '1', '6', '3', '2', '6', '3', '3', '6', '3', '4',
  '6', '3', '5', '6', '3', '6', '6', '3', '7', '6', '3', '8', '6', '3', '9',
  '6', '4', '0', '6', '4', '1', '6', '4', '2', '6', '4', '3', '6', '4', '4',
  '6', '4', '5', '6', '4', '6', '6', '4', '7', '6', '4', '8', '6', '4', '9',
  '6', '5', '0', '6', '5', '1', '6', '5', '2', '6', '5', '3', '6', '5', '4',
  '6', '5', '5', '6', '5', '6', '6', '5', '7', '6', '5', '8', '6', '5', '9',
  '6', '6', '0', '6', '6', '1', '6', '6', '2', '6', '6', '3', '6', '6', '4',
  '6', '6', '5', '6', '6', '6', '6', '6', '7', '6', '6', '8', '6', '6', '9',
  '6', '7', '0', '6', '7', '1', '6', '7', '2', '6', '7', '3', '6', '7', '4',
  '6', '7', '5', '6', '7', '6', '6', '7', '7', '6', '7', '8', '6', '7', '9',
  '6', '8', '0', '6', '8', '1', '6', '8', '2', '6', '8', '3', '6', '8', '4',
  '6', '8', '5', '6', '8', '6', '6', '8', '7', '6', '8', '8', '6', '8', '9',
  '6', '9', '0', '6', '9', '1', '6', '9', '2', '6', '9', '3', '6', '9', '4',
  '6', '9', '5', '6', '9', '6', '6', '9', '7', '6', '9', '8', '6', '9', '9',
  '7', '0', '0', '7', '0', '1', '7', '0', '2', '7', '0', '3', '7', '0', '4',
  '7', '0', '5', '7', '0', '6', '7', '0', '7', '7', '0', '8', '7', '0', '9',
  '7', '1', '0', '7', '1', '1', '7', '1', '2', '7', '1', '3', '7', '1', '4',
  '7', '1', '5', '7', '1', '6', '7', '1', '7', '7', '1', '8', '7', '1', '9',
  '7', '2', '0', '7', '2', '1', '7', '2', '2', '7', '2', '3', '7', '2', '4',
  '7', '2', '5', '7', '2', '6', '7', '2', '7', '7', '2', '8', '7', '2', '9',
  '7', '3', '0', '7', '3', '1', '7', '3', '2', '7', '3', '3', '7', '3', '4',
  '7', '3', '5', '7', '3', '6', '7', '3', '7', '7', '3', '8', '7', '3', '9',
  '7', '4', '0', '7', '4', '1', '7', '4', '2', '7', '4', '3', '7', '4', '4',
  '7', '4', '5', '7', '4', '6', '7', '4', '7', '7', '4', '8', '7', '4', '9',
  '7', '5', '0', '7', '5', '1', '7', '5', '2', '7', '5', '3', '7', '5', '4',
  '7', '5', '5', '7', '5', '6', '7', '5', '7', '7', '5', '8', '7', '5', '9',
  '7', '6', '0', '7', '6', '1', '7', '6', '2', '7', '6', '3', '7', '6', '4',
  '7', '6', '5', '7', '6', '6', '7', '6', '7', '7', '6', '8', '7', '6', '9',
  '7', '7', '0', '7', '7', '1', '7', '7', '2', '7', '7', '3', '7', '7', '4',
  '7', '7', '5', '7', '7', '6', '7', '7', '7', '7', '7', '8', '7', '7', '9',
  '7', '8', '0', '7', '8', '1', '7', '8', '2', '7', '8', '3', '7', '8', '4',
  '7', '8', '5', '7', '8', '6', '7', '8', '7', '7', '8', '8', '7', '8', '9',
  '7', '9', '0', '7', '9', '1', '7', '9', '2', '7', '9', '3', '7', '9', '4',
  '7', '9', '5', '7', '9', '6', '7', '9', '7', '7', '9', '8', '7', '9', '9',
  '8', '0', '0', '8', '0', '1', '8', '0', '2', '8', '0', '3', '8', '0', '4',
  '8', '0', '5', '8', '0', '6', '8', '0', '7', '8', '0', '8', '8', '0', '9',
  '8', '1', '0', '8', '1', '1', '8', '1', '2', '8', '1', '3', '8', '1', '4',
  '8', '1', '5', '8', '1', '6', '8', '1', '7', '8', '1', '8', '8', '1', '9',
  '8', '2', '0', '8', '2', '1', '8', '2', '2', '8', '2', '3', '8', '2', '4',
  '8', '2', '5', '8', '2', '6', '8', '2', '7', '8', '2', '8', '8', '2', '9',
  '8', '3', '0', '8', '3', '1', '8', '3', '2', '8', '3', '3', '8', '3', '4',
  '8', '3', '5', '8', '3', '6', '8', '3', '7', '8', '3', '8', '8', '3', '9',
  '8', '4', '0', '8', '4', '1', '8', '4', '2', '8', '4', '3', '8', '4', '4',
  '8', '4', '5', '8', '4', '6', '8', '4', '7', '8', '4', '8', '8', '4', '9',
  '8', '5', '0', '8', '5', '1', '8', '5', '2', '8', '5', '3', '8', '5', '4',
  '8', '5', '5', '8', '5', '6', '8', '5', '7', '8', '5', '8', '8', '5', '9',
  '8', '6', '0', '8', '6', '1', '8', '6', '2', '8', '6', '3', '8', '6', '4',
  '8', '6', '5', '8', '6', '6', '8', '6', '7', '8', '6', '8', '8', '6', '9',
  '8', '7', '0', '8', '7', '1', '8', '7', '2', '8', '7', '3', '8', '7', '4',
  '8', '7', '5', '8', '7', '6', '8', '7', '7', '8', '7', '8', '8', '7', '9',
  '8', '8', '0', '8', '8', '1', '8', '8', '2', '8', '8', '3', '8', '8', '4',
  '8', '8', '5', '8', '8', '6', '8', '8', '7', '8', '8', '8', '8', '8', '9',
  '8', '9', '0', '8', '9', '1', '8', '9', '2', '8', '9', '3', '8', '9', '4',
  '8', '9', '5', '8', '9', '6', '8', '9', '7', '8', '9', '8', '8', '9', '9',
  '9', '0', '0', '9', '0', '1', '9', '0', '2', '9', '0', '3', '9', '0', '4',
  '9', '0', '5', '9', '0', '6', '9', '0', '7', '9', '0', '8', '9', '0', '9',
  '9', '1', '0', '9', '1', '1', '9', '1', '2', '9', '1', '3', '9', '1', '4',
  '9', '1', '5', '9', '1', '6', '9', '1', '7', '9', '1', '8', '9', '1', '9',
  '9', '2', '0', '9', '2', '1', '9', '2', '2', '9', '2', '3', '9', '2', '4',
  '9', '2', '5', '9', '2', '6', '9', '2', '7', '9', '2', '8', '9', '2', '9',
  '9', '3', '0', '9', '3', '1', '9', '3', '2', '9', '3', '3', '9', '3', '4',
  '9', '3', '5', '9', '3', '6', '9', '3', '7', '9', '3', '8', '9', '3', '9',
  '9', '4', '0', '9', '4', '1', '9', '4', '2', '9', '4', '3', '9', '4', '4',
  '9', '4', '5', '9', '4', '6', '9', '4', '7', '9', '4', '8', '9', '4', '9',
  '9', '5', '0', '9', '5', '1', '9', '5', '2', '9', '5', '3', '9', '5', '4',
  '9', '5', '5', '9', '5', '6', '9', '5', '7', '9', '5', '8', '9', '5', '9',
  '9', '6', '0', '9', '6', '1', '9', '6', '2', '9', '6', '3', '9', '6', '4',
  '9', '6', '5', '9', '6', '6', '9', '6', '7', '9', '6', '8', '9', '6', '9',
  '9', '7', '0', '9', '7', '1', '9', '7', '2', '9', '7', '3', '9', '7', '4',
  '9', '7', '5', '9', '7', '6', '9', '7', '7', '9', '7', '8', '9', '7', '9',
  '9', '8', '0', '9', '8', '1', '9', '8', '2', '9', '8', '3', '9', '8', '4',
  '9', '8', '5', '9', '8', '6', '9', '8', '7', '9', '8', '8', '9', '8', '9',
  '9', '9', '0', '9', '9', '1', '9', '9', '2', '9', '9', '3', '9', '9', '4',
  '9', '9', '5', '9', '9', '6', '9', '9', '7', '9', '9', '8', '9', '9', '9',
};

static uint32_t convert_num[128] = {
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 2, 3, 4, 5, 6, 7,
  8, 9, 0, 0, 0, 0, 0, 0,
  0, 10, 11, 12, 13, 14, 15, 16,
  17, 18, 19, 20, 21, 22, 23, 24,
  25, 26, 27, 28, 29, 30, 31, 32,
  33, 34, 35, 0, 0, 0, 0, 0,
  0, 10, 11, 12, 13, 14, 15, 16,
  17, 18, 19, 20, 21, 22, 23, 24,
  25, 26, 27, 28, 29, 30, 31, 32,
  33, 34, 35, 0, 0, 0, 0, 0,
};

static unsigned char uppercase[256] = {
  0, 1, 2, 3, 4, 5, 6, 7,
  8, 9, 10, 11, 12, 13, 14, 15,
  16, 17, 18, 19, 20, 21, 22, 23,
  24, 25, 26, 27, 28, 29, 30, 31,
  32, 33, 34, 35, 36, 37, 38, 39,
  40, 41, 42, 43, 44, 45, 46, 47,
  48, 49, 50, 51, 52, 53, 54, 55,
  56, 57, 58, 59, 60, 61, 62, 63,
  64, 65, 66, 67, 68, 69, 70, 71,
  72, 73, 74, 75, 76, 77, 78, 79,
  80, 81, 82, 83, 84, 85, 86, 87,
  88, 89, 90, 91, 92, 93, 94, 95,
  96, 65, 66, 67, 68, 69, 70, 71,
  72, 73, 74, 75, 76, 77, 78, 79,
  80, 81, 82, 83, 84, 85, 86, 87,
  88, 89, 90, 123, 124, 125, 126, 127,
  128, 129, 130, 131, 132, 133, 134, 135,
  136, 137, 138, 139, 140, 141, 142, 143,
  144, 145, 146, 147, 148, 149, 150, 151,
  152, 153, 154, 155, 156, 157, 158, 159,
  160, 161, 162, 163, 164, 165, 166, 167,
  168, 169, 170, 171, 172, 173, 174, 175,
  176, 177, 178, 179, 180, 181, 182, 183,
  184, 185, 186, 187, 188, 189, 190, 191,
  192, 193, 194, 195, 196, 197, 198, 199,
  200, 201, 202, 203, 204, 205, 206, 207,
  208, 209, 210, 211, 212, 213, 214, 215,
  216, 217, 218, 219, 220, 221, 222, 223,
  224, 225, 226, 227, 228, 229, 230, 231,
  232, 233, 234, 235, 236, 237, 238, 239,
  240, 241, 242, 243, 244, 245, 246, 247,
  248, 249, 250, 251, 252, 253, 254, 255
};

static uint32_t valid_num[256] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
  0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

#define	DECIMAL_POINT	(*lconv->decimal_point)

static struct lconv default_lconv = {.decimal_point = "." };

static struct lconv *lconv = &default_lconv;

/** \brief init_fast_convert
 *
 * \b Description
 *
 * Update decimal_point from locale at startup
 */
static void __attribute__((constructor))
  init_fast_convert (void)
{
  lconv = localeconv ();
}

#ifndef __GNUC__

static const unsigned char fast_clz[256] = {
  0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8
};

/** \brief calc_clz32
 *
 * \b Description
 *
 * Calulate leading zero's
 *
 * \param x 32 bits input
 * \returns Number of leading zero's
 */

static unsigned int
calc_clz32 (uint32_t x)
{
  unsigned int n;

  return (n = (x >> 24)) ? 8 - fast_clz[n] :
    (n = (x >> 16)) ? 16 - fast_clz[n] :
    (n = (x >> 8)) ? 24 - fast_clz[n] :
    32 - fast_clz[x];
}

/** \brief calc_clz64
 *
 * \b Description
 *
 * Calulate leading zero's
 *
 * \param x 64 bits input
 * \returns Number of leading zero's
 */

static unsigned int
calc_clz64 (uint64_t x)
{
  unsigned int n, t;

  return ((t = (x >> 32))) ?
    ((n = (t >> 24)) ? 8 - fast_clz[n] :
     (n = (t >> 16)) ? 16 - fast_clz[n] :
     (n = (t >> 8)) ? 24 - fast_clz[n] :
     32 - fast_clz[t]) :
    ((n = (x >> 24)) ? 40 - fast_clz[n] :
     (n = (x >> 16)) ? 48 - fast_clz[n] :
     (n = (x >> 8)) ? 56 - fast_clz[n] :
     64 - fast_clz[x]);
}
#endif

static unsigned int
log10_32 (uint32_t v)
{
  static const unsigned int tab32[32] = {
    1, 1, 1, 1, 1, 1, 2, 2, 2, 3,
    3, 3, 3, 4, 4, 4, 5, 5, 5, 6,
    6, 6, 6, 7, 7, 7, 8, 8, 8, 9,
    9, 9
  };
  unsigned int n;

  if (v == 0) {
    return 1;
  }
#ifdef __GNUC__
  n = tab32[31 ^ __builtin_clz (v)];
#else
  n = tab32[31 ^ calc_clz32 (v)];
#endif
  return n + (v >= ipowers32[n]);
}

static unsigned int
log10_64 (uint64_t v)
{
  static const unsigned int tab64[64] = {
    1, 1, 1, 1, 1, 1, 2, 2, 2, 3,
    3, 3, 3, 4, 4, 4, 5, 5, 5, 6,
    6, 6, 6, 7, 7, 7, 8, 8, 8, 9,
    9, 9, 10, 10, 10, 11, 11, 11, 12, 12,
    12, 13, 13, 13, 14, 14, 14, 15, 15, 15,
    16, 16, 16, 16, 17, 17, 17, 18, 18, 18,
    18, 19, 19, 19
  };
  unsigned int n;

  if (v == 0) {
    return 1;
  }
#ifdef WIN
  n = tab64[63 ^ __builtin_clzll (v)];
#else
#ifdef __GNUC__
#if __WORDSIZE == 64
  n = tab64[63 ^ __builtin_clzl (v)];
#else
  n = tab64[63 ^ __builtin_clzll (v)];
#endif
#else
  n = tab64[63 ^ calc_clz64 (v)];
#endif
#endif
  return n + (v >= ipowers64[n]);
}

#if __WORDSIZE != 64
static void
do_uint (uint32_t v, char *p, unsigned int len)
{
  p += len;
  while (v >= 1000) {
    uint32_t d = v / 1000;

    p = (char *) memcpy (p - 3, &num3[(v - d * 1000) * 3], 3);
    v = d;
  }
  if (v >= 100) {
    memcpy (p - 3, &num3[v * 3], 3);
  }
  else if (v >= 10) {
    memcpy (p - 2, &num3[v * 3] + 1, 2);
  }
  else {
    memcpy (p - 1, &num3[v * 3] + 2, 1);
  }
}

static const uint64_t
div_10 (const uint64_t n)
{
  /* n * 14757395258967641293 >> 64 >> 3 */
  uint32_t n_lo;
  uint32_t n_hi;
  uint32_t c_lo;
  uint32_t c_hi;
  uint64_t res;
  uint64_t tmp;

  n_lo = n;
  n_hi = n >> 32;
  c_lo = 3435973837u;
  c_hi = 3435973836u;
  res = ((uint64_t) n_lo * c_lo) >> 32;
  tmp = res += (uint64_t) n_lo *c_hi;
  res += (uint64_t) n_hi *c_lo;
  res = (res >> 32) + (res < tmp ? UINT64_C (4294967296) : 0);
  res += (uint64_t) n_hi *c_hi;
  return res >> 3;
}

static const uint64_t
div_100 (const uint64_t n)
{
  /* (n >> 2) * 2951479051793528259 >> 64 >> 2 */
  uint32_t n_lo;
  uint32_t n_hi;
  uint32_t c_lo;
  uint32_t c_hi;
  uint64_t res;
  uint64_t tmp;

  tmp = n >> 2;
  n_lo = tmp;
  n_hi = tmp >> 32;
  c_lo = 1546188227u;
  c_hi = 687194767u;
  res = ((uint64_t) n_lo * c_lo) >> 32;
  tmp = res += (uint64_t) n_lo *c_hi;
  res += (uint64_t) n_hi *c_lo;
  res = (res >> 32) + (res < tmp ? UINT64_C (4294967296) : 0);
  res += (uint64_t) n_hi *c_hi;
  return res >> 2;
}

static const uint64_t
div_1000 (const uint64_t n)
{
  /* (n >> 3) * 2361183241434822607 >> 64 >> 4 */
  uint32_t n_lo;
  uint32_t n_hi;
  uint32_t c_lo;
  uint32_t c_hi;
  uint64_t res;
  uint64_t tmp;

  tmp = n >> 3;
  n_lo = tmp;
  n_hi = tmp >> 32;
  c_lo = 3813930959u;
  c_hi = 549755813u;
  res = ((uint64_t) n_lo * c_lo) >> 32;
  tmp = res += (uint64_t) n_lo *c_hi;
  res += (uint64_t) n_hi *c_lo;
  res = (res >> 32) + (res < tmp ? UINT64_C (4294967296) : 0);
  res += (uint64_t) n_hi *c_hi;
  return res >> 4;
}

static const uint64_t
div_1000000000 (const uint64_t n)
{
  /* (n >> 9) * 19342813113834067 >> 64 >> 11 */
  uint32_t n_lo;
  uint32_t n_hi;
  uint32_t c_lo;
  uint32_t c_hi;
  uint64_t res;
  uint64_t tmp;

  tmp = n >> 9;
  n_lo = tmp;
  n_hi = tmp >> 32;
  c_lo = 2694535763u;
  c_hi = 4503599u;
  res = ((uint64_t) n_lo * c_lo) >> 32;
  tmp = res += (uint64_t) n_lo *c_hi;
  res += (uint64_t) n_hi *c_lo;
  res = (res >> 32) + (res < tmp ? UINT64_C (4294967296) : 0);
  res += (uint64_t) n_hi *c_hi;
  return res >> 11;
}
#endif

/** \brief mul_10_add
 *
 * \b Description
 *
 * Multiply 128 bit number * 10 and add a
 *
 * \param m pointer to 64 bit msb
 * \param l pointer to 64 bit lsb
 * \param a value to add
 */

static void
mul_10_add (uint64_t * m, uint64_t * l, uint32_t a)
{
  uint64_t hi;
  uint64_t lo;

  *m = (*m << 1) + (*l >> 63);
  *l <<= 1;
  hi = (*m << 2) + (*l >> 62);
  lo = *l << 2;
  *l += lo;
  *m += hi + (*l < lo);
  *l += a;
  *m += (*l < a);
}

/** \brief mul_96
 *
 * \b Description
 *
 * Calulate n * m >> 96
 *
 * \param n Long integer
 * \param m Buffer to print to
 * \returns n * m >> 96
 */

#if defined (__GNUC__) && defined(__x86_64__)
static const uint64_t
mul_96 (const uint64_t n, const uint64_t m, const uint32_t o, uint32_t * low)
{
  uint64_t hi1, lo1, hi2, lo2;

__asm__ ("mulq %3": "=d" (hi1), "=a" (lo1): "%a" (n), "rm" (m):"cc");
__asm__ ("mulq %3": "=d" (hi2), "=a" (lo2): "%a" (n), "rm" ((uint64_t) o):"cc");
  lo2 >>= 32;
  lo1 += lo2;
  lo1 = (lo1 >> 32) + ((lo1 < lo2) ? UINT64_C (4294967296) : 0);
  *low = hi2 += lo1;
  return hi1 + (hi2 >> 32) + ((hi2 < lo1) ? UINT64_C (4294967296) : 0);
}
#else
static const uint64_t
mul_96 (const uint64_t n, const uint64_t m, const uint32_t o, uint32_t * low)
{
  uint32_t n_lo;
  uint32_t n_hi;
  uint32_t c_lo;
  uint32_t c_hi;
  uint64_t res;
  uint64_t tmp;

  n_lo = n;
  n_hi = n >> 32;
  c_lo = m;
  c_hi = m >> 32;
  res = ((uint64_t) n_lo * o) >> 32;
  tmp = res += (uint64_t) n_hi *o;
  res += (uint64_t) n_lo *c_lo;
  res = (res >> 32) + ((res < tmp) ? UINT64_C (4294967296) : 0);
  tmp = res += (uint64_t) n_lo *c_hi;
  *low = res += (uint64_t) n_hi *c_lo;
  res = (res >> 32) + ((res < tmp) ? UINT64_C (4294967296) : 0);
  res += (uint64_t) n_hi *c_hi;
  return res;
}
#endif

/** \brief mul_64
 *
 * \b Description
 *
 * Calulate n * m >> 64
 *
 * \param n Long integer
 * \param m Buffer to print to
 * \returns n * m >> 64
 */

#if defined (__GNUC__) && defined(__x86_64__)
static const uint64_t
mul_64 (const uint64_t n, const uint64_t m)
{
  uint64_t hi, lo;
__asm__ ("mulq %3": "=d" (hi), "=a" (lo): "%a" (n), "rm" (m):"cc");
  return hi;
}
#else
static const uint64_t
mul_64 (const uint64_t n, const uint64_t m)
{
  uint32_t n_lo;
  uint32_t n_hi;
  uint32_t c_lo;
  uint32_t c_hi;
  uint64_t res;
  uint64_t tmp;

  n_lo = n;
  n_hi = n >> 32;
  c_lo = m;
  c_hi = m >> 32;
  res = ((uint64_t) n_lo * c_lo) >> 32;
  tmp = res += (uint64_t) n_lo *c_hi;
  res += (uint64_t) n_hi *c_lo;
  res = (res >> 32) + ((res < tmp) ? UINT64_C (4294967296) : 0);
  res += (uint64_t) n_hi *c_hi;
  return res;
}
#endif

/** \brief mul_56
 *
 * \b Description
 *
 * Calulate n * m >> 56
 *
 * \param n Long integer
 * \param m Buffer to print to
 * \returns n * m >> 56
 */

#if defined (__GNUC__) && defined(__x86_64__)
static const uint64_t
mul_56 (const uint32_t n, const uint64_t m, uint32_t * low)
{
  uint64_t hi, lo;

__asm__ ("mulq %3": "=d" (hi), "=a" (lo): "%a" (n), "rm" (m):"cc");
  *low = lo >> 24;
  return (hi << 8) + (lo >> 56);
}
#else
static const uint64_t
mul_56 (const uint32_t n, const uint64_t m, uint32_t * low)
{
  uint64_t res;

  res = n * (m & UINT64_C (0xFFFFFFFF));
  *low = (res >> 24) & 0xFFu;
  res >>= 32;
  res += n * (m >> 32);
  *low |= res << 8;
  return res >> 24;
}
#endif

/** \brief fast_sint32
 *
 * \b Description
 *
 * Convert signed integer to string
 *
 * \param v Integer
 * \param str Buffer to print to
 * \returns lenght string
 */

unsigned int
fast_sint32 (int32_t v, char *str)
{
  unsigned int is_signed = v < 0;
  unsigned int len;
  char *p;
  uint32_t j = (uint32_t) v;

  if (is_signed) {
    j = ~j + 1;
    *str++ = '-';
  }
  len = log10_32 (j);
  p = str + len;
  *p = '\0';
  while (j >= 1000) {
    uint32_t d = j / 1000;

    p = (char *) memcpy (p - 3, &num3[(j - d * 1000) * 3], 3);
    j = d;
  }
  if (j >= 100) {
    memcpy (p - 3, &num3[j * 3], 3);
  }
  else if (j >= 10) {
    memcpy (p - 2, &num3[j * 3] + 1, 2);
  }
  else {
    memcpy (p - 1, &num3[j * 3] + 2, 1);
  }
  return len + is_signed;
}

/** \brief fast_sint64
 *
 * \b Description
 *
 * Convert signed long integer to string
 *
 * \param v Long integer
 * \param str Buffer to print to
 * \returns lenght string
 */

unsigned int
fast_sint64 (int64_t v, char *str)
{
  unsigned int is_signed = v < 0;
  unsigned int len;
  char *p;
  uint64_t j = (uint64_t) v;

  if (is_signed) {
    j = ~j + 1;
    *str++ = '-';
  }
  len = log10_64 (j);
  p = str + len;
  *p = '\0';
#if __WORDSIZE == 64
  while (j >= 1000) {
    uint64_t d = j / 1000;

    p = (char *) memcpy (p - 3, &num3[(j - d * 1000) * 3], 3);
    j = d;
  }
  if (j >= 100) {
    memcpy (p - 3, &num3[j * 3], 3);
  }
  else if (j >= 10) {
    memcpy (p - 2, &num3[j * 3] + 1, 2);
  }
  else {
    memcpy (p - 1, &num3[j * 3] + 2, 1);
  }
#else
  if (j >= 1000000000) {
    uint64_t d = div_1000000000 (j);
    uint32_t m = j - d * 1000000000;

    do_uint (m + 1000000000, p - 10, 10);
    j = d;
    p -= 9;
    if (j >= 1000000000) {
      d = div_1000000000 (j);
      m = j - d * 1000000000;

      do_uint (m + 1000000000, p - 10, 10);
      j = d;
      p -= 9;
    }
  }
  do_uint (j, str, p - str);
#endif
  return len + is_signed;
}

/** \brief fast_uint32
 *
 * \b Description
 *
 * Convert unsigned integer to string
 *
 * \param v Integer
 * \param str Buffer to print to
 * \returns lenght string
 */

unsigned int
fast_uint32 (uint32_t v, char *str)
{
  unsigned int len = log10_32 (v);
  char *p = str + len;

  *p = '\0';
  while (v >= 1000) {
    uint32_t d = v / 1000;

    p = (char *) memcpy (p - 3, &num3[(v - d * 1000) * 3], 3);
    v = d;
  }
  if (v >= 100) {
    memcpy (p - 3, &num3[v * 3], 3);
  }
  else if (v >= 10) {
    memcpy (p - 2, &num3[v * 3] + 1, 2);
  }
  else {
    memcpy (p - 1, &num3[v * 3] + 2, 1);
  }
  return len;
}

/** \brief fast_uint64
 *
 * \b Description
 *
 * Convert unsigned long integer to string
 *
 * \param v Long integer
 * \param str Buffer to print to
 * \returns lenght string
 */

unsigned int
fast_uint64 (uint64_t v, char *str)
{
  unsigned int len = log10_64 (v);
  char *p = str + len;

  *p = '\0';
#if __WORDSIZE == 64
  while (v >= 1000) {
    uint64_t d = v / 1000;

    p = (char *) memcpy (p - 3, &num3[(v - d * 1000) * 3], 3);
    v = d;
  }
  if (v >= 100) {
    memcpy (p - 3, &num3[v * 3], 3);
  }
  else if (v >= 10) {
    memcpy (p - 2, &num3[v * 3] + 1, 2);
  }
  else {
    memcpy (p - 1, &num3[v * 3] + 2, 1);
  }
#else
  if (v >= 1000000000) {
    uint64_t d = div_1000000000 (v);
    uint32_t m = v - d * 1000000000;

    do_uint (m + 1000000000, p - 10, 10);
    v = d;
    p -= 9;
    if (v >= 1000000000) {
      d = div_1000000000 (v);
      m = v - d * 1000000000;

      do_uint (m + 1000000000, p - 10, 10);
      v = d;
      p -= 9;
    }
  }
  do_uint (v, str, p - str);
#endif
  return len;
}

/** \brief fast_base_sint32
 *
 * \b Description
 *
 * Convert signed integer to string with base and upper case
 *
 * \param v Integer
 * \param str Buffer to print to
 * \param base Base to use
 * \param upper Use uppecase
 * \returns lenght string
 */

unsigned int
fast_base_sint32 (int32_t s, char *str, int base, int upper)
{
  const char *d = upper ? "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    : "0123456789abcdefghijklmnopqrstuvwxyz";
  char tmp[100];
  char *p = &tmp[sizeof (tmp) - 1];
  unsigned int is_signed = s < 0;
  uint32_t u = (uint32_t) s;

  *p = '\0';
  if (is_signed) {
    u = ~u + 1;
  }
  switch (base) {
    DO_BASE (2);
    DO_BASE (3);
    DO_BASE (4);
    DO_BASE (5);
    DO_BASE (6);
    DO_BASE (7);
    DO_BASE (8);
    DO_BASE (9);
  case 10:
    return fast_sint32 (s, str);
    DO_BASE (11);
    DO_BASE (12);
    DO_BASE (13);
    DO_BASE (14);
    DO_BASE (15);
    DO_BASE (16);
    DO_BASE (17);
    DO_BASE (18);
    DO_BASE (19);
    DO_BASE (20);
    DO_BASE (21);
    DO_BASE (22);
    DO_BASE (23);
    DO_BASE (24);
    DO_BASE (25);
    DO_BASE (26);
    DO_BASE (27);
    DO_BASE (28);
    DO_BASE (29);
    DO_BASE (30);
    DO_BASE (31);
    DO_BASE (32);
    DO_BASE (33);
    DO_BASE (34);
    DO_BASE (35);
    DO_BASE (36);
  }
  if (is_signed) {
    *--p = '-';
  }
  memcpy (str, p, &tmp[sizeof (tmp)] - p);
  return &tmp[sizeof (tmp) - 1] - p;
}

/** \brief fast_base_sint64
 *
 * \b Description
 *
 * Convert signed long integer to string with base and upper case
 *
 * \param v Long integer
 * \param str Buffer to print to
 * \param base Base to use
 * \param upper Use uppecase
 * \returns lenght string
 */
unsigned int
fast_base_sint64 (int64_t s, char *str, int base, int upper)
{
  const char *d = upper ? "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    : "0123456789abcdefghijklmnopqrstuvwxyz";
  char tmp[100];
  char *p = &tmp[sizeof (tmp) - 1];
  unsigned int is_signed = s < 0;
  uint64_t u = (uint64_t) s;

  *p = '\0';
  if (is_signed) {
    u = ~u + 1;
  }
  switch (base) {
    DO_BASE (2);
    DO_BASE (3);
    DO_BASE (4);
    DO_BASE (5);
    DO_BASE (6);
    DO_BASE (7);
    DO_BASE (8);
    DO_BASE (9);
  case 10:
    return fast_sint64 (s, str);
    DO_BASE (11);
    DO_BASE (12);
    DO_BASE (13);
    DO_BASE (14);
    DO_BASE (15);
    DO_BASE (16);
    DO_BASE (17);
    DO_BASE (18);
    DO_BASE (19);
    DO_BASE (20);
    DO_BASE (21);
    DO_BASE (22);
    DO_BASE (23);
    DO_BASE (24);
    DO_BASE (25);
    DO_BASE (26);
    DO_BASE (27);
    DO_BASE (28);
    DO_BASE (29);
    DO_BASE (30);
    DO_BASE (31);
    DO_BASE (32);
    DO_BASE (33);
    DO_BASE (34);
    DO_BASE (35);
    DO_BASE (36);
  }
  if (is_signed) {
    *--p = '-';
  }
  memcpy (str, p, &tmp[sizeof (tmp)] - p);
  return &tmp[sizeof (tmp) - 1] - p;
}

/** \brief fast_uint32
 *
 * \b Description
 *
 * Convert unsigned base_integer to string with base and upper case
 *
 * \param v Integer
 * \param str Buffer to print to
 * \param base Base to use
 * \param upper Use uppecase
 * \returns lenght string
 */

unsigned int
fast_base_uint32 (uint32_t u, char *str, int base, int upper)
{
  const char *d = upper ? "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    : "0123456789abcdefghijklmnopqrstuvwxyz";
  char tmp[100];
  char *p = &tmp[sizeof (tmp) - 1];

  *p = '\0';
  switch (base) {
    DO_BASE (2);
    DO_BASE (3);
    DO_BASE (4);
    DO_BASE (5);
    DO_BASE (6);
    DO_BASE (7);
    DO_BASE (8);
    DO_BASE (9);
  case 10:
    return fast_uint32 (u, str);
    DO_BASE (11);
    DO_BASE (12);
    DO_BASE (13);
    DO_BASE (14);
    DO_BASE (15);
    DO_BASE (16);
    DO_BASE (17);
    DO_BASE (18);
    DO_BASE (19);
    DO_BASE (20);
    DO_BASE (21);
    DO_BASE (22);
    DO_BASE (23);
    DO_BASE (24);
    DO_BASE (25);
    DO_BASE (26);
    DO_BASE (27);
    DO_BASE (28);
    DO_BASE (29);
    DO_BASE (30);
    DO_BASE (31);
    DO_BASE (32);
    DO_BASE (33);
    DO_BASE (34);
    DO_BASE (35);
    DO_BASE (36);
  }
  memcpy (str, p, &tmp[sizeof (tmp)] - p);
  return &tmp[sizeof (tmp) - 1] - p;
}

/** \brief fast_base_uint64
 *
 * \b Description
 *
 * Convert unsigned long integer to string with base and upper case
 *
 * \param v Long integer
 * \param str Buffer to print to
 * \param base Base to use
 * \param upper Use uppecase
 * \returns lenght string
 */

unsigned int
fast_base_uint64 (uint64_t u, char *str, int base, int upper)
{
  const char *d = upper ? "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    : "0123456789abcdefghijklmnopqrstuvwxyz";
  char tmp[100];
  char *p = &tmp[sizeof (tmp) - 1];

  *p = '\0';
  switch (base) {
    DO_BASE (2);
    DO_BASE (3);
    DO_BASE (4);
    DO_BASE (5);
    DO_BASE (6);
    DO_BASE (7);
    DO_BASE (8);
    DO_BASE (9);
  case 10:
    return fast_uint64 (u, str);
    DO_BASE (11);
    DO_BASE (12);
    DO_BASE (13);
    DO_BASE (14);
    DO_BASE (15);
    DO_BASE (16);
    DO_BASE (17);
    DO_BASE (18);
    DO_BASE (19);
    DO_BASE (20);
    DO_BASE (21);
    DO_BASE (22);
    DO_BASE (23);
    DO_BASE (24);
    DO_BASE (25);
    DO_BASE (26);
    DO_BASE (27);
    DO_BASE (28);
    DO_BASE (29);
    DO_BASE (30);
    DO_BASE (31);
    DO_BASE (32);
    DO_BASE (33);
    DO_BASE (34);
    DO_BASE (35);
    DO_BASE (36);
  }
  memcpy (str, p, &tmp[sizeof (tmp)] - p);
  return &tmp[sizeof (tmp) - 1] - p;
}

/** \brief fast_strtos32
 *
 * \b Description
 *
 * Convert string to signed integer
 *
 * \param str String to convert from
 * \param endptr optional endptr
 * \returns converted string
 */

int32_t
fast_strtos32 (const char *str, char **endptr, int base)
{
  uint32_t n;
  unsigned int sign = 0;
  unsigned char *cp = (unsigned char *) str;
#define	M(x)	(2147483647u / (x))
  static const uint32_t maxp[37] = {
    0, 0, M (2), M (3), M (4), M (5), M (6), M (7), M (8), M (9), M (10),
    M (11), M (12), M (13), M (14), M (15), M (16), M (17), M (18),
    M (19), M (20), M (21), M (22), M (23), M (24), M (25), M (26),
    M (27), M (28), M (29), M (30), M (31), M (32), M (33), M (34),
    M (35), M (36)
  };
#undef M
#define	M(x)	(2147483648u / (x))
  static const uint32_t maxn[37] = {
    0, 0, M (2), M (3), M (4), M (5), M (6), M (7), M (8), M (9), M (10),
    M (11), M (12), M (13), M (14), M (15), M (16), M (17), M (18),
    M (19), M (20), M (21), M (22), M (23), M (24), M (25), M (26),
    M (27), M (28), M (29), M (30), M (31), M (32), M (33), M (34),
    M (35), M (36)
  };
#undef M
#define	R(x)	(2147483647u % (x))
  static const uint32_t remp[37] = {
    0, 0, R (2), R (3), R (4), R (5), R (6), R (7), R (8), R (9), R (10),
    R (11), R (12), R (13), R (14), R (15), R (16), R (17), R (18),
    R (19), R (20), R (21), R (22), R (23), R (24), R (25), R (26),
    R (27), R (28), R (29), R (30), R (31), R (32), R (33), R (34),
    R (35), R (36)
  };
#undef R
#define	R(x)	(2147483648u % (x))
  static const uint32_t remn[37] = {
    0, 0, R (2), R (3), R (4), R (5), R (6), R (7), R (8), R (9), R (10),
    R (11), R (12), R (13), R (14), R (15), R (16), R (17), R (18),
    R (19), R (20), R (21), R (22), R (23), R (24), R (25), R (26),
    R (27), R (28), R (29), R (30), R (31), R (32), R (33), R (34),
    R (35), R (36)
  };
#undef R

  while (isspace (*cp)) {
    cp++;
  }
  if (*cp == '+') {
    cp++;
  }
  else if (*cp == '-') {
    sign = 1;
    cp++;
  }
  n = 0;
  if (UNLIKELY (base) &&
      (base != 10 || *cp < '1') &&
      (base != 8 || *cp != '0' || cp[1] < '0' || cp[1] > '7') &&
      (base != 16 || *cp != '0' || (cp[1] != 'x' && cp[1] != 'X'))) {
    if (base >= 2 && base <= 36) {
      uint32_t max = sign ? maxn[base] : maxp[base];
      uint32_t rem = sign ? remn[base] : remp[base];

      if (base <= 10) {
	unsigned char u = *cp;
	unsigned char end = '0' + base;

	if (u >= '0' && u < end) {
	  do {
	    uint32_t v = convert_num[u];

	    if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	      break;
	    }
	    n = n * base + v;
	    u = *++cp;
	  } while (u >= '0' && u < end);
	}
	else {
	  cp = (unsigned char *) str;
	}
      }
      else {
	unsigned char u = uppercase[*cp];
	unsigned char end = 'A' + base - 10;

	if (valid_num[u] && u < end) {
	  do {
	    uint32_t v = convert_num[u];

	    if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	      break;
	    }
	    n = n * base + v;
	    u = uppercase[*++cp];
	  } while (valid_num[u] && u < end);
	}
	else {
	  cp = (unsigned char *) str;
	}
      }
    }
    else {
      cp = (unsigned char *) str;
    }
  }
  else if (*cp == '0') {
    cp++;
    if ((*cp == 'x' || *cp == 'X') && isxdigit (cp[1])) {
      uint32_t max = sign ? maxn[16] : maxp[16];
      uint32_t rem = sign ? remn[16] : remp[16];
      unsigned char u = *++cp;

      do {
	uint32_t v = convert_num[u];

	if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	  break;
	}
	n = n * 16 + v;
	u = *++cp;
      } while (isxdigit (u));
    }
    else {
      uint32_t max = sign ? maxn[8] : maxp[8];
      uint32_t rem = sign ? remn[8] : remp[8];
      unsigned char u = *cp;

      while (u >= '0' && u <= '7') {
	uint32_t v = convert_num[u];

	if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	  break;
	}
	n = n * 8 + v;
	u = *++cp;
      }
    }
  }
  else if (isdigit (*cp)) {
    uint32_t max = sign ? maxn[10] : maxp[10];
    uint32_t rem = sign ? remn[10] : remp[10];
    unsigned char u = *cp;

    do {
      uint32_t v = convert_num[u];

      if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	break;
      }
      n = n * 10 + v;
      u = *++cp;
    } while (isdigit (u));
  }
  else {
    cp = (unsigned char *) str;
  }
  if (endptr) {
    *endptr = (char *) cp;
  }
  return sign ? -(int32_t) n : (int32_t) n;
}

/** \brief fast_strtos64
 *
 * \b Description
 *
 * Convert string to signed long
 *
 * \param str String to convert from
 * \param endptr optional endptr
 * \returns converted string
 */

int64_t
fast_strtos64 (const char *str, char **endptr, int base)
{
  uint64_t n;
  unsigned int sign = 0;
  unsigned char *cp = (unsigned char *) str;
#define	M(x)	(UINT64_C(9223372036854775807) / (x))
  static const uint64_t maxp[37] = {
    0, 0, M (2), M (3), M (4), M (5), M (6), M (7), M (8), M (9), M (10),
    M (11), M (12), M (13), M (14), M (15), M (16), M (17), M (18),
    M (19), M (20), M (21), M (22), M (23), M (24), M (25), M (26),
    M (27), M (28), M (29), M (30), M (31), M (32), M (33), M (34),
    M (35), M (36)
  };
#undef M
#define	M(x)	(UINT64_C(9223372036854775808) / (x))
  static const uint64_t maxn[37] = {
    0, 0, M (2), M (3), M (4), M (5), M (6), M (7), M (8), M (9), M (10),
    M (11), M (12), M (13), M (14), M (15), M (16), M (17), M (18),
    M (19), M (20), M (21), M (22), M (23), M (24), M (25), M (26),
    M (27), M (28), M (29), M (30), M (31), M (32), M (33), M (34),
    M (35), M (36)
  };
#undef M
#define	R(x)	(UINT64_C(9223372036854775807) % (x))
  static const uint64_t remp[37] = {
    0, 0, R (2), R (3), R (4), R (5), R (6), R (7), R (8), R (9), R (10),
    R (11), R (12), R (13), R (14), R (15), R (16), R (17), R (18),
    R (19), R (20), R (21), R (22), R (23), R (24), R (25), R (26),
    R (27), R (28), R (29), R (30), R (31), R (32), R (33), R (34),
    R (35), R (36)
  };
#undef R
#define	R(x)	(UINT64_C(9223372036854775808) % (x))
  static const uint64_t remn[37] = {
    0, 0, R (2), R (3), R (4), R (5), R (6), R (7), R (8), R (9), R (10),
    R (11), R (12), R (13), R (14), R (15), R (16), R (17), R (18),
    R (19), R (20), R (21), R (22), R (23), R (24), R (25), R (26),
    R (27), R (28), R (29), R (30), R (31), R (32), R (33), R (34),
    R (35), R (36)
  };
#undef R

  while (isspace (*cp)) {
    cp++;
  }
  if (*cp == '+') {
    cp++;
  }
  else if (*cp == '-') {
    sign = 1;
    cp++;
  }
  n = 0;
  if (UNLIKELY (base) &&
      (base != 10 || *cp < '1') &&
      (base != 8 || *cp != '0' || cp[1] < '0' || cp[1] > '7') &&
      (base != 16 || *cp != '0' || (cp[1] != 'x' && cp[1] != 'X'))) {
    if (base >= 2 && base <= 36) {
      uint64_t max = sign ? maxn[base] : maxp[base];
      uint64_t rem = sign ? remn[base] : remp[base];

      if (base <= 10) {
	unsigned char u = *cp;
	unsigned char end = '0' + base;

	if (u >= '0' && u < end) {
	  do {
	    uint32_t v = convert_num[u];

	    if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	      break;
	    }
	    n = n * base + v;
	    u = *++cp;
	  } while (u >= '0' && u < end);
	}
	else {
	  cp = (unsigned char *) str;
	}
      }
      else {
	unsigned char u = uppercase[*cp];
	unsigned char end = 'A' + base - 10;

	if (valid_num[u] && u < end) {
	  do {
	    uint32_t v = convert_num[u];

	    if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	      break;
	    }
	    n = n * base + v;
	    u = uppercase[*++cp];
	  } while (valid_num[u] && u < end);
	}
	else {
	  cp = (unsigned char *) str;
	}
      }
    }
    else {
      cp = (unsigned char *) str;
    }
  }
  else if (*cp == '0') {
    cp++;
    if ((*cp == 'x' || *cp == 'X') && isxdigit (cp[1])) {
      uint64_t max = sign ? maxn[16] : maxp[16];
      uint64_t rem = sign ? remn[16] : remp[16];
      unsigned char u = *++cp;

      do {
	uint32_t v = convert_num[u];

	if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	  break;
	}
	n = n * 16 + v;
	u = *++cp;
      } while (isxdigit (u));
    }
    else {
      uint64_t max = sign ? maxn[8] : maxp[8];
      uint64_t rem = sign ? remn[8] : remp[8];
      unsigned char u = *cp;

      while (u >= '0' && u <= '7') {
	uint32_t v = convert_num[u];

	if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	  break;
	}
	n = n * 8 + v;
	u = *++cp;
      }
    }
  }
  else if (isdigit (*cp)) {
    uint64_t max = sign ? maxn[10] : maxp[10];
    uint64_t rem = sign ? remn[10] : remp[10];
    unsigned char u = *cp;

    do {
      uint32_t v = convert_num[u];

      if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	break;
      }
      n = n * 10 + v;
      u = *++cp;
    } while (isdigit (u));
  }
  else {
    cp = (unsigned char *) str;
  }
  if (endptr) {
    *endptr = (char *) cp;
  }
  return sign ? -(int64_t) n : (int64_t) n;
}

/** \brief fast_strtou32
 *
 * \b Description
 *
 * Convert string to unsigned int
 *
 * \param str String to convert from
 * \param endptr optional endptr
 * \returns converted string
 */

uint32_t
fast_strtou32 (const char *str, char **endptr, int base)
{
  uint32_t n;
  unsigned char *cp = (unsigned char *) str;
#define	M(x)	(4294967295u / (x))
  static const uint32_t maxp[37] = {
    0, 0, M (2), M (3), M (4), M (5), M (6), M (7), M (8), M (9), M (10),
    M (11), M (12), M (13), M (14), M (15), M (16), M (17), M (18),
    M (19), M (20), M (21), M (22), M (23), M (24), M (25), M (26),
    M (27), M (28), M (29), M (30), M (31), M (32), M (33), M (34),
    M (35), M (36)
  };
#undef M
#define	R(x)	(4294967295u % (x))
  static const uint32_t remp[37] = {
    0, 0, R (2), R (3), R (4), R (5), R (6), R (7), R (8), R (9), R (10),
    R (11), R (12), R (13), R (14), R (15), R (16), R (17), R (18),
    R (19), R (20), R (21), R (22), R (23), R (24), R (25), R (26),
    R (27), R (28), R (29), R (30), R (31), R (32), R (33), R (34),
    R (35), R (36)
  };
#undef R

  while (isspace (*cp)) {
    cp++;
  }
  n = 0;
  if (UNLIKELY (base) &&
      (base != 10 || *cp < '1') &&
      (base != 8 || *cp != '0' || cp[1] < '0' || cp[1] > '7') &&
      (base != 16 || *cp != '0' || (cp[1] != 'x' && cp[1] != 'X'))) {
    if (base >= 2 && base <= 36) {
      uint32_t max = maxp[base];
      uint32_t rem = remp[base];

      if (base <= 10) {
	unsigned char u = *cp;
	unsigned char end = '0' + base;

	if (u >= '0' && u < end) {
	  do {
	    uint32_t v = convert_num[u];

	    if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	      break;
	    }
	    n = n * base + v;
	    u = *++cp;
	  } while (u >= '0' && u < end);
	}
	else {
	  cp = (unsigned char *) str;
	}
      }
      else {
	unsigned char u = uppercase[*cp];
	unsigned char end = 'A' + base - 10;

	if (valid_num[u] && u < end) {
	  do {
	    uint32_t v = convert_num[u];

	    if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	      break;
	    }
	    n = n * base + v;
	    u = uppercase[*++cp];
	  } while (valid_num[u] && u < end);
	}
	else {
	  cp = (unsigned char *) str;
	}
      }
    }
    else {
      cp = (unsigned char *) str;
    }
  }
  else if (*cp == '0') {
    cp++;
    if ((*cp == 'x' || *cp == 'X') && isxdigit (cp[1])) {
      uint32_t max = maxp[16];
      uint32_t rem = remp[16];
      unsigned char u = *++cp;

      do {
	uint32_t v = convert_num[u];

	if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	  break;
	}
	n = n * 16 + v;
	u = *++cp;
      } while (isxdigit (u));
    }
    else {
      uint32_t max = maxp[8];
      uint32_t rem = remp[8];
      unsigned char u = *cp;

      while (u >= '0' && u <= '7') {
	uint32_t v = convert_num[u];

	if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	  break;
	}
	n = n * 8 + v;
	u = *++cp;
      }
    }
  }
  else if (isdigit (*cp)) {
    uint32_t max = maxp[10];
    uint32_t rem = remp[10];
    unsigned char u = *cp;

    do {
      uint32_t v = convert_num[u];

      if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	break;
      }
      n = n * 10 + v;
      u = *++cp;
    } while (isdigit (u));
  }
  else {
    cp = (unsigned char *) str;
  }
  if (endptr) {
    *endptr = (char *) cp;
  }
  return n;
}

/** \brief fast_strtou64
 *
 * \b Description
 *
 * Convert string to unsigned long
 *
 * \param str String to convert from
 * \param endptr optional endptr
 * \returns converted string
 */

uint64_t
fast_strtou64 (const char *str, char **endptr, int base)
{
  uint64_t n;
  unsigned char *cp = (unsigned char *) str;
#define	M(x)	(UINT64_C(18446744073709551615) / (x))
  static const uint64_t maxp[37] = {
    0, 0, M (2), M (3), M (4), M (5), M (6), M (7), M (8), M (9), M (10),
    M (11), M (12), M (13), M (14), M (15), M (16), M (17), M (18),
    M (19), M (20), M (21), M (22), M (23), M (24), M (25), M (26),
    M (27), M (28), M (29), M (30), M (31), M (32), M (33), M (34),
    M (35), M (36)
  };
#undef M
#define	R(x)	(UINT64_C(18446744073709551615) % (x))
  static const uint64_t remp[37] = {
    0, 0, R (2), R (3), R (4), R (5), R (6), R (7), R (8), R (9), R (10),
    R (11), R (12), R (13), R (14), R (15), R (16), R (17), R (18),
    R (19), R (20), R (21), R (22), R (23), R (24), R (25), R (26),
    R (27), R (28), R (29), R (30), R (31), R (32), R (33), R (34),
    R (35), R (36)
  };
#undef R

  while (isspace (*cp)) {
    cp++;
  }
  n = 0;
  if (UNLIKELY (base) &&
      (base != 10 || *cp < '1') &&
      (base != 8 || *cp != '0' || cp[1] < '0' || cp[1] > '7') &&
      (base != 16 || *cp != '0' || (cp[1] != 'x' && cp[1] != 'X'))) {
    if (base >= 2 && base <= 36) {
      uint64_t max = maxp[base];
      uint64_t rem = remp[base];

      if (base <= 10) {
	unsigned char u = *cp;
	unsigned char end = '0' + base;

	if (u >= '0' && u < end) {
	  do {
	    uint32_t v = convert_num[u];

	    if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	      break;
	    }
	    n = n * base + v;
	    u = *++cp;
	  } while (u >= '0' && u < end);
	}
	else {
	  cp = (unsigned char *) str;
	}
      }
      else {
	unsigned char u = uppercase[*cp];
	unsigned char end = 'A' + base - 10;

	if (valid_num[u] && u < end) {
	  do {
	    uint32_t v = convert_num[u];

	    if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	      break;
	    }
	    n = n * base + v;
	    u = uppercase[*++cp];
	  } while (valid_num[u] && u < end);
	}
	else {
	  cp = (unsigned char *) str;
	}
      }
    }
    else {
      cp = (unsigned char *) str;
    }
  }
  else if (*cp == '0') {
    cp++;
    if ((*cp == 'x' || *cp == 'X') && isxdigit (cp[1])) {
      uint64_t max = maxp[16];
      uint64_t rem = remp[16];
      unsigned char u = *++cp;

      do {
	uint32_t v = convert_num[u];

	if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	  break;
	}
	n = n * 16 + v;
	u = *++cp;
      } while (isxdigit (u));
    }
    else {
      uint64_t max = maxp[8];
      uint64_t rem = remp[8];
      unsigned char u = *cp;

      while (u >= '0' && u <= '7') {
	uint32_t v = convert_num[u];

	if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	  break;
	}
	n = n * 8 + v;
	u = *++cp;
      }
    }
  }
  else if (isdigit (*cp)) {
    uint64_t max = maxp[10];
    uint64_t rem = remp[10];
    unsigned char u = *cp;

    do {
      uint32_t v = convert_num[u];

      if (UNLIKELY (n >= max) && (n > max || (n == max && v > rem))) {
	break;
      }
      n = n * 10 + v;
      u = *++cp;
    } while (isdigit (u));
  }
  else {
    cp = (unsigned char *) str;
  }
  if (endptr) {
    *endptr = (char *) cp;
  }
  return n;
}

/** \brief fast_ftoa
 *
 * \b Description
 *
 * Convert float to ascii
 *
 * \param v float value
 * \param size precision
 * \param line pointer to result
 * \returns lenght string
 */

unsigned int
fast_ftoa (float v, int size, char *line)
{
  uint32_t q;
  int exp;
  int r;
  unsigned int l;
  char *s = line;
  union
  {
    float f;
    uint32_t u;
  } f;
  uint32_t lo;
  uint64_t qq;

  if (UNLIKELY (size <= 0 || size > PREC_FLT_NR)) {
    size = PREC_FLT_NR;
  }
  size--;

  f.f = v;
  if (f.u & 0x80000000) {
    *s++ = '-';
  }
  exp = (int) ((f.u >> 23) & 0xFF);
  if (UNLIKELY (exp == 0xFF)) {
    if ((f.u & 0x007FFFFF) == 0) {
      strcpy (s, "inf");
    }
    else {
      strcpy (s, "nan");
    }
    return (s + 3) - line;
  }
  else if (LIKELY (exp)) {
    exp += 25;
  }
  else {
    if (UNLIKELY ((f.u & 0x007FFFFF) == 0)) {
      strcpy (s, "0");
      return (s + 1) - line;
    }
    f.f *= 33554432.0;		/* 2^25 */
    exp = (int) ((f.u >> 23) & 0xFF);
  }
  q = (f.u & 0x007FFFFF) + 0x00800000;
  qq = mul_56 (q << 8, fpowers2[exp].mul, &lo);
  exp = fpowers2[exp].exp;
  if (size != PREC_FLT_NR - 1) {
    int n = size - (qq >= ipowers64[PREC_FLT_NR])
      - (qq >= ipowers64[PREC_FLT_NR + 1])
      - (qq >= ipowers64[PREC_FLT_NR + 2]);
    uint64_t d = qq / ipowers64[PREC_FLT_NR - 1 - n];
    uint64_t m = qq - d * ipowers64[PREC_FLT_NR - 1 - n];
    uint64_t h = ipowers64[PREC_FLT_NR - 2 - n] * 5;

    if (UNLIKELY (m == h)) {
      q = d + (lo != 0 || (d & 1));
    }
    else {
      q = d + (m > h);
    }
    exp += PREC_FLT_NR - n - 1;
    if (UNLIKELY (q >= ipowers32[size + 1])) {
      q = (q + 5) / 10;
      exp++;
    }
  }
  else if (qq >= ipowers64[PREC_FLT_NR + 2]) {
#if __WORDSIZE == 64
    uint64_t d = qq / 1000;
#else
    uint64_t d = div_1000 (qq);
#endif
    uint64_t m = qq - d * 1000;

    if (UNLIKELY (m == 500)) {
      q = d + (lo != 0 || (d & 1));
    }
    else {
      q = d + (m > 500);
    }
    exp += 3;
  }
  else if (qq >= ipowers64[PREC_FLT_NR + 1]) {
#if __WORDSIZE == 64
    uint64_t d = qq / 100;
#else
    uint64_t d = div_100 (qq);
#endif
    uint64_t m = qq - d * 100;

    if (UNLIKELY (m == 50)) {
      q = d + (lo != 0 || (d & 1));
    }
    else {
      q = d + (m > 50);
    }
    if (UNLIKELY (q >= ipowers32[PREC_FLT_NR])) {
      /* fp = 423718298 (9.99999999e-24) */
      q = (q + 5) / 10;
      exp++;
    }
    exp += 2;
  }
  else if (qq >= ipowers64[PREC_FLT_NR]) {
    /* Never reached */
#if __WORDSIZE == 64
    uint64_t d = qq / 10;
#else
    uint64_t d = div_10 (qq);
#endif
    uint64_t m = qq - d * 10;

    if (UNLIKELY (m == 5)) {
      q = d + (lo != 0 || (d & 1));
    }
    else {
      q = d + (m > 5);
    }
    exp++;
  }
  else {
    /* Never reached */
    q = qq;
  }

  r = 0;
  while ((q % 10) == 0) {
    r++;
    exp++;
    q /= 10;
  }
  if (exp >= 0 && exp <= r) {
    q *= ipowers32[exp];
    s += fast_uint32 (q, s);
  }
  else if (exp < 0 && exp > -(size - r + 5)) {
    exp += size - r;
    if (exp < 0) {
      *s++ = '0';
      *s++ = DECIMAL_POINT;
      while (++exp < 0) {
	*s++ = '0';
      }
      s += fast_uint32 (q, s);
    }
    else {
      l = fast_uint32 (q, s + 1);
      while (l-- != 0) {
	s[0] = s[1];
	s++;
	if (exp-- == 0) {
	  *s++ = DECIMAL_POINT;
	  s += l;
	  break;
	}
      }
    }
  }
  else {
    l = fast_uint32 (q, s + 1);
    s[0] = s[1];
    s++;
    if (l > 1) {
      *s = DECIMAL_POINT;
      s += l;
    }
    exp += size - r;

    if (exp) {
      *s++ = 'e';
      if (exp < 0) {
	*s++ = '-';
	exp = -exp;
      }
      else {
	*s++ = '+';
      }
      *s++ = '0' + (exp / 10);
      *s++ = '0' + (exp % 10);
    }
  }
  *s = '\0';
  return s - line;
}

/** \brief fast_dtoa
 *
 * \b Description
 *
 * Convert double to ascii
 *
 * \param v float value
 * \param size precision
 * \param line pointer to result
 * \returns lenght string
 */

unsigned int
fast_dtoa (double v, int size, char *line)
{
  uint64_t q;
  int exp;
  int r;
  unsigned int l;
  char *s = line;
  union
  {
    double d;
    uint64_t u;
  } d;
  uint32_t lo;

  if (UNLIKELY (size <= 0 || size > PREC_DBL_NR)) {
    size = PREC_DBL_NR;
  }
  size--;

  d.d = v;
  if (d.u & UINT64_C (0x8000000000000000)) {
    *s++ = '-';
  }
  exp = (int) ((d.u >> 52) & 0x7FF);
  if (UNLIKELY (exp == 0x7FF)) {
    if ((d.u & UINT64_C (0x000FFFFFFFFFFFFF)) == 0) {
      strcpy (s, "inf");
    }
    else {
      strcpy (s, "nan");
    }
    return (s + 3) - line;
  }
  else if (LIKELY (exp)) {
    exp += 54;
  }
  else {
    if (UNLIKELY ((d.u & UINT64_C (0x000FFFFFFFFFFFFF)) == 0)) {
      strcpy (s, "0");
      return (s + 1) - line;
    }
    d.d *= 18014398509481984.0;	/* 2^54 */
    exp = (int) ((d.u >> 52) & 0x7FF);
  }
  q = (d.u & UINT64_C (0x000FFFFFFFFFFFFF)) + UINT64_C (0x0010000000000000);
  q = mul_96 (q << 11, dpowers2[exp].mul1, dpowers2[exp].mul2, &lo);
  exp = dpowers2[exp].exp;
  if (size != PREC_DBL_NR - 1) {
    int n = size - (q >= ipowers64[PREC_DBL_NR])
      - (q >= ipowers64[PREC_DBL_NR + 1])
      - (q >= ipowers64[PREC_DBL_NR + 2]);
    uint64_t d = q / ipowers64[PREC_DBL_NR - 1 - n];
    uint64_t m = q - d * ipowers64[PREC_DBL_NR - 1 - n];
    uint64_t h = ipowers64[PREC_DBL_NR - 2 - n] * 5;

    if (UNLIKELY (m == h)) {
      q = d + (lo != 0 || (d & 1));
    }
    else {
      q = d + (m > h);
    }
    exp += PREC_DBL_NR - n - 1;
    if (UNLIKELY (q >= ipowers64[size + 1])) {
#if __WORDSIZE == 64
      q = (q + 5) / 10;
#else
      q = div_10 (q + 5);
#endif
      exp++;
    }
  }
  else if (q >= ipowers64[PREC_DBL_NR + 2]) {
#if __WORDSIZE == 64
    uint64_t d = q / 1000;
#else
    uint64_t d = div_1000 (q);
#endif
    uint64_t m = q - d * 1000;

    if (UNLIKELY (m == 500)) {
      q = d + (lo != 0 || (d & 1));
    }
    else {
      q = d + (m > 500);
    }
    exp += 3;
  }
  else if (q >= ipowers64[PREC_DBL_NR + 1]) {
#if __WORDSIZE == 64
    uint64_t d = q / 100;
#else
    uint64_t d = div_100 (q);
#endif
    uint64_t m = q - d * 100;

    if (UNLIKELY (m == 50)) {
      q = d + (lo != 0 || (d & 1));
    }
    else {
      q = d + (m > 50);
    }
    exp += 2;
  }
  else if (q >= ipowers64[PREC_DBL_NR]) {
#if __WORDSIZE == 64
    uint64_t d = q / 10;
#else
    uint64_t d = div_10 (q);
#endif
    uint64_t m = q - d * 10;

    if (UNLIKELY (m == 5)) {
      q = d + (lo != 0 || (d & 1));
    }
    else {
      q = d + (m > 5);
    }
    exp++;
  }

  r = 0;
#if __WORDSIZE == 64
  while ((q % 10) == 0) {
    r++;
    exp++;
    q /= 10;
  }
#else
  while ((q - div_10 (q) * 10) == 0) {
    r++;
    exp++;
    q = div_10 (q);
  }
#endif
  if (exp >= 0 && exp <= r) {
    q *= ipowers64[exp];
    s += fast_uint64 (q, s);
  }
  else if (exp < 0 && exp > -(size - r + 5 + (exp >= 100))) {
    exp += size - r;
    if (exp < 0) {
      *s++ = '0';
      *s++ = DECIMAL_POINT;
      while (++exp < 0) {
	*s++ = '0';
      }
      s += fast_uint64 (q, s);
    }
    else {
      l = fast_uint64 (q, s + 1);
      while (l-- != 0) {
	s[0] = s[1];
	s++;
	if (exp-- == 0) {
	  *s++ = DECIMAL_POINT;
	  s += l;
	  break;
	}
      }
    }
  }
  else {
    l = fast_uint64 (q, s + 1);
    s[0] = s[1];
    s++;
    if (l > 1) {
      *s = DECIMAL_POINT;
      s += l;
    }
    exp += size - r;

    if (exp) {
      *s++ = 'e';
      if (exp < 0) {
	*s++ = '-';
	exp = -exp;
      }
      else {
	*s++ = '+';
      }
      if (exp >= 100) {
	*s++ = '0' + (exp / 100);
	exp %= 100;
      }
      *s++ = '0' + (exp / 10);
      *s++ = '0' + (exp % 10);
    }
  }
  *s = '\0';
  return s - line;
}

/** \brief fast_strtof
 *
 * \b Description
 *
 * Convert string to float
 *
 * \param str string to convert
 * \param endptr optional endptr
 * \returns converted float value
 */

float
fast_strtof (const char *str, char **endptr)
{
  char *cp = (char *) str;
  int sign = 0;
  int esign = 0;
  int exp;
  int tmp;
  int c;
  uint64_t n;
  union
  {
    uint32_t u;
    float f;
  } tf;

  while (isspace (*cp)) {
    cp++;
  }
  if (*cp == '+') {
    cp++;
  }
  else if (*cp == '-') {
    sign = 1;
    cp++;
  }
  if (*cp == 'n' || *cp == 'N') {
    if (strncasecmp (cp, "nan", 3) == 0) {
      cp += strlen ("nan");
      if (endptr) {
	*endptr = cp;
      }
      if (*cp == '(') {
	cp++;
	while (isalpha (*cp) || isdigit (*cp) || *cp == '_') {
	  cp++;
	}
	if (*cp == ')') {
	  if (endptr) {
	    *endptr = cp + 1;
	  }
	}
      }
      tf.u = sign ? 0xFFC00000 : 0x7FC00000;
      return tf.f;
    }
    if (endptr) {
      *endptr = (char *) str;
    }
    return 0.0;
  }
  if (*cp == 'i' || *cp == 'I') {
    if (strncasecmp (cp, "inf", 3) == 0) {
      cp += strlen ("inf");
      if (strncasecmp (cp, "inity", strlen ("inity")) == 0) {
	cp += strlen ("inity");
      }
      if (endptr) {
	*endptr = cp;
      }
      tf.u = sign ? 0xFF800000 : 0x7F800000;
      return tf.f;
    }
    if (endptr) {
      *endptr = (char *) str;
    }
    return 0.0;
  }
  if (*cp == '0' && (cp[1] == 'x' || cp[1] == 'X')) {
    if (!isxdigit (cp[2]) && (cp[2] != DECIMAL_POINT || !isxdigit (cp[3]))) {
      if (endptr) {
	*endptr = &cp[1];
      }
      return 0.0 * (sign ? -1.0 : 1.0);
    }
    cp += 2;
    n = 0;
    exp = 0;
    c = 0;
    while (isxdigit (*cp)) {
      if (c < 16) {
	n = n * 16 + convert_num[*cp & 0xFFu];
      }
      else {
	exp += 4;
      }
      cp++;
      if (n) {
	c++;
      }
    }
    if (*cp == DECIMAL_POINT) {
      cp++;
      while (isxdigit (*cp)) {
	if (c < 16) {
	  n = n * 16 + convert_num[*cp & 0xFFu];
	  exp -= 4;
	}
	cp++;
	if (n) {
	  c++;
	}
      }
    }
    if (*cp == 'p' || *cp == 'P') {
      cp++;
      if (*cp == '+') {
	cp++;
      }
      else if (*cp == '-') {
	esign = 1;
	cp++;
      }
      tmp = 0;
      c = 0;
      while (isdigit (*cp)) {
	if (c < 5) {
	  tmp = tmp * 10 + (*cp - '0');
	}
	cp++;
	if (tmp) {
	  c++;
	}
      }
      if (esign) {
	exp -= tmp;
      }
      else {
	exp += tmp;
      }
    }
    if (endptr) {
      *endptr = cp;
    }
    /* sets ERANGE and returns HUGE_VAL on error */
    return ldexpf (n, exp) * (sign ? -1.0 : 1.0);
  }
  if (!isdigit (*cp) && (*cp != DECIMAL_POINT || !isdigit (cp[1]))) {
    if (endptr) {
      *endptr = (char *) str;
    }
    return 0.0;
  }
  n = 0;
  exp = 0;
  c = 0;
  while (isdigit (*cp)) {
    if (c < 19) {
      n = n * 10 + (*cp - '0');
    }
    else {
      exp++;
    }
    cp++;
    if (n) {
      c++;
    }
  }
  if (*cp == DECIMAL_POINT) {
    cp++;
    while (isdigit (*cp)) {
      if (c < 19) {
	n = n * 10 + (*cp - '0');
	exp--;
      }
      cp++;
      if (n) {
	c++;
      }
    }
  }
  if (*cp == 'e' || *cp == 'E') {
    cp++;
    if (*cp == '+') {
      cp++;
    }
    else if (*cp == '-') {
      esign = 1;
      cp++;
    }
    tmp = 0;
    c = 0;
    while (isdigit (*cp)) {
      if (c < 5) {
	tmp = tmp * 10 + (*cp - '0');
      }
      cp++;
      if (tmp) {
	c++;
      }
    }
    if (esign) {
      exp -= tmp;
    }
    else {
      exp += tmp;
    }
  }
  if (endptr) {
    *endptr = cp;
  }
  if (n && (exp >= -64 && exp <= 39)) {
#ifdef WIN
    unsigned int s = __builtin_clzll (n);
#else
#ifdef __GNUC__
#if __WORDSIZE == 64
    unsigned int s = __builtin_clzl (n);
#else
    unsigned int s = __builtin_clzll (n);
#endif
#else
    unsigned int s = calc_clz64 (n);
#endif
#endif

    /* sets ERANGE and returns HUGE_VAL on error */
    return ldexpf (mul_64 (n << s, fpowers10[exp + 64].mul) + 5,
		   fpowers10[exp + 64].exp + 32 - s) * (sign ? -1.0 : 1.0);
  }
  return 0.0 * (sign ? -1.0 : 1.0);
}

/** \brief fast_strtod
 *
 * \b Description
 *
 * Convert string to double
 *
 * \param str string to convert
 * \param endptr optional endptr
 * \returns converted double value
 */

double
fast_strtod (const char *str, char **endptr)
{
  char *cp = (char *) str;
  int sign = 0;
  int esign = 0;
  int exp;
  int tmp;
  int c;
  uint64_t n1;
  uint64_t n2;
  union
  {
    uint64_t u;
    double d;
  } td;

  while (isspace (*cp)) {
    cp++;
  }
  if (*cp == '+') {
    cp++;
  }
  else if (*cp == '-') {
    sign = 1;
    cp++;
  }
  if (*cp == 'n' || *cp == 'N') {
    if (strncasecmp (cp, "nan", 3) == 0) {
      cp += strlen ("nan");
      if (endptr) {
	*endptr = cp;
      }
      if (*cp == '(') {
	cp++;
	while (isalpha (*cp) || isdigit (*cp) || *cp == '_') {
	  cp++;
	}
	if (*cp == ')') {
	  if (endptr) {
	    *endptr = cp + 1;
	  }
	}
      }
      td.u =
	sign ? UINT64_C (0xFFF8000000000000) : UINT64_C (0x7FF8000000000000);
      return td.d;
    }
    if (endptr) {
      *endptr = (char *) str;
    }
    return 0.0;
  }
  if (*cp == 'i' || *cp == 'I') {
    if (strncasecmp (cp, "inf", 3) == 0) {
      cp += strlen ("inf");
      if (strncasecmp (cp, "inity", strlen ("inity")) == 0) {
	cp += strlen ("inity");
      }
      if (endptr) {
	*endptr = cp;
      }
      td.u =
	sign ? UINT64_C (0xFFF0000000000000) : UINT64_C (0x7FF0000000000000);
      return td.d;
    }
    if (endptr) {
      *endptr = (char *) str;
    }
    return 0.0;
  }
  if (*cp == '0' && (cp[1] == 'x' || cp[1] == 'X')) {
    if (!isxdigit (cp[2]) && (cp[2] != DECIMAL_POINT || !isxdigit (cp[3]))) {
      if (endptr) {
	*endptr = &cp[1];
      }
      return 0.0 * (sign ? -1.0 : 1.0);
    }
    cp += 2;
    n1 = 0;
    exp = 0;
    c = 0;
    while (isxdigit (*cp)) {
      if (c < 16) {
	n1 = n1 * 16 + convert_num[*cp & 0xFFu];
      }
      else {
	exp += 4;
      }
      cp++;
      if (n1) {
	c++;
      }
    }
    if (*cp == DECIMAL_POINT) {
      cp++;
      while (isxdigit (*cp)) {
	if (c < 16) {
	  n1 = n1 * 16 + convert_num[*cp & 0xFFu];
	  exp -= 4;
	}
	cp++;
	if (n1) {
	  c++;
	}
      }
    }
    if (*cp == 'p' || *cp == 'P') {
      cp++;
      if (*cp == '+') {
	cp++;
      }
      else if (*cp == '-') {
	esign = 1;
	cp++;
      }
      tmp = 0;
      c = 0;
      while (isdigit (*cp)) {
	if (c < 5) {
	  tmp = tmp * 10 + (*cp - '0');
	}
	cp++;
	if (tmp) {
	  c++;
	}
      }
      if (esign) {
	exp -= tmp;
      }
      else {
	exp += tmp;
      }
    }
    if (endptr) {
      *endptr = cp;
    }
    /* sets ERANGE and returns HUGE_VAL on error */
    return ldexp (n1, exp) * (sign ? -1.0 : 1.0);
  }
  if (!isdigit (*cp) && (*cp != DECIMAL_POINT || !isdigit (cp[1]))) {
    if (endptr) {
      *endptr = (char *) str;
    }
    return 0.0;
  }
  n1 = 0;
  n2 = 0;
  exp = 0;
  c = 0;
  while (c < 19 && isdigit (*cp)) {
    n2 = n2 * 10 + (*cp - '0');
    cp++;
    if (n2) {
      c++;
    }
  }
  while (isdigit (*cp)) {
    if (c < 38) {
      mul_10_add (&n1, &n2, *cp - '0');
    }
    else {
      exp++;
    }
    cp++;
    if (n1 || n2) {
      c++;
    }
  }
  if (*cp == DECIMAL_POINT) {
    cp++;
    while (c < 19 && isdigit (*cp)) {
      n2 = n2 * 10 + (*cp - '0');
      exp--;
      cp++;
      if (n2) {
	c++;
      }
    }
    while (isdigit (*cp)) {
      if (c < 38) {
	mul_10_add (&n1, &n2, *cp - '0');
	exp--;
      }
      cp++;
      if (n1 || n2) {
	c++;
      }
    }
  }
  if (*cp == 'e' || *cp == 'E') {
    cp++;
    if (*cp == '+') {
      cp++;
    }
    else if (*cp == '-') {
      esign = 1;
      cp++;
    }
    tmp = 0;
    c = 0;
    while (isdigit (*cp)) {
      if (c < 5) {
	tmp = tmp * 10 + (*cp - '0');
      }
      cp++;
      if (tmp) {
	c++;
      }
    }
    if (esign) {
      exp -= tmp;
    }
    else {
      exp += tmp;
    }
  }
  if (endptr) {
    *endptr = cp;
  }
  c = 64;
  if (n1 == 0) {
    n1 = n2;
    n2 = 0;
    c = 0;
  }
  if (n1 && (exp >= -362 && exp <= 309)) {
#ifdef WIN
    unsigned int s = __builtin_clzll (n1);
#else
#ifdef __GNUC__
#if __WORDSIZE == 64
    unsigned int s = __builtin_clzl (n1);
#else
    unsigned int s = __builtin_clzll (n1);
#endif
#else
    unsigned int s = calc_clz64 (n1);
#endif
#endif

    n1 = (n1 << s) | (n2 >> (64 - s));
    n2 <<= s;

    uint32_t lo;
    uint64_t r = mul_96 (n1, dpowers10[exp + 362].mul1,
			 dpowers10[exp + 362].mul2 + 1, &lo);

    if (n2) {
      uint32_t lo2;
      uint64_t l = mul_96 (n2, dpowers10[exp + 362].mul1,
			   dpowers10[exp + 362].mul2 + 1, &lo2) >> 32;
      lo += l;
      if (lo < l) {
	r++;
	if (r == 0) {
	  r--;
	}
      }
    }

    return ldexp (r, dpowers10[exp + 362].exp + 64 - s + c) *
      (sign ? -1.0 : 1.0);
  }
  return 0.0 * (sign ? -1.0 : 1.0);
}

#endif				/* __FAST_STDIO_H */
