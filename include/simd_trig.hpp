#ifndef MINTRIG_HPP
#define MINTRIG_HPP

#include <immintrin.h>

//Most ninja tricks used here:
//http://fastcpp.blogspot.fr/2011/03/changing-sign-of-float-values-using-sse.html
//http://www.songho.ca/misc/sse/sse.html
//http://markplusplus.wordpress.com/2007/03/14/fast-sse-select-operation/
//http://www.masmforum.com/board/index.php?PHPSESSID=786dd40408172108b65a5a36b09c88c0&topic=9515.0
//http://cbloomrants.blogspot.fr/2010/11/11-20-10-function-approximation-by_20.html
//http://assemblyrequired.crashworks.org/2009/10/16/timing-square-root/
//http://nghiaho.com/?p=997
//http://www.researchgate.net/publication/3321724_Efficient_approximations_for_the_arctangent_function
//http://www.ganssle.com/approx/approx.pdf
//http://forum.allaboutcircuits.com/newsgroups/viewtopic.php?t=68185

//static constants are slower, than standard constants. To be replaced with constexpr
    const __m256 AVX_SIGNMASK_PS =  _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));

namespace FMA{

    __m256 sqrt(__m256 squared);
    __m256 sqrt_ps(__m256 squared);

    __m256 sqrt_2pi_ps(__m256 squared);
    __m256 cos_2pi_poly_ps(__m256 x);
    __m256 cos_2pi_ps(__m256 angle);
    __m256 sin_2pi_poly_ps(__m256 x);
    __m256 sin_2pi_ps(__m256 angle);
    void sincos_2pi_ps(__m256, __m256*, __m256*);
}


__m256 FMA::sqrt_ps(__m256 squared){return _mm256_sqrt_ps(squared);}


// FMA single precision implementation of sin(2πx) and cos(2πx)
inline __m256 FMA::cos_2pi_poly_ps(__m256 x) {
    const __m256 c1 = _mm256_set1_ps(9.99996748424514608493916615939939119e-01f);
    const __m256 c2 = _mm256_set1_ps(-1.97363223756305008622054184241628274e+01f);
    const __m256 c3 = _mm256_set1_ps(6.46807901818389923284666266409658356e+01f);
    const __m256 c4 = _mm256_set1_ps(-7.84122201283542750695287203100567467e+01f);
    const __m256 x2 = _mm256_mul_ps(x, x);

    // Using FMA instructions for more efficient computation
    __m256 result = _mm256_fmadd_ps(c4, x2, c3); // c3 + c4 * x2
    result = _mm256_fmadd_ps(result, x2, c2);    // c2 + (c3 + c4 * x2) * x2
    result = _mm256_fmadd_ps(result, x2, c1);    // c1 + (c2 + (c3 + c4 * x2) * x2) * x2

    return result;
}

inline __m256 FMA::sin_2pi_poly_ps(__m256 x) {
    const __m256 c1 = _mm256_set1_ps(6.283167562848116546596361730703918648688e+00f); //constexpr
    const __m256 c2 = _mm256_set1_ps(-4.133751762600560324647974921076836107669e+01f);
    const __m256 c3 = _mm256_set1_ps(8.135167773428457047365429407228838464935e+01f);
    const __m256 c4 = _mm256_set1_ps(-7.108735817320177431942998552177234798743e+01f);
    const __m256 x2 = _mm256_mul_ps(x, x);

    // Using FMA instructions for odd powers
    __m256 result = _mm256_fmadd_ps(c4, x2, c3); // c3 + c4 * x^2
    result = _mm256_fmadd_ps(result, x2, c2);    // c2 + (c3 + c4 * x^2) * x^2
    result = _mm256_fmadd_ps(result, x2, c1);    // c1 + (c2 + (c3 + c4 * x^2) * x^2) * x^2

    return _mm256_mul_ps(result, x); // Multiply by x to apply odd powers
}

__m256 FMA::cos_2pi_ps(__m256 angle){
    //clamp to the range [0, 1]
    angle=_mm256_andnot_ps(AVX_SIGNMASK_PS,angle); //take absolute value
    angle = _mm256_sub_ps(angle, _mm256_mul_ps(_mm256_floor_ps(angle), _mm256_set1_ps(1.0f)));  //fmod(angle, 1)

    __m256 cosangle=angle;
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.25f), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(0.5),angle))));
    cosangle=_mm256_xor_ps(cosangle,_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.5f), _CMP_GE_OQ), AVX_SIGNMASK_PS));
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.75f), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(1),angle))));

    __m256 result=FMA::cos_2pi_poly_ps(cosangle);

    result=_mm256_xor_ps(result,_mm256_and_ps(_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.25f), _CMP_GE_OQ),_mm256_cmp_ps(angle,_mm256_set1_ps(0.75f), _CMP_LT_OQ)), AVX_SIGNMASK_PS));
    return result;
}

__m256 FMA::sin_2pi_ps(__m256 angle) {
    angle = _mm256_sub_ps(angle, _mm256_floor_ps(angle));

    __m256 sinangle = angle;

    sinangle = _mm256_xor_ps(sinangle, _mm256_and_ps(_mm256_cmp_ps(angle, _mm256_set1_ps(0.25f), _CMP_GE_OQ), _mm256_xor_ps(sinangle, _mm256_sub_ps(_mm256_set1_ps(0.5f), angle))));
    sinangle = _mm256_xor_ps(sinangle,_mm256_and_ps(_mm256_cmp_ps(angle, _mm256_set1_ps(0.5f), _CMP_GE_OQ), AVX_SIGNMASK_PS));
    sinangle = _mm256_xor_ps(sinangle, _mm256_and_ps(_mm256_cmp_ps(angle, _mm256_set1_ps(0.75f), _CMP_GE_OQ),_mm256_xor_ps(sinangle, _mm256_sub_ps(_mm256_set1_ps(1.0f), angle))));

    __m256 result = FMA::sin_2pi_poly_ps(sinangle);

    result = _mm256_xor_ps(result, _mm256_and_ps(_mm256_cmp_ps(angle, _mm256_set1_ps(0.5f), _CMP_GE_OQ), AVX_SIGNMASK_PS));

    return result;
}


void FMA::sincos_2pi_ps(__m256 angle, __m256* sin, __m256* cos) {

    angle = _mm256_sub_ps(angle, _mm256_floor_ps(angle));

    __m256 cosangle = angle;
    cosangle = _mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle, _mm256_set1_ps(0.25f), _CMP_GE_OQ), _mm256_xor_ps(cosangle, _mm256_sub_ps(_mm256_set1_ps(0.5f), angle))));
    cosangle = _mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle, _mm256_set1_ps(0.5f), _CMP_GE_OQ), AVX_SIGNMASK_PS));
    cosangle = _mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle, _mm256_set1_ps(0.75f), _CMP_GE_OQ), _mm256_xor_ps(cosangle, _mm256_sub_ps(_mm256_set1_ps(1.0f), angle))));
    __m256 sinangle = cosangle;

    __m256 cos_result = FMA::cos_2pi_poly_ps(cosangle);
    cos_result = _mm256_xor_ps(cos_result, _mm256_and_ps(_mm256_and_ps(_mm256_cmp_ps(angle, _mm256_set1_ps(0.25f), _CMP_GE_OQ), _mm256_cmp_ps(angle, _mm256_set1_ps(0.75f), _CMP_LT_OQ)),AVX_SIGNMASK_PS));

    __m256 sin_result = FMA::sin_2pi_poly_ps(sinangle);
    sin_result = _mm256_xor_ps(sin_result, _mm256_and_ps(_mm256_cmp_ps(angle, _mm256_set1_ps(0.5f), _CMP_GE_OQ), AVX_SIGNMASK_PS));

    *cos = cos_result;
    *sin = sin_result;
return;
}



#endif // MINTRIG_HPP
