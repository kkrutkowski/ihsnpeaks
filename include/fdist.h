#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <float.h>
#include <math.h>

#define locEPS (1000.0 * DBL_EPSILON)

// hyp2f1 to be implemented

long double psi(long double xx)
{
    //                    Evaluation of the digamma function
    //
    //                          -----------
    //
    //    Psi(xx) is assigned the value 0 when the digamma function cannot
    //    be computed.
    //
    //    The main computation involves evaluation of rational chebyshev
    //    approximations published in math. Comp. 27, 123-127(1973) By
    //    cody, strecok and thacher.
    //
    //    ----------------------------------------------------------------
    //    Psi was written at Argonne National Laboratory for the FUNPACK
    //    package of special function subroutines. Psi was modified by
    //    A.H. Morris (nswc).

    long double aug, den, dx0, sgn, upper, w, x, xmax1, xmx0, xsmall, z;
    long double p1[7] = {0.895385022981970e-02, 0.477762828042627e+01,
                    0.142441585084029e+03, 0.118645200713425e+04,
                    0.363351846806499e+04, 0.413810161269013e+04,
                    0.130560269827897e+04};
    long double q1[6] = {0.448452573429826e+02, 0.520752771467162e+03,
                    0.221000799247830e+04, 0.364127349079381e+04,
                    0.190831076596300e+04, 0.691091682714533e-05};
    long double p2[4] = {-0.212940445131011e+01, -0.701677227766759e+01,
                    -0.448616543918019e+01, -0.648157123766197e+00};
    long double q2[4] = {0.322703493791143e+02, 0.892920700481861e+02,
                    0.546117738103215e+02, 0.777788548522962e+01};
    int nq, i;
    dx0 = 1.461632144968362341262659542325721325;
    xmax1 = 4503599627370496.0;
    xsmall = 1e-9;
    x = xx;
    aug = 0.0;

    if (x < 0.5) {
        if (fabs(x) <= xsmall) {
            if (x == 0.) {
                return 0.0;
            }
            aug = -1./x;
        } else {
            // 10
            w = -x;
            sgn = M_PI / 4;
            if (w <= 0.) {
                w = -w;
                sgn = -sgn;
            }
            // 20
            if (w >= xmax1) {
                return 0.0;
            }
            w -= (int)w;
            nq = (int)(w*4.0);
            w = 4.*(w - 0.25*nq);

            if (nq % 2 == 1) {
                w = 1. - w;
            }
            z = (M_PI / 4.)*w;

            if ((nq / 2) % 2 == 1) {
                sgn = -sgn;
            }
            if ((((nq + 1) / 2) % 2) == 1) {
                aug = sgn * (tan(z)*4.);
            } else {
                if (z == 0.) {
                    return 0.0;
                }
                aug = sgn * (4./tan(z));
            }
        }
        x = 1 - x;
    }

    if (x <= 3.0) {
        // 50
        den = x;
        upper = p1[0]*x;
        for (i = 0; i < 5; i++)
        {
            den = (den + q1[i])*x;
            upper = (upper + p1[i+1])*x;
        }
        den = (upper + p1[6]) / (den + q1[5]);
        xmx0 = x - dx0;
        return (den * xmx0) + aug;
    } else {
        // 70
        if (x < xmax1) {
            w = 1. / (x*x);
            den = w;
            upper = p2[0]*w;

            for (i = 0; i < 3; i++) {
                den = (den + q2[i])*w;
                upper = (upper + p2[i+1])*w;
            }
            aug += upper / (den + q2[3]) - 0.5/x;
        }
        return aug + log(x);
    }
}

long double alnrel(long double a)
{
    //    Evaluation of the function ln(1 + a)

    long double p[3] = {-0.129418923021993e+01, 0.405303492862024e+00, -0.178874546012214e-01};
    long double q[3] = {-0.162752256355323e+01, 0.747811014037616e+00, -0.845104217945565e-01};
    long double t, t2, w;

    if (fabs(a) > 0.375) {
        return logl(1. + a);
    } else {
        t = a / (a + 2.);
        t2 = t*t;
        w = ((p[2]*t2 + p[1])*t2 + p[0])*t2 + 1.;
        w /= ((q[2]*t2 + q[1])*t2 + q[0])*t2 + 1.;
        return 2.0*t*w;
    }
}

long double gamln1(long double a)
{
    //    Evaluation of ln(gamma(1 + a)) for -0.2 <= A <= 1.25

    long double bot, top, w, x;

    const long double p[7] = { .577215664901533e+00,  .844203922187225e+00,
                         -.168860593646662e+00, -.780427615533591e+00,
                         -.402055799310489e+00, -.673562214325671e-01,
                         -.271935708322958e-02};
    const long double q[6] = {.288743195473681e+01, .312755088914843e+01,
                         .156875193295039e+01, .361951990101499e+00,
                         .325038868253937e-01, .667465618796164e-03};
    const long double r[6] = {.422784335098467e+00, .848044614534529e+00,
                         .565221050691933e+00, .156513060486551e+00,
                         .170502484022650e-01, .497958207639485e-03};
    const long double s[5] = {.124313399877507e+01, .548042109832463e+00,
                         .101552187439830e+00, .713309612391000e-02,
                         .116165475989616e-03};

    if (a < 0.6) {
        top = ((((((p[6]
                   )*a+p[5]
                  )*a+p[4]
                 )*a+p[3]
                )*a+p[2]
               )*a+p[1]
              )*a+p[0];
        bot = ((((((q[5]
                   )*a+q[4]
                  )*a+q[3]
                 )*a+q[2]
                )*a+q[1]
               )*a+q[0]
              )*a+1.;
        w = top/bot;
        return -a*w;
    } else {
        x = (a - 0.5) - 0.5;
        top = (((((r[5]
                  )*x+r[4]
                 )*x+r[3]
                )*x+r[2]
               )*x+r[1]
              )*x+r[0];
        bot = (((((s[4]
                  )*x+s[3]
                 )*x+s[2]
                )*x+s[1]
               )*x+s[0]
              )*x+1.;
        w = top/bot;
        return x*w;
    }
}

long double gamln(long double a)
{
    //    Evaluation of ln(gamma(a)) for positive a

    long double t, w, d = .418938533204673;
    int i,n;
    const long double c[6] = {.833333333333333e-01, -.277777777760991e-02,
                         .793650666825390e-03, -.595202931351870e-03,
                         .837308034031215e-03, -.165322962780713e-02};

    if (a <= 0.8) { return gamln1(a) - logl(a); }

    if (a <= 2.25) {
        t = (a-0.5) - 0.5;
        return gamln1(t);
    }

    if (a < 10) {
        n = (int)(a - 1.25);
        t = a;
        w = 1.0;
        for (i = 0; i < n; i++)
        {
            t -= 1.0;
            w *= t;
        }
        return gamln1(t-1.) + logl(w);
    }
    t = pow(1/a, 2);
    w = (((((c[5]*t+c[4])*t+c[3])*t+c[2])*t+c[1])*t+c[0])/a;
    return (d + w) + (a-0.5)*(logl(a) - 1.);
}

long double algdiv(long double a, long double b)
{
    //         Computation of ln(gamma(b)/gamma(a+b)) when b >= 8
    //
    //                             --------
    //
    //         In this algorithm, del(x) is the function defined by
    //         Ln(gamma(x)) = (x - 0.5)*ln(x) - x + 0.5*ln(2*pi) + del(x).
    //

    long double c, d, h, s11, s3, s5, s7, s9, t, u, v, w, x, x2;
    long double carr[6] = {0.833333333333333e-01, -0.277777777760991e-02,
                      0.793650666825390e-03, -0.595202931351870e-03,
                      0.837308034031215e-03, -0.165322962780713e-02};

    if (a > b) {
        h = b / a;
        c = 1./(1. + h);
        x = h/(1. + h);
        d = a + (b - 0.5);
    } else {
        h = a / b;
        c = h/(1. + h);
        x = 1./(1. + h);
        d = b + (a - 0.5);
    }
    // Set sn = (1 - x**n)/(1 - x)
    x2 = x*x;
    s3 = 1. + (x + x2);
    s5 = 1. + (x + x2*s3);
    s7 = 1. + (x + x2*s5);
    s9 = 1. + (x + x2*s7);
    s11 = 1. + (x + x2*s9);

    // Set w = del(b) - del(a + b)
    t = pow((1. / b), 2);
    w = (((((carr[5]*s11
            )*t + carr[4]*s9
           )*t + carr[3]*s7
          )*t + carr[2]*s5
         )*t + carr[1]*s3
        )*t + carr[0];
    w *= c / b;
    // Combine the results
    u = d * alnrel(a / b);
    v = a * (logl(b) - 1.);
    return (u > v ? (w - v) - u : (w - u) - v);
}

long double bcorr(long double a0, long double b0)
{
    //    Evaluation of  del(a0) + del(b0) - del(a0 + b0)  where
    //    ln(gamma(a)) = (a - 0.5)*ln(a) - a + 0.5*ln(2*pi) + del(a).
    //    It is assumed that a0 >= 8 And b0 >= 8.

    long double a,b,c,h,s11,s3,s5,s7,s9,t,w,x,x2;
    long double carr[6] = {0.833333333333333e-01, -0.277777777760991e-02,
                      0.793650666825390e-03, -0.595202931351870e-03,
                      0.837308034031215e-03, -0.165322962780713e-02};

    a = fmin(a0, b0);
    b = fmax(a0, b0);
    h = a / b;
    c = h/(1. + h);
    x = 1./(1. + h);
    x2 = x*x;
    //  Set sn = (1 - x**n)/(1 - x)
    s3 = 1. + (x + x2);
    s5 = 1. + (x + x2*s3);
    s7 = 1. + (x + x2*s5);
    s9 = 1. + (x + x2*s7);
    s11 = 1. + (x + x2*s9);
    // Set w = del(b) - del(a + b)
    t = pow((1. / b), 2);
    w = (((((carr[5]*s11
            )*t + carr[4]*s9
           )*t + carr[3]*s7
          )*t + carr[2]*s5
         )*t + carr[1]*s3
        )*t + carr[0];
    w *= c / b;
    // Compute  del(a) + w
    t = pow((1. / a), 2);
    return ((((((carr[5])*t + carr[4]
               )*t + carr[3]
              )*t + carr[2]
             )*t + carr[1]
            )*t + carr[0]
           )/a + w;
}

long double gsumln(long double a, long double b)
{
    //     Evaluation of the function ln(gamma(a + b))
    //     for 1 <= A <= 2  And  1 <= B <= 2

    long double x;

    x = a + b - 2;
    if (x <= 0.25) {
        return gamln1(1. + x);
    }

    if (x <= 1.25) {
        return gamln1(x) + alnrel(x);
    }

    return gamln1(x - 1.) + logl(x*(1. + x));
}

long double betaln(long double a0, long double b0)
{
    //    Evaluation of the logarithm of the beta function

    long double a, b, c, h, u, v, w, z;
    long double e = 0.918938533204673;
    int i, n;

    a = fmin(a0, b0);
    b = fmax(a0, b0);

    if (a >= 8.0) {
        w = bcorr(a, b);
        h = a / b;
        c = h/(1. + h);
        u = -(a - 0.5)*logl(c);
        v = b*alnrel(h);
        if (u > v) {
            return (((-0.5*logl(b)+e)+w)-v) - u;
        } else {
            return (((-0.5*logl(b)+e)+w)-u) - v;
        }
    }
    if (a < 1) {
        if (b > 8) {
            return gamln(a) + algdiv(a,b);
        } else {
            return gamln(a) + (gamln(b) - gamln(a+b));
        }
    }

    if (a <= 2) {
        if (b <= 2) {
            return gamln(a) + gamln(b) - gsumln(a, b);
        }
        if (b >= 8) {
            return gamln(a) + algdiv(a, b);
        }
        w = 0.;
    }

    if (a > 2) {
        if (b <= 1000) {
            n = (int)(a - 1.);
            w = 1.;
            for (i = 0; i < n; i++) {
                a -= 1.0;
                h = a / b;
                w *= h/(1.+h);
            }
            w = logl(w);
            if (b >= 8.0) {
                return w + gamln(a) + algdiv(a, b);
            }
        } else {
            n = (int)(a - 1.);
            w = 1.0;
            for (i = 0; i < n; i++) {
                a -= 1.0;
                w *= a/(1. + (a/b));
            }
            return (logl(w) - n*logl(b)) + (gamln(a) + algdiv(a, b));
        }
    }
    n = (int)(b - 1.);
    z = 1.0;
    for (i = 0; i < n; i++) {
        b -= 1.0;
        z *= b / (a + b);
    }
    return w + logl(z) + (gamln(a) + gamln(b) - gsumln(a, b));
}


void sf_complex_log(double zr, double zi, double *log1_r, double *log1_i) {
    double complex z = zr + zi * I;
    double complex log_z = clog(z);

    *log1_r = creal(log_z);
    *log1_i = cimag(log_z);
}

void sf_complex_logsin(double zr, double zi, double *lszr, double *lszi) {
    double complex z = zr + zi * I;
    double complex sin_z = csin(z);
    double complex log_sin_z = clog(sin_z);

    *lszr = creal(log_sin_z);
    *lszi = cimag(log_sin_z);
}

/* coefficients for gamma=7, kmax=8  Lanczos method */
static double lanczos_7_c[9] = {
    0.99999999999980993227684700473478,
    676.520368121885098567009190444019,
    -1259.13921672240287047156078755283,
    771.3234287776530788486528258894,
    -176.61502916214059906584551354,
    12.507343278686904814458936853,
    -0.13857109526572011689554707,
    9.984369578019570859563e-6,
    1.50563273514931155834e-7
};

double gsl_sf_angle_restrict_symm(double theta) {
    // synthetic extended precision constants
    const double P1 = 4 * 7.8539812564849853515625e-01;
    const double P2 = 4 * 3.7748947079307981766760e-08;
    const double P3 = 4 * 2.6951514290790594840552e-15;
    const double TwoPi = 2 * (P1 + P2 + P3);

    const double y = (double)((theta > 0) - (theta < 0)) * 2 * floor(fabs(theta) / TwoPi);
    double r = ((theta - y * P1) - y * P2) - y * P3;

    if (r > M_PI) {
        r = (((r - 2 * P1) - 2 * P2) - 2 * P3);  // r-TwoPi
    } else if (r < -M_PI) {
        r = (((r + 2 * P1) + 2 * P2) + 2 * P3); // r+TwoPi
    }

    if (fabs(theta) > 0.0625 / DBL_EPSILON) {
        return NAN;
    }

    return r;
}

void lngamma_lanczos_complex(double zr, double zi, double *yr, double *yi) {
    int k;
    double log1_r, log1_i;
    double logAg_r, logAg_i;
    double Ag_r, Ag_i;

    zr -= 1.0; // Lanczos writes z! instead of Gamma(z)

    Ag_r = lanczos_7_c[0];
    Ag_i = 0.0;
    for (k = 1; k <= 8; k++) {
        double R = zr + k;
        double J = zi;
        double a = lanczos_7_c[k] / (R * R + J * J);
        Ag_r += a * R;
        Ag_i -= a * J;
    }

    sf_complex_log(zr + 7.5, zi, &log1_r, &log1_i);
    sf_complex_log(Ag_r, Ag_i, &logAg_r, &logAg_i);

    // (z+0.5)*log(z+7.5) - (z+7.5) + LogRootTwoPi_ + log(Ag(z))
    *yr = (zr + 0.5) * log1_r - zi * log1_i - (zr + 7.5) + log(sqrt(2 * M_PI)) + logAg_r;
    *yi = zi * log1_r + (zr + 0.5) * log1_i - zi + logAg_i;
}

double hyperg_2F1_series(double a, double b, double c, double x) {
    double sum_pos = 1.0;
    double sum_neg = 0.0;
    double del_pos = 1.0;
    double del_neg = 0.0;
    double del = 1.0;
    double del_prev;
    double k = 0.0;
    int i = 0;

    if (fabs(c) < DBL_EPSILON) {
        return 0.0; /* FIXME: ?? */
    }

    do {
        if (++i > 30000) {
            return sum_pos - sum_neg;
        }
        del_prev = del;
        del *= (a + k) * (b + k) * x / ((c + k) * (k + 1.0));  /* Gauss series */

        if (del > 0.0) {
            del_pos = del;
            sum_pos += del;
        } else if (del == 0.0) {
            /* Exact termination (a or b was a negative integer). */
            del_pos = 0.0;
            del_neg = 0.0;
            break;
        } else {
            del_neg = -del;
            sum_neg -= del;
        }

        if (fabs(del_prev / (sum_pos - sum_neg)) < DBL_EPSILON &&
            fabs(del / (sum_pos - sum_neg)) < DBL_EPSILON) {
            break;
        }

        k += 1.0;
    } while (fabs((del_pos + del_neg) / (sum_pos - sum_neg)) > DBL_EPSILON);

    return sum_pos - sum_neg;
}

double hyperg_2F1_luke(double a, double b, double c, double xin) {
    const double RECUR_BIG = 1.0e+50;
    const int nmax = 20000;
    int n = 3;
    const double x = -xin;
    const double x3 = x * x * x;
    const double t0 = a * b / c;
    const double t1 = (a + 1.0) * (b + 1.0) / (2.0 * c);
    const double t2 = (a + 2.0) * (b + 2.0) / (2.0 * (c + 1.0));
    double F = 1.0;
    double prec;

    double Bnm3 = 1.0;                                  /* B0 */
    double Bnm2 = 1.0 + t1 * x;                         /* B1 */
    double Bnm1 = 1.0 + t2 * x * (1.0 + t1 / 3.0 * x);  /* B2 */

    double Anm3 = 1.0;                                                      /* A0 */
    double Anm2 = Bnm2 - t0 * x;                                            /* A1 */
    double Anm1 = Bnm1 - t0 * (1.0 + t2 * x) * x + t0 * t1 * (c / (c + 1.0)) * x * x;   /* A2 */

    while (1) {
        double npam1 = n + a - 1;
        double npbm1 = n + b - 1;
        double npcm1 = n + c - 1;
        double npam2 = n + a - 2;
        double npbm2 = n + b - 2;
        double npcm2 = n + c - 2;
        double tnm1 = 2 * n - 1;
        double tnm3 = 2 * n - 3;
        double tnm5 = 2 * n - 5;
        double n2 = n * n;
        double F1 = (3.0 * n2 + (a + b - 6) * n + 2 - a * b - 2 * (a + b)) / (2 * tnm3 * npcm1);
        double F2 = -(3.0 * n2 - (a + b + 6) * n + 2 - a * b) * npam1 * npbm1 / (4 * tnm1 * tnm3 * npcm2 * npcm1);
        double F3 = (npam2 * npam1 * npbm2 * npbm1 * (n - a - 2) * (n - b - 2)) / (8 * tnm3 * tnm3 * tnm5 * (n + c - 3) * npcm2 * npcm1);
        double E = -npam1 * npbm1 * (n - c - 1) / (2 * tnm3 * npcm2 * npcm1);

        double An = (1.0 + F1 * x) * Anm1 + (E + F2 * x) * x * Anm2 + F3 * x3 * Anm3;
        double Bn = (1.0 + F1 * x) * Bnm1 + (E + F2 * x) * x * Bnm2 + F3 * x3 * Bnm3;
        double r = An / Bn;

        prec = fabs((F - r) / F);
        F = r;

        if (prec < DBL_EPSILON || n > nmax) break;

        if (fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
            An /= RECUR_BIG;
            Bn /= RECUR_BIG;
            Anm1 /= RECUR_BIG;
            Bnm1 /= RECUR_BIG;
            Anm2 /= RECUR_BIG;
            Bnm2 /= RECUR_BIG;
            Anm3 /= RECUR_BIG;
            Bnm3 /= RECUR_BIG;
        } else if (fabs(An) < 1.0 / RECUR_BIG || fabs(Bn) < 1.0 / RECUR_BIG) {
            An *= RECUR_BIG;
            Bn *= RECUR_BIG;
            Anm1 *= RECUR_BIG;
            Bnm1 *= RECUR_BIG;
            Anm2 *= RECUR_BIG;
            Bnm2 *= RECUR_BIG;
            Anm3 *= RECUR_BIG;
            Bnm3 *= RECUR_BIG;
        }

        n++;
        Bnm3 = Bnm2;
        Bnm2 = Bnm1;
        Bnm1 = Bn;
        Anm3 = Anm2;
        Anm2 = Anm1;
        Anm1 = An;
    }

    return F;
}

double hyperg_2F1_reflect(double a, double b, double c, double x) {
    const double d = c - a - b;
    const int intd = floor(d + 0.5);
    const int d_integer = (fabs(d - intd) < locEPS);

    if (d_integer) {
        const double ln_omx = log(1.0 - x);
        const double ad = fabs(d);
        double sgn_2;
        double F1, F2;
        double d1, d2;
        double lng_c, lng_ad2, lng_bd2;

        if (d >= 0.0) {
            d1 = d;
            d2 = 0.0;
        } else {
            d1 = 0.0;
            d2 = d;
        }

        lng_ad2 = gamln(a + d2);
        lng_bd2 = gamln(b + d2);
        lng_c = gamln(c);

        /* Evaluate F1. */
        if (ad < DBL_EPSILON) {
            /* d = 0 */
            F1 = 0.0;
        } else {
            double lng_ad = gamln(ad);
            double lng_ad1 = gamln(a + d1);
            double lng_bd1 = gamln(b + d1);

            /* Gamma functions in the denominator are ok.
             * Proceed with evaluation.
             */
            int i;
            double sum1 = 1.0;
            double term = 1.0;
            double ln_pre1_val = lng_ad + lng_c + d2 * ln_omx - lng_ad1 - lng_bd1;

            /* Do F1 sum. */
            for (i = 1; i < ad; i++) {
                int j = i - 1;
                term *= (a + d2 + j) * (b + d2 + j) / (1.0 + d2 + j) / i * (1.0 - x);
                sum1 += term;
            }

            F1 = exp(ln_pre1_val) * sum1;
        } /* end F1 evaluation */

        /* Evaluate F2. */
        if (1) { // Assume gamma functions are ok
            const int maxiter = 2000;
            double psi_1 = -M_E;
            double psi_1pd = psi(1.0 + ad);
            double psi_apd1 = psi(a + d1);
            double psi_bpd1 = psi(b + d1);

            double psi_val = psi_1 + psi_1pd - psi_apd1 - psi_bpd1 - ln_omx;
            double fact = 1.0;
            double sum2_val = psi_val;
            double ln_pre2_val = lng_c + d1 * ln_omx - lng_ad2 - lng_bd2;

            int j;

            /* Do F2 sum. */
            for (j = 1; j < maxiter; j++) {
                /* values for psi functions use recurrence; Abramowitz+Stegun 6.3.5 */
                double term1 = 1.0 / (double)j + 1.0 / (ad + j);
                double term2 = 1.0 / (a + d1 + j - 1.0) + 1.0 / (b + d1 + j - 1.0);
                double delta = 0.0;
                psi_val += term1 - term2;
                fact *= (a + d1 + j - 1.0) * (b + d1 + j - 1.0) / ((ad + j) * j) * (1.0 - x);
                delta = fact * psi_val;
                sum2_val += delta;
                if (fabs(delta) < DBL_EPSILON * fabs(sum2_val)) break;
            }

            F2 = exp(ln_pre2_val) * sum2_val;
        } else {
            /* Gamma functions in the denominator not ok.
             * So the F2 term is zero.
             */
            F2 = 0.0;
        } /* end F2 evaluation */

        sgn_2 = ((intd % 2 != 0) ? -1.0 : 1.0);
        return F1 + sgn_2 * F2;
    } else {
        /* d not an integer */

        double pre1, pre2;
        double sgn1, sgn2;
        double F1, F2;

        /* These gamma functions appear in the denominator, so we
         * catch their harmless domain errors and set the terms to zero.
         */
        double ln_g1ca = gamln(c - a);
        double ln_g1cb = gamln(c - b);
        double ln_g2a = gamln(a);
        double ln_g2b = gamln(b);

        double ln_gc = gamln(c);
        double ln_gd = gamln(d);
        double ln_gmd = gamln(-d);

        sgn1 = 1.0; // Assume positive signs
        sgn2 = 1.0; // Assume positive signs

        if (1) { // Assume gamma functions are ok
            double ln_pre1_val = ln_gc + ln_gd - ln_g1ca - ln_g1cb;
            double ln_pre2_val = ln_gc + ln_gmd - ln_g2a - ln_g2b + d * log(1.0 - x);

            pre1 = exp(ln_pre1_val) * sgn1;
            pre2 = exp(ln_pre2_val) * sgn2;
        } else {
            pre1 = 0.0;
            pre2 = 0.0;
        }

        F1 = hyperg_2F1_series(a, b, 1.0 - d, 1.0 - x);
        F2 = hyperg_2F1_series(c - a, c - b, 1.0 + d, 1.0 - x);

        return pre1 * F1 + pre2 * F2;
    }
}

double pow_omx(double x, double p) {
    double ln_omx;
    if (fabs(x) < sqrt(sqrt(DBL_EPSILON))) {
        ln_omx = -x * (1.0 + x * (1.0 / 2.0 + x * (1.0 / 3.0 + x / 4.0 + x * x / 5.0)));
    } else {
        ln_omx = log(1.0 - x);
    }
    return exp(p * ln_omx);
}

double sf_hyperg_2F1(double a, double b, double c, double x) {
    const double d = c - a - b;
    const double rinta = floor(a + 0.5);
    const double rintb = floor(b + 0.5);
    const double rintc = floor(c + 0.5);
    const int a_neg_integer = (a < 0.0 && fabs(a - rinta) < locEPS);
    const int b_neg_integer = (b < 0.0 && fabs(b - rintb) < locEPS);
    const int c_neg_integer = (c < 0.0 && fabs(c - rintc) < locEPS);

    /* Handle x == 1.0 RJM */
    if (fabs(x - 1.0) < locEPS && (c - a - b) > 0 && c != 0 && !c_neg_integer) {
        double lngamc = gamln(c);
        double lngamcab = gamln(c - a - b);
        double lngamca = gamln(c - a);
        double lngamcb = gamln(c - b);

        return exp(lngamc + lngamcab - lngamca - lngamcb);
    }

    if (x < -1.0 || 1.0 <= x) {
        return NAN;
    }

    if (c_neg_integer) {
        /* If c is a negative integer, then either a or b must be a
           negative integer of smaller magnitude than c to ensure
           cancellation of the series. */
        if (!(a_neg_integer && a > c + 0.1) && !(b_neg_integer && b > c + 0.1)) {
            return NAN;
        }
    }

    if (fabs(c - b) < locEPS || fabs(c - a) < locEPS) {
        return pow_omx(x, d);  /* (1-x)^(c-a-b) */
    }

    if (a >= 0.0 && b >= 0.0 && c >= 0.0 && x >= 0.0 && x < 0.995) {
        /* Series has all positive definite
         * terms and x is not close to 1.
         */
        return hyperg_2F1_series(a, b, c, x);
    }

    if (fabs(a) < 10.0 && fabs(b) < 10.0) {
        /* a and b are not too large, so we attempt
         * variations on the series summation.
         */
        if (a_neg_integer) {
            return hyperg_2F1_series(rinta, b, c, x);
        }
        if (b_neg_integer) {
            return hyperg_2F1_series(a, rintb, c, x);
        }

        if (x < -0.25) {
            return hyperg_2F1_luke(a, b, c, x);
        } else if (x < 0.5) {
            return hyperg_2F1_series(a, b, c, x);
        } else {
            if (fabs(c) > 10.0) {
                return hyperg_2F1_series(a, b, c, x);
            } else {
                return hyperg_2F1_reflect(a, b, c, x);
            }
        }
    } else {
        /* Either a or b or both large.
         * Introduce some new variables ap,bp so that bp is
         * the larger in magnitude.
         */
        double ap, bp;
        if (fabs(a) > fabs(b)) {
            bp = a;
            ap = b;
        } else {
            bp = b;
            ap = a;
        }

        if (x < 0.0) {
            // What the hell, maybe Luke will converge.
            return hyperg_2F1_luke(a, b, c, x);
        }

        if (DBL_MAX * (fabs(ap), 1.0) * fabs(bp) * fabs(x) < 2.0 * fabs(c)) {
            // If c is large enough or x is small enough,
            //we can attempt the series anyway.
            return hyperg_2F1_series(a, b, c, x);
        }

        /* We give up. */
        return NAN;
    }
}


// betaln_hyp function implementation
long double betaln_hyp(long double x, long double a, long double b) {
    if (x <= 0) {return -INFINITY;}  // ln(0) = -inf
    if (x >= 1) {return 0.0;}        // ln(1) = 0

    // Compute the hypergeometric function 1F1(a + b, a + 1, x)
    long double hyp2f1 = sf_hyperg_2F1(a + b, 1.0, a + 1, x);
    // Compute the logarithm of the beta function
    long double ln_beta = betaln(a, b);

    // Combine terms
    long double result = logl(hyp2f1) + a * logl(x) + b * logl(1 - x) - logl(a) - ln_beta;
    return result;
}

// logfdtrc function implementation
long double logfdtrc(long double x, int ia, int ib) {
    // Check for domain errors
    if (ia < 1 || ib < 1 || x < 0.0) {
        fprintf(stderr, "fdtrc domain error: ia < 1, ib < 1, or x < 0\n");
        exit(1);
    }

    long double a = ia;
    long double b = ib;
    long double w = b / (b + a * x);

    // Call betaln with corrected argument order
    return betaln_hyp(w, 0.5 * b, 0.5 * a);
}


// Main function
int main(int argc, char* argv[]) {
    // Check for command-line arguments
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <N> <R>\n", argv[0]);
        fprintf(stderr, "N - number of measurements\n");
        fprintf(stderr, "R = chi^2_0 / chi^2_{ref} \n");
        return 1;
    }

    // Parse command-line arguments
    int N = atoi(argv[1]);  // Number of measurements
    long double R = atof(argv[2]);  // Ratio of chi-squared values

    // Degrees of freedom
    int N1 = N - 2;  // Degrees of freedom for numerator
    int N2 = (2 * N) - 4;  // Degrees of freedom for denominator

    // Compute log-SF using logfdtrc implementation
    long double log_sf = logfdtrc(R, N1, N2);

    // Output the result
    printf("ln(p): %Lf\n", log_sf);

    return 0;
}
