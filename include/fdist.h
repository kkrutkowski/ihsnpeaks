#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>

#include "hyperg_2F1.c"

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
