#include "wgs84.h"

#include <math.h>


/*
 * WGS84 specifications:
 *   Earth ellipsoid semi-major axis, a = 6378137.0 m
 *   Earth inverse ellipsoid flattening, 1/f = 298.257223563
 *
 * From these, all other parameters can be derived:
 *   Earth ellipsoid semi-minor axis, b = a * (1-f)
 *   First eccentricity squared, e^2 = f * (2-f)
 *   Second eccentricity squared, (e')^2 = e^2 / ( (1-f)^2 )
 */

const double wgs84_a = 6378137.0;
const double wgs84_inverse_f = 298.257223563;
const double wgs84_b = 6.356752314245179e6;


/*
 * Private constants and functions
 */

static const double wgs84_e2 = 6.694379990141316e-03;           /* 1st ecc. squared */
static const double wgs84_ep2 = 6.739496742276434e-03;          /* 2nd ecc. squared */
static const double wgs84_one_minus_e2 = 9.933056200098587e-01;
static const double wgs84_one_minus_f = 9.966471893352525e-01;


static double sqr(double x) { return x * x; }


/*
 * Iterative method, typically converges in less than 5 iterations for terrestrial locaions
 *
 * Hofmann-Wellenhof, B., Lichtenegger, H., and Collins, J.
 * Global Positioning System, Theory and Practice, 3rd ed.
 * New York, Springer-Verlag Wien, 1994
 */

int wgs84_ecef2geodetic_iter(double x, double y, double z, double *phi, double *lambda, double *h)
{
    int k = 0;
    double r, phi_k, cos_phi_k, sin_phi_k, R_k, h_k, phi_kp1;

    *lambda = atan2(y, x);
    if (*lambda >= M_PI)
    {
        *lambda -= 2.0 * M_PI; /* atan2() can return +M_PI, so we need to normalize to -M_PI */
    }

    r = sqrt(x * x + y * y);

    phi_kp1 = atan2(z, r * wgs84_one_minus_e2);
    for (k = 1; k <= 10; k++) /* typically 5 iterations gives centimeter accuracy for terrestrial locations */
    {
        phi_k = phi_kp1;
        cos_phi_k = cos(phi_k);
        sin_phi_k = sin(phi_k);
        R_k = wgs84_a / sqrt(1.0 - wgs84_e2 * sin_phi_k * sin_phi_k);
        if (fabs(phi_k) <= M_PI_4) { /* Guard against /0 by using best form of h */
            h_k = r / cos_phi_k - R_k;
        }
        else
        {
            h_k = z / sin_phi_k - wgs84_one_minus_e2 * R_k;
        }

        phi_kp1 = atan2(z, r * (1.0 - wgs84_e2 * R_k / (R_k + h_k)));
        if (fabs(phi_kp1 - phi_k) < 1e-10) break;
    }

    phi_k = phi_kp1;
    cos_phi_k = cos(phi_k);
    sin_phi_k = sin(phi_k);
    R_k = wgs84_a / sqrt(1.0 - wgs84_e2 * sin_phi_k * sin_phi_k);
    if (fabs(phi_k) <= M_PI_4) { /* Guard against /0 by using best form of h */
        h_k = r / cos_phi_k - R_k;
    }
    else
    {
        h_k = z / sin_phi_k - wgs84_one_minus_e2 * R_k;
    }

    *phi = phi_k;
    *h = h_k;

    return k;
}


/*
 * Bowring's method, single iteration is sufficiently accurate for most general applications
 *
 * Hofmann-Wellenhof, B., Lichtenegger, H., and Collins, J.
 * Global Positioning System, Theory and Practice, 3rd ed.
 * New York, Springer-Verlag Wien, 1994
 */

void wgs84_ecef2geodetic_bowring(double x, double y, double z, double *phi, double *lambda, double *h)
{
    double r, numer, denom, norm, cos_theta, sin_theta, cos_phi, sin_phi, R;

    *lambda = atan2(y, x);
    if (*lambda >= M_PI)
    {
        *lambda -= 2.0 * M_PI; /* atan2() can return +M_PI, so we need to normalize to -M_PI */
    }
    
    r = sqrt(x * x + y * y);

    /*
     * theta = atan2(z * wgs84_a, r * wgs84_b);
     * ... or
     * theta = atan2(z, r * wgs84_one_minus_f);
     * cos_theta = cos(theta);
     * sin_theta = sin(theta);
     */
    numer = z;
    denom = r * wgs84_one_minus_f;
    norm = sqrt(numer * numer + denom * denom);
    sin_theta = numer / norm;
    cos_theta = denom / norm;
    
    /*
     * *phi = atan2(z + wgs84_ep2 * wgs84_b * sin_theta * sin_theta * sin_theta,
     *              r - wgs84_e2 * wgs84_a * cos_theta * cos_theta * cos_theta);
     * cos_phi = cos(*phi);
     * sin_phi = sin(*phi);
     */
    numer = z + wgs84_ep2 * wgs84_b * sin_theta * sin_theta * sin_theta;
    denom = r - wgs84_e2 * wgs84_a * cos_theta * cos_theta * cos_theta;
    norm = sqrt(numer * numer + denom * denom);
    sin_phi = numer / norm;
    cos_phi = denom / norm;
    *phi = atan2(sin_phi, cos_phi);
    
    R = wgs84_a / sqrt(1.0 - wgs84_e2 * sin_phi * sin_phi);
    if (fabs(*phi) <= M_PI_4) { /* Guard against /0 by using best form of h */
        *h = r / cos_phi - R;
    }
    else
    {
        *h = z / sin_phi - wgs84_one_minus_e2 * R;
    }
}


/*
 * Heikkinen's method: Optimized to compute fewest terms, but requires cbrt() (use pow() if you must)
 *
 * Elliot D. Kaplan, Ed.
 * Understanding GPS, Principles and Applications
 * Boston, MA, Academic Press, Inc., 1996
 */

void wgs84_ecef2geodetic_heikkinen(double xu, double yu, double zu, double *phi, double *lambda, double *h)
{
    static const double a2 = 4.068063159076900e13;      /* a^2 */
    static const double b2 = 4.040829998466145e13;      /* b^2 */
    static const double e4 = 4.481472345240445e-05;     /* e^4 */
    static const double a2_e4 = 1.823091254609461e9;    /* e^2 * (a^2 - b^2) == a^2 * e^4 */
    
    double zu2, one_minus_e2_x_zu2,
           r2, r,
           F, G, G2,
           c, s, k, P, Q, one_plus_Q,
           w, D, r0, r_minus_e2r0_2,
           U, V, aV, z0;
    
    zu2 = zu * zu;
    one_minus_e2_x_zu2 = wgs84_one_minus_e2 * zu2;

    r2 = xu * xu + yu * yu;
    r = sqrt(r2);
    F = 54.0 * b2 * zu2;
    G = r2 + one_minus_e2_x_zu2 - a2_e4;
    G2 = G * G;
    c = e4 * F * r2 / (G2 * G);
    s = cbrt(1.0 + c + sqrt(c * (c + 2.0)));
    /* s = pow(1.0 + c + sqrt(c * (c + 2.0)), 1.0/3.0); */
    k = s + 1.0/s + 1.0;
    P = F / (3.0 * k * k * G2);
    Q = sqrt(1.0 + 2.0 * e4 * P);
    one_plus_Q = 1.0 + Q;
    w = -P * wgs84_e2 * r / one_plus_Q;
    D = 0.5 * a2 * (1.0 + 1.0/Q) - P * one_minus_e2_x_zu2 / (Q * one_plus_Q) - 0.5 * P * r2;
    r0 = D >= 0.0 ? w + sqrt(D) : w; /* D can be numerically small negative at the poles */
    r_minus_e2r0_2 = sqr(r - wgs84_e2 * r0);
    U = sqrt(r_minus_e2r0_2 + zu2);
    V = sqrt(r_minus_e2r0_2 + one_minus_e2_x_zu2);
    aV = wgs84_a * V;
    z0 = (b2 * zu) / aV;
    *h = U * (1.0 - b2 / aV);
    *phi = atan2((zu + wgs84_ep2 * z0), r);
    *lambda = atan2(yu, xu);
    
    if (*lambda >= M_PI)
    {
        *lambda -= 2.0 * M_PI; /* atan2() can return +M_PI, so we need to normalize to -M_PI */
    }
}


/*
 * Olson's method: Same accuracy as Heikkinen but requires less computations, does not need
 * a cube root operation, and is 63% quicker.
 *
 * Code is taken directly from the IEEE AES article referenced below, but with constants
 * made static const so they aren't computed at runtime, and with a finishing change to
 * normalize longitude +M_PI to -M_PI.
 *
 * Donald K. Olson
 * Converting Earth-centered, Earth-fixed coordinates to geodetic coordinates
 * IEEE Trans. on Aerospace and Electronic Systems, Vol. 32, Issue 1, January 1996
 */

void wgs84_ecef2geodetic_olson(double x, double y, double z, double *lat, double *lon, double *ht)
{
    static const double a = 6378137.0; /*wgs-84*/
    static const double e2 = 6.6943799901377997e-3;
    static const double a1 = 4.2697672707157535e+4;
    static const double a2 = 1.8230912546075455e+9;
    static const double a3 = 1.4291722289812413e+2;
    static const double a4 = 4.5577281365188637e+9;
    static const double a5 = 4.2840589930055659e+4;
    static const double a6 = 9.9330562000986220e-1;
    /* a1 = a*e2, a2 = a1*a1, a3 = a1*e2/2, */
    /* a4 = (5/2)*a2, a5 = a1+a3, a6 = 1-e2 */
    double zp,w2,w,z2,r2,r,s2,c2,s,c,ss;
    double g,rg,rf,u,v,m,f,p;
    zp = fabs(z);
    w2 = x*x+y*y;
    w = sqrt(w2);
    z2 = z*z;
    r2 = w2+z2;
    r = sqrt(r2);
    if (r < 100000.)
    {
        *lat = 0.;
        *lon = 0.;
        *ht = -1.e7;
        return;
    }
    *lon = atan2(y,x);
    s2 = z2/r2;
    c2 = w2/r2;
    u = a2/r;
    v = a3-a4/r;
    if (c2 > .3)
    {
        s = (zp/r)*(1.+c2*(a1+u+s2*v)/r);
        *lat = asin(s);
        ss = s*s;
        c = sqrt(1.-ss);
    }
    else
    {
        c = (w/r)*(1.-s2*(a5-u-c2*v)/r);
        *lat = acos(c);
        ss = 1.-c*c;
        s = sqrt(ss);
    }
    g = 1.-e2*ss;
    rg = a/sqrt(g);
    rf = a6*rg;
    u = w-rg*c;
    v = zp-rf*s;
    f = c*u+s*v;
    m = c*v-s*u;
    p = m/(rf/g+f);
    *lat = *lat+p;
    *ht = f+m*p/2.;
    if (z < 0.)
        *lat = -*lat;
    
    if (*lon >= M_PI)
    {
        *lon -= 2.0 * M_PI; /* atan2() can return +M_PI, so we need to normalize to -M_PI */
    }
}


/*
 * Standard method to convert geodetic to ECEF
 *
 * Elliot D. Kaplan, Ed.
 * Understanding GPS, Principles and Applications
 * Boston, MA, Academic Press, Inc., 1996
 *
 * Hofmann-Wellenhof, B., Lichtenegger, H., and Collins, J.
 * Global Positioning System, Theory and Practice, 3rd ed.
 * New York, Springer-Verlag Wien, 1994
 */

void wgs84_geodetic2ecef(double phi, double lambda, double h, double *x, double *y, double *z)
{
    double cos_phi, sin_phi, cos_lambda, sin_lambda, R;
    
    cos_phi = cos(phi);
    sin_phi = sin(phi);
    cos_lambda = cos(lambda);
    sin_lambda = sin(lambda);

    R = wgs84_a / sqrt(1.0 - wgs84_e2 * sin_phi * sin_phi);

    *x = (R + h) * cos_phi * cos_lambda;
    *y = (R + h) * cos_phi * sin_lambda;
    *z = (R * wgs84_one_minus_e2 + h) * sin_phi;
}
