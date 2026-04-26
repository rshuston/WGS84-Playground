#ifndef wgs84_h
#define wgs84_h

#ifdef __cplusplus
extern "C" {
#endif

extern const double wgs84_a;            /* Semi-major axis, meters */
extern const double wgs84_inverse_f;
extern const double wgs84_b;            /* Semi-minor axis, meters */

/*
 * All lengths are in meters.
 * All angles are in radians.
 */

int wgs84_ecef2geodetic_iter(double x, double y, double z, double *phi, double *lambda, double *h);
void wgs84_ecef2geodetic_heikkinen(double xu, double yu, double zu, double *phi, double *lambda, double *h);
void wgs84_ecef2geodetic_olson(double x, double y, double z, double *lat, double *lon, double *ht);

void wgs84_geodetic2ecef(double phi, double lambda, double h, double *x, double *y, double *z);

#ifdef __cplusplus
}
#endif

#endif /* wgs84_h */
