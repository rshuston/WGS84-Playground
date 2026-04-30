#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "../WGS84/wgs84.h"


#define deg2rad(x)  ((x) * M_PI / 180.0)
#define rad2deg(x)  ((x) * 180.0 / M_PI)


const double lambda_reference_deg = 45.0;

const double phi_min_deg = -90.0;
const double phi_max_deg = 90.0;
const double phi_step_deg = 0.1;

const double h_min_m = -20000.0;    // Mariana trench is -10984 m
const double h_max_m = 100000.0;    // Airliner altitudes are 11600 m, GPS orbit height is 20180 km
const double h_step_m = 1000.0;


extern "C" {
    typedef int (*c_convfunc_t)(double, double, double, double*, double*, double*);
}


double evaluateConversionFunction(c_convfunc_t conv_f, unsigned numPoints, double ecef[], double geodetic[]);
double computePeakEcefError(unsigned numPoints, double ecef[], double geodetic[]);


int main(int argc, const char * argv[]) {
    unsigned N_heights = static_cast<unsigned>((h_max_m - h_min_m) / h_step_m + 1);
    unsigned N_latitudes = static_cast<unsigned>((phi_max_deg - phi_min_deg) / phi_step_deg + 1);
    unsigned N_points = N_heights * N_latitudes;
    printf("Number of heights:   %u\n", N_heights);
    printf("Number of latitudes: %u\n", N_latitudes);
    printf("Number of points:    %u\n", N_points);
    
    double *ecef_0 = (double *)std::malloc(N_points * 3 * sizeof(double));
    double *geodetic = (double *)std::malloc(N_points * 3 * sizeof(double));

    for (unsigned i = 0; i < N_heights; i++)
    {
        double h_m = h_min_m + i * h_step_m;
        for (unsigned j = 0; j < N_latitudes; j++)
        {
            double phi_deg = phi_min_deg + j * phi_step_deg;
            unsigned k = i * N_latitudes + j;
            wgs84_geodetic2ecef(deg2rad(phi_deg), deg2rad(lambda_reference_deg), h_m,
                                &ecef_0[3 * k], &ecef_0[3 * k + 1], &ecef_0[3 * k + 2]);
        }
    }
    
    double duration_ms, peakEcefError_m;
    
    duration_ms = evaluateConversionFunction(wgs84_ecef2geodetic_iter,
                                             N_points, ecef_0, geodetic);
    peakEcefError_m = computePeakEcefError(N_points, ecef_0, geodetic);
    printf("Iterative time: %f ms, peak error: %e m\n", duration_ms, peakEcefError_m);
    
    duration_ms = evaluateConversionFunction(wgs84_ecef2geodetic_bowring,
                                             N_points, ecef_0, geodetic);
    peakEcefError_m = computePeakEcefError(N_points, ecef_0, geodetic);
    printf("Bowring time:   %f ms, peak error: %e m\n", duration_ms, peakEcefError_m);
    
    duration_ms = evaluateConversionFunction(wgs84_ecef2geodetic_heikkinen,
                                             N_points, ecef_0, geodetic);
    peakEcefError_m = computePeakEcefError(N_points, ecef_0, geodetic);
    printf("Heikkinen time: %f ms, peak error: %e m\n", duration_ms, peakEcefError_m);
    
    duration_ms = evaluateConversionFunction(wgs84_ecef2geodetic_olson,
                                             N_points, ecef_0, geodetic);
    peakEcefError_m = computePeakEcefError(N_points, ecef_0, geodetic);
    printf("Olson time:     %f ms, peak error: %e m\n", duration_ms, peakEcefError_m);

    std::free(ecef_0);
    ecef_0 = nullptr;
    std::free(geodetic);
    geodetic = nullptr;

    return EXIT_SUCCESS;
}


double evaluateConversionFunction(c_convfunc_t conv_f, unsigned numPoints, double ecef[], double geodetic[])
{
    std::chrono::steady_clock::time_point start, end;
    std::chrono::duration<double, std::milli> duration;
    double duration_ms;

    start = std::chrono::high_resolution_clock::now();
    for (unsigned k = 0; k < numPoints; k++)
    {
        unsigned k3 = 3 * k;
        conv_f(ecef[k3], ecef[k3 + 1], ecef[k3 + 2],
               &geodetic[k3], &geodetic[k3 + 1], &geodetic[k3 + 2]);
    }
    end = std::chrono::high_resolution_clock::now();
    
    duration = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
    duration_ms = duration.count();
    
    return duration_ms;
}


double computePeakEcefError(unsigned numPoints, double ecef[], double geodetic[])
{
    double peakError = -1.0;
    
    for (unsigned k = 0; k < numPoints; k++)
    {
        unsigned k3 = 3 * k;
        double x, y, z;
        wgs84_geodetic2ecef(geodetic[k3], geodetic[k3 + 1], geodetic[k3 + 2],
                            &x, &y, &z);
        double d0 = sqrt(ecef[k3] * ecef[k3] + ecef[k3 + 1] * ecef[k3 + 1] + ecef[k3 + 2] * ecef[k3 + 2]);
        double d1 = sqrt(x * x + y * y + z * z);
        double error = fabs(d1 - d0);
        if (error > peakError)
        {
            peakError = error;
        }
    }

    return peakError;
}
