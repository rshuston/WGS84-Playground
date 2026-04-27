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
const double phi_step_deg = 0.25;

const double h_min_m = -20000.0;    // Mariana trench is -10984 m
const double h_max_m = 100000.0;    // Airliner altitudes are 11600 m, GPS orbit height is 20180 km
const double h_step_m = 20000.0;


int main(int argc, const char * argv[]) {
    unsigned N_heights = static_cast<unsigned>((h_max_m - h_min_m) / h_step_m + 1);
    unsigned N_latitudes = static_cast<unsigned>((phi_max_deg - phi_min_deg) / phi_step_deg + 1);
    unsigned N_points = N_heights * N_latitudes;
    printf("Number of heights:   %u\n", N_heights);
    printf("Number of latitudes: %u\n", N_latitudes);
    printf("Number of points:    %u\n", N_points);
    
    double *EcefPoints = (double *)std::malloc(N_points * 3 * sizeof(double));
    double *GeodeticPoints = (double *)std::malloc(N_points * 3 * sizeof(double));

    for (int i = 0; i < N_heights; i++)
    {
        double h_m = h_min_m + i * h_step_m;
//        printf("h = %g m\n", h_m);
        for (int j = 0; j < N_latitudes; j++)
        {
            double phi_deg = phi_min_deg + j * phi_step_deg;
//            printf("  phi = %g deg\n", phi_deg);
            int k = i * N_latitudes + j;
            wgs84_geodetic2ecef(deg2rad(phi_deg), deg2rad(lambda_reference_deg), h_m,
                                &EcefPoints[3 * k], &EcefPoints[3 * k + 1], &EcefPoints[3 * k + 2]);
        }
    }
//    for (int k = 0; k < N_points; k++)
//    {
//        int k3 = 3 * k;
//        printf("p[%02d] = (%.1f, %.1f, %.1f)\n", k, EcefPoints[k3], EcefPoints[k3+1], EcefPoints[k3+2]);
//    }

    std::chrono::steady_clock::time_point start, end;
    std::chrono::duration<double, std::milli> duration;
    double duration_ms;
    
    start = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < N_points; k++)
    {
        int k3 = 3 * k;
        wgs84_ecef2geodetic_iter(EcefPoints[k3], EcefPoints[k3 + 1], EcefPoints[k3 + 2],
                                 &GeodeticPoints[k3], &GeodeticPoints[k3 + 1], &GeodeticPoints[k3 + 2]);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
    duration_ms = duration.count();
    printf("Iterative time elapsed: %f ms\n", duration_ms);
    
    start = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < N_points; k++)
    {
        int k3 = 3 * k;
        wgs84_ecef2geodetic_heikkinen(EcefPoints[k3], EcefPoints[k3 + 1], EcefPoints[k3 + 2],
                                      &GeodeticPoints[k3], &GeodeticPoints[k3 + 1], &GeodeticPoints[k3 + 2]);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
    duration_ms = duration.count();
    printf("Heikkinen time elapsed: %f ms\n", duration_ms);
    
    start = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < N_points; k++)
    {
        int k3 = 3 * k;
        wgs84_ecef2geodetic_olson(EcefPoints[k3], EcefPoints[k3 + 1], EcefPoints[k3 + 2],
                                  &GeodeticPoints[k3], &GeodeticPoints[k3 + 1], &GeodeticPoints[k3 + 2]);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
    duration_ms = duration.count();
    printf("Olson time elapsed:     %f ms\n", duration_ms);

    std::free(EcefPoints);
    EcefPoints = nullptr;
    std::free(GeodeticPoints);
    GeodeticPoints = nullptr;

    return EXIT_SUCCESS;
}
