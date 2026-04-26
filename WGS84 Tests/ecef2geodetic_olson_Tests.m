#import <XCTest/XCTest.h>

#import <math.h>

#import "../WGS84/wgs84.h"


#define deg2rad(x)  ((x) * M_PI / 180.0)
#define rad2deg(x)  ((x) * 180.0 / M_PI)


@interface ecef2geodetic_olson_Tests : XCTestCase

@end

@implementation ecef2geodetic_olson_Tests

- (void)test_ecef2geodetic_Equator_GreenwichMeridian {
    double latitude = 2.0 * M_PI;
    double longitude = 2.0 * M_PI;
    double height = 2.0 * wgs84_a;

    wgs84_ecef2geodetic_olson(wgs84_a + 123.0, 0.0, 0.0, &latitude, &longitude, &height);
    
    XCTAssertEqualWithAccuracy(0.0, latitude, 0.001);
    XCTAssertEqualWithAccuracy(0.0, longitude, 0.001);
    XCTAssertEqualWithAccuracy(123.0, height, 0.001);
}

- (void)test_ecef2geodetic_Equator_InternationalDateLine {
    double latitude = 2.0 * M_PI;
    double longitude = 2.0 * M_PI;
    double height = 2.0 * wgs84_a;

    wgs84_ecef2geodetic_olson(-(wgs84_a + 123.0), 0.0, 0.0, &latitude, &longitude, &height);
    
    XCTAssertEqualWithAccuracy(0.0, latitude, 0.001);
    XCTAssertEqualWithAccuracy(-M_PI, longitude, 0.001);
    XCTAssertEqualWithAccuracy(123.0, height, 0.001);
}

- (void)test_ecef2geodetic_NorthPole {
    double latitude = 2.0 * M_PI;
    double longitude = 2.0 * M_PI;
    double height = 2.0 * wgs84_a;

    wgs84_ecef2geodetic_olson(0.0, 0.0, wgs84_b + 123.0, &latitude, &longitude, &height);

    XCTAssertEqualWithAccuracy(M_PI_2, latitude, 0.001);
    XCTAssertEqualWithAccuracy(0.0, longitude, 0.001);
    XCTAssertEqualWithAccuracy(123.0, height, 0.001);
}

- (void)test_ecef2geodetic_SouthPole {
    double latitude = 2.0 * M_PI;
    double longitude = 2.0 * M_PI;
    double height = 2.0 * wgs84_a;

    wgs84_ecef2geodetic_olson(0.0, 0.0, -(wgs84_b + 123.0), &latitude, &longitude, &height);

    XCTAssertEqualWithAccuracy(-M_PI_2, latitude, 0.001);
    XCTAssertEqualWithAccuracy(0.0, longitude, 0.001);
    XCTAssertEqualWithAccuracy(123.0, height, 0.001);
}

- (void)test_ecef2geodetic_SydneyAustralia {
    double latitude = 2.0 * M_PI;
    double longitude = 2.0 * M_PI;
    double height = 2.0 * wgs84_a;

    // From mapcoordinates.net
    wgs84_ecef2geodetic_olson(-4646012.8, 2553292.4, -3534517.0, &latitude, &longitude, &height);

    XCTAssertEqualWithAccuracy(-33.869844, rad2deg(latitude), 1e-6);
    XCTAssertEqualWithAccuracy(151.208285, rad2deg(longitude), 1e-6);
    XCTAssertEqualWithAccuracy(87.0, height, 1);
}

- (void)test_ecef2geodetic_SomewhereSomeplace {
    double latitude = 2.0 * M_PI;
    double longitude = 2.0 * M_PI;
    double height = 2.0 * wgs84_a;

    // Numerical values obtained from:
    // Transformation of Cartesian to Geodetic Coordinates without Iterations
    // Rey-Jer You
    // Journal of Surveying Engineering, Vol. 126, No. 1, February, 2000

    wgs84_ecef2geodetic_olson(-2262330.973, 3918472.189, 4494419.477, &latitude, &longitude, &height);

    XCTAssertEqualWithAccuracy(45.0, rad2deg(latitude), 1e-6);
    XCTAssertEqualWithAccuracy(120.0, rad2deg(longitude), 1e-6);
    XCTAssertEqualWithAccuracy(10000.0, height, 1);
}

@end
