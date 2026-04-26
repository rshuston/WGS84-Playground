#import <XCTest/XCTest.h>

#import <math.h>

#import "../WGS84/wgs84.h"


#define deg2rad(x)  ((x) * M_PI / 180.0)
#define rad2deg(x)  ((x) * 180.0 / M_PI)


@interface geodetic2ecef_Tests : XCTestCase

@end

@implementation geodetic2ecef_Tests

- (void)test_geodetic2ecef_Equator_GreenwichMeridian {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    wgs84_geodetic2ecef(0.0, 0.0, 123.0, &x, &y, &z);

    XCTAssertEqualWithAccuracy(wgs84_a + 123.0, x, 0.001);
    XCTAssertEqualWithAccuracy(0.0, y, 0.001);
    XCTAssertEqualWithAccuracy(0.0, z, 0.001);
}

- (void)test_geodetic2ecef_Equator_InternationalDateline {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    wgs84_geodetic2ecef(0.0, -M_PI, 123.0, &x, &y, &z);

    XCTAssertEqualWithAccuracy(-(wgs84_a + 123.0), x, 0.001);
    XCTAssertEqualWithAccuracy(0.0, y, 0.001);
    XCTAssertEqualWithAccuracy(0.0, z, 0.001);
}

- (void)test_geodetic2ecef_NorthPole {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    wgs84_geodetic2ecef(M_PI_2, 0.0, 123.0, &x, &y, &z);

    XCTAssertEqualWithAccuracy(0.0, x, 0.001);
    XCTAssertEqualWithAccuracy(0.0, y, 0.001);
    XCTAssertEqualWithAccuracy(wgs84_b + 123.0, z, 0.001);
}

- (void)test_geodetic2ecef_SouthPole {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    wgs84_geodetic2ecef(-M_PI_2, 0.0, 123.0, &x, &y, &z);

    XCTAssertEqualWithAccuracy(0.0, x, 0.001);
    XCTAssertEqualWithAccuracy(0.0, y, 0.001);
    XCTAssertEqualWithAccuracy(-(wgs84_b + 123.0), z, 0.001);
}

- (void)test_geodetic2ecef_SydneyAustralia {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    // From mapcoordinates.net
    double latitude = deg2rad(-33.869844);
    double longitude = deg2rad(151.208285);
    double height = 87.0;  // meters
    
    wgs84_geodetic2ecef(latitude, longitude, height, &x, &y, &z);

    XCTAssertEqualWithAccuracy(-4646012.8, x, 0.1);
    XCTAssertEqualWithAccuracy(2553292.4, y, 0.1);
    XCTAssertEqualWithAccuracy(-3534517.0, z, 0.1);
}

- (void)test_geodetic2ecef_SomewhereSomeplace {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    // Numerical values obtained from:
    // Transformation of Cartesian to Geodetic Coordinates without Iterations
    // Rey-Jer You
    // Journal of Surveying Engineering, Vol. 126, No. 1, February, 2000
    double latitude = deg2rad(45.0);
    double longitude = deg2rad(120.0);
    double height = 10000.0;  // meters

    wgs84_geodetic2ecef(latitude, longitude, height, &x, &y, &z);

    XCTAssertEqualWithAccuracy(-2262330.973, x, 0.001);
    XCTAssertEqualWithAccuracy(3918472.189, y, 0.001);
    XCTAssertEqualWithAccuracy(4494419.477, z, 0.001);
}

@end
