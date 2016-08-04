#include <map>
#include "Catch-master/single_include/catch.hpp"
#include "PhasespacePointCloud.hh"
#include "PhasespaceCoord.hh"
#include "PhasespacePoint.hh"



TEST_CASE("PhasespacePointCloud distance calculation"){

    PhasespacePointCloud cloud;
    cloud.RegisterPhasespaceCoord("c1", 10, false);
    cloud.RegisterPhasespaceCoord("c2", 15, false);

    PhasespacePoint point1;
    point1.SetCoordinate("c1", 3);
    point1.SetCoordinate("c2", 5);

    PhasespacePoint point2;
    point2.SetCoordinate("c1", 8);
    point2.SetCoordinate("c2", 2);

    SECTION("Check exception if unarranged coordinates are used"){
        REQUIRE_THROWS(cloud.CalcPhasespaceDistance(&point1, &point2));
    }

    cloud.ArrangePointCoordinates(point1);
    cloud.ArrangePointCoordinates(point2);

    // Distance calculation should yield:
    // sqrt( (3 - 8)^2 / (10. * 10.) + (5 - 2)^2 / (15. * 15.)) = 0.538516480713
    double refResult = 0.538516480713;
    double distance = cloud.CalcPhasespaceDistance(&point1, &point2);
    REQUIRE(distance == Approx(refResult).epsilon(1E-6));

    distance = cloud.CalcPhasespaceDistance(&point2, &point1);
    REQUIRE(distance == Approx(refResult).epsilon(1E-6));

    distance = cloud.CalcPhasespaceDistance(&point1, &point1);
    REQUIRE(distance == Approx(0).epsilon(1E-10));
}



TEST_CASE("PhasespacePointCloud distance calculation 2Pi circular coordinate"){
    PhasespacePointCloud cloud;
    cloud.RegisterPhasespaceCoord("c1", PhasespacePointCloud::Pi, PhasespacePointCloud::IS_2PI_CIRCULAR);

    PhasespacePoint point1;
    point1.SetCoordinate("c1", 0.1*PhasespacePointCloud::Pi);

    PhasespacePoint point2;
    point2.SetCoordinate("c1", 1.9 * PhasespacePointCloud::Pi);

    cloud.ArrangePointCoordinates(point1);
    cloud.ArrangePointCoordinates(point2);

    // Distance should yield 0.2, not 1.8
    double distance = cloud.CalcPhasespaceDistance(&point1, &point2);
    REQUIRE(distance == Approx(0.2).epsilon(1E-6));
}