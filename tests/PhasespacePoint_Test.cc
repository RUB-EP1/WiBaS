#include "Catch-master/single_include/catch.hpp"
#include "PhasespacePoint.hh"



TEST_CASE("PhasespacePoint Tests"){

    PhasespacePoint point;
    point.SetCoordinate("c3", 30);
    point.SetCoordinate("c2", 20);
    point.SetCoordinate("c1", 10);

    PhasespaceCoord c1(0, 1, false);
    PhasespaceCoord c2(1, 2, false);
    PhasespaceCoord c3(2, 3, false);

    std::map< std::string, PhasespaceCoord > coordNameMap;
    coordNameMap.insert( std::pair< std::string, PhasespaceCoord >("c1", c1));
    coordNameMap.insert( std::pair< std::string, PhasespaceCoord >("c2", c2));
    coordNameMap.insert( std::pair< std::string, PhasespaceCoord >("c3", c3));

    SECTION( "Check coordinate value mapping" ) {
        REQUIRE(point.coordValueMap.size() == 3);
        REQUIRE(point.coordValueMap.at("c2") == 20);
    }

    SECTION( "Check arranged value vector"){
        point.ArrangeCoordinates(coordNameMap);
        REQUIRE(point.coordValueVector.size() == 3);
        REQUIRE(point.GetCoordValue(0) == 10);
        REQUIRE(point.GetCoordValue(2) == 30);
    }
}
