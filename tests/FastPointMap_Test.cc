#define CATCH_CONFIG_MAIN
#include "Catch-master/single_include/catch.hpp"
#include "PhasespacePoint.hh"
#include "FastPointMap.hh"



TEST_CASE("FastPointMap Tests"){

    PhasespacePoint p1;
    PhasespacePoint p2;
    FastPointMap fpm1(&p1, 1);
    FastPointMap fpm2(&p2, 2);
    FastPointMap comp(NULL, 0);

    REQUIRE(comp(fpm1, fpm2) == true);
    REQUIRE(comp(fpm2, fpm1) == false);
    REQUIRE(comp(fpm1, fpm1) == false);
}