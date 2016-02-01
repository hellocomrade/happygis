#include "latlon2Dist.h"

/*
 *http://www.movable-type.co.uk/scripts/latlong.html
 *
 * Haversine Great-Circle distance
 *
 * https://github.com/joto/osmium/blob/master/include/osmium/geometry/haversine.hpp
 */
double lonlatDistHaversine(double lon1, double lat1, double lon2, double lat2)
{
    double lonArc = (lon2 - lon1) * RadianPerDegree;
    double latArc = (lat2 - lat1) * RadianPerDegree;
    double lonh = sin(lonArc * 0.5);
    double lath = sin(latArc * 0.5);
    return 2.0 * QuadraticMeanRadiusOfEarth * asin(sqrt(lath * lath + cos(lat1 * RadianPerDegree) * cos(lat2 * RadianPerDegree) * lonh * lonh));
}

