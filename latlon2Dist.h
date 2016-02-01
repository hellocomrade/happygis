#ifndef LATLON2DIST_H
#define LATLON2DIST_H
#include <cmath>
//PI/180
const double RadianPerDegree = 0.01745329251994329576923690768489;
/*
 * https://en.wikipedia.org/wiki/Earth_radius
 *
 * and WGS84
 * GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]
 *
 * Quadratic Mean = sqrt((a^2+b^2)/2)
 *
 * Given a = 6378137 and b = a * (1 - 1/298.257223563)
 */
const double QuadraticMeanRadiusOfEarth =  6367453.6344937783740582933718864;

double lonlatDistHaversine(double lon1, double lat1, double lon2, double lat2);
#endif // LATLON2DIST_H
