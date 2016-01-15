#include "morton.h"


int64_t naiveEncodeGeohash(double lat, double lng, uint8_t bitLen) {
    assert(bitLen <= 64);
    double minLat = llRange.ymin,  maxLat = llRange.ymax;
    double minLng = llRange.xmin, maxLng = llRange.xmax;
    double midpoint = 0;
    int64_t result = 0;
    int8_t bit=0;
    for (int i = 0; i < bitLen; ++i)
    {
        if (0 == (i&1)) // even bit: bisect longitude
        {
            midpoint = (minLng + maxLng) / 2;
            if (lng < midpoint)
            {
                bit = 0; // push a zero bit
                maxLng = midpoint;   // shrink range downwards
            }
            else
            {
                bit = 1; // push a one bit
                minLng = midpoint;  // shrink range upwards
            }
        }
        else // odd bit: bisect latitude
        {
            midpoint = (minLat + maxLat) / 2;
            if (lat < midpoint)
            {
                bit = 0;                 // push a zero bit
                maxLat = midpoint;            // shrink range downwards
            }
            else
            {
                bit = 1;     // push a one bit
                minLat = midpoint;            // shrink range upwards
            }
        }
        result <<= 1;
        result += bit;
  }
  return result;
}
Envelope naiveDecodeGeohash(const int64_t code, const uint8_t bitLen)
{
    assert(bitLen <= 64);
    double minLat = llRange.ymin,  maxLat = llRange.ymax;
    double minLng = llRange.xmin, maxLng = llRange.xmax;
    Envelope env(minLng,minLat,maxLng,maxLat);
    int8_t pos=bitLen-1;
    for(int i=0; i < bitLen; ++i)
    {
        if(0 == (i&1))
        {
            if((code & (1ll<<(pos-i))) == 0) //lng < midpoint
                env.xmax = (env.xmax + env.xmin) / 2;
            else
                env.xmin = (env.xmax + env.xmin) / 2;
        }
        else
        {
            if((code & (1ll<<(pos-i))) == 0) //lat < midpoint
                env.ymax = (env.ymax + env.ymin)/2;
            else
                env.ymin = (env.ymax + env.ymin)/2;
        }
    }
    return env;
}
