#include "morton.h"
/*
 * http://devblogs.nvidia.com/parallelforall/thinking-parallel-part-iii-tree-construction-gpu/
 * http://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
 * https://www.factual.com/blog/how-geohashes-work
*/
uint64_t naiveEncodeGeohash(double lng, double lat, uint8_t bitLen, Envelope bbox) {
    assert(bitLen <= 64);
    double minLat = bbox.ymin,  maxLat = bbox.ymax;
    double minLng = bbox.xmin, maxLng = bbox.xmax;
    double midpoint = 0;
    uint64_t result = 0;
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
Envelope naiveDecodeGeohash(const int64_t code, const uint8_t bitLen, Envelope bbox)
{
    assert(bitLen <= 64);
    double minLat = bbox.ymin,  maxLat = bbox.ymax;
    double minLng = bbox.xmin, maxLng = bbox.xmax;
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
static uint64_t naiveMoveByOne(uint64_t val, const uint8_t start, uint8_t bitLen, int8_t flag)
{
    assert(1 == start || 0 == start);
    uint64_t mask = (1LL << start);
    //we could have an overflow here,
    //say given bitLen = 2, index 0 move west(return 2), or index 2 move east(return 0)
    for(int i = start; i < bitLen; i += 2)
    {
        if(0 == (val & mask))
        {
            val |= mask;
            if(0 == flag) // addition break here
                break;
        }
        else//val has 1 set at pos i
        {
            val &= (val ^ mask);
            if(1 == flag) //substraction break here
                break;
        }
        mask <<= 2;
    }
    return val;
}

uint64_t naiveMoveXByOne(uint64_t hash, const uint8_t bitLen, bool moveEast)
{
    return naiveMoveByOne(hash & EvenBitsMask, 1, bitLen, moveEast ? 0 : 1) | (hash & OddBitsMask);
}
uint64_t naiveMoveYByOne(uint64_t hash, const uint8_t bitLen, bool moveNorth)
{
    return naiveMoveByOne(hash & OddBitsMask, 0, bitLen, moveNorth ? 0 : 1) | (hash & EvenBitsMask);
}
