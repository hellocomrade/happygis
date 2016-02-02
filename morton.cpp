#include "morton.h"
//http://gis.stackexchange.com/questions/18330/would-it-be-possible-to-use-geohash-for-proximity-searches

/*
From: https://github.com/davetroy/geohash-js

SUMMARY
This is a basic demonstration of how the GeoHash algorithm can be used to generate bounding box searches without the use
of specialized spatial indexing approaches.

BACKGROUND
The Geohash algorithm was first described by Gustavo Niemeyer in February 2008.  By interleaving latitude and longitude
information in a bitwise fashion, a composite value is generated that provides a high resolution geographic point, and is
well suited for storage or transmission as a character string.

Geohash also has the property that as the number of digits decreases (from the right), accuracy degrades.  This property
can be used to do bounding box searches, as points near to one another will share similar Geohash prefixes.

However, because a given point may appear at the edge of a given Geohash bounding box, it is necessary to generate a list
of Geohash values in order to perform a true proximity search around a point.  Because the Geohash algorithm uses a base-32
numbering system, it is possible to derive the Geohash values surrounding any other given Geohash value using a simple
lookup table.
*/

/*
 *
 * When conducting the bbox query, we'd like to find a 'precision' or the number of bits that could generate a hash to cover
 * a bound that is large enough to cover a diameter that is 2 times of the given radius.
 *
 * This can only be done on a projected coordinate system under Euclidean space and The most popular case is
 * Google Web Mercator.
 *
 * We have to keep in mind: every time we split/bisect the plane, the diameter of the bbox will be cut to half. So, in the worst
 * case, we can NOT split the plane if the radius is larger than wmRange.xmax
 *
 * Below is the impl from: https://github.com/yinqiwen/ardb/
 * which doesn't make sense to me, for (radius > 0.5 * wmRange.xmax and radius < wmRange.xmax),
 * the plane will still be split, then radius * 2 will be greater than the diameter of the new bound!
 * static uint8_t estimate_geohash_steps_by_radius(double range_meters)
    {
        uint8_t step = 1;
        double v = range_meters;
        while (v < wmRange.xmax)
        {
            v *= 2;
            step++;
        }
        step--;
        return step;
    }

Here is an accepted version at github for redis:

https://github.com/antirez/redis/blob/846da5b22e7877614694b4bda94fe0fdca34ada9/deps/geohash-int/geohash_helper.c

and here is the redis impl:

https://github.com/antirez/redis/blob/unstable/src/geo.c#L281
https://github.com/antirez/redis/blob/unstable/deps/geohash-int/geohash_helper.c

 */
uint8_t guessNumberOfBits(double radius)
{
    uint8_t num = 0;
    if(radius > 0)
    {
        while( (radius *= 2) < wmRange.xmax)
            num += 2; //num is always an odd number since both lon and lat need identical numbers of bit
    }
    return num;
}
uint8_t estimate_geohash_steps_by_radius(double range_meters)
{
    uint8_t step = 1;
    double v = range_meters;
    while (v < wmRange.xmax)
    {
        v *= 2;
        step++;
    }
    step--;
    return step;
}

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
            if((code & (1ULL<<(pos-i))) == 0) //lng < midpoint
                env.xmax = (env.xmax + env.xmin) / 2;
            else
                env.xmin = (env.xmax + env.xmin) / 2;
        }
        else
        {
            if((code & (1ULL<<(pos-i))) == 0) //lat < midpoint
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
    uint64_t mask = (1ULL << start);
    //we could have an overflow here due to the size of the bit sequence,
    //say given bitLen = 2, index 0 move west(return 2), or index 2 move east(return 0)
    uint64_t oldVal = val;
    for(int i = start; i < bitLen; i += 2)
    {
        if(0 == (val & mask))
        {
            val |= mask;
            if(0 == flag) // addition break here
                return val;
        }
        else//val has 1 set at pos i
        {
            val &= (val ^ mask);
            if(1 == flag) //substraction break here
                return val;
        }
        mask <<= 2;
    }
    return oldVal;
}

uint64_t naiveMoveXByOne(uint64_t hash, const uint8_t bitLen, bool moveEast)
{
    return naiveMoveByOne(hash & EvenBitsMask, 1, bitLen, moveEast ? 0 : 1) | (hash & OddBitsMask);
}
uint64_t naiveMoveYByOne(uint64_t hash, const uint8_t bitLen, bool moveNorth)
{
    return naiveMoveByOne(hash & OddBitsMask, 0, bitLen, moveNorth ? 0 : 1) | (hash & EvenBitsMask);
}
uint64_t naiveNeighbor(uint64_t hash, uint8_t bitLen, GH_DIRECTION dir)
{
    uint64_t val;
    switch(dir)
    {
    case GH_DIRECTION::NORTH:
        val = naiveMoveYByOne(hash, bitLen);
        break;
    case GH_DIRECTION::SOUTH:
        val = naiveMoveYByOne(hash, bitLen, false);
        break;
    case GH_DIRECTION::EAST:
        val = naiveMoveXByOne(hash, bitLen);
        break;
    case GH_DIRECTION::WEST:
        val = naiveMoveXByOne(hash, bitLen, false);
        break;
    case GH_DIRECTION::NORTHEAST:
        val = (naiveMoveByOne(hash & OddBitsMask, 0, bitLen, 0) | naiveMoveByOne(hash & EvenBitsMask, 1, bitLen, 0));
        break;
    case GH_DIRECTION::NORTHWEST:
        val = (naiveMoveByOne(hash & OddBitsMask, 0, bitLen, 0) | naiveMoveByOne(hash & EvenBitsMask, 1, bitLen, 1));
        break;
    case GH_DIRECTION::SOUTHEAST:
        val = (naiveMoveByOne(hash & OddBitsMask, 0, bitLen, 1) | naiveMoveByOne(hash & EvenBitsMask, 1, bitLen, 0));
        break;
    case GH_DIRECTION::SOUTHWEST:
        val = (naiveMoveByOne(hash & OddBitsMask, 0, bitLen, 1) | naiveMoveByOne(hash & EvenBitsMask, 1, bitLen, 1));
        break;
    default:
        val = hash;
        break;
    }
    return val;
}
/* Compute the sorted set scores min (inclusive), max (exclusive) we should
 * query in order to retrieve all the elements inside the specified area
 * 'hash'. The two scores are returned by reference in *min and *max. */
void naiveCellMaxMin(uint64_t hash, uint8_t bitLen, uint64_t *min, uint64_t *max, uint8_t targetBitLen)
{
    /* We want to compute the sorted set scores that will include all the
     * elements inside the specified Geohash 'hash', which has as many
     * bits as specified by targetBitLen.
     *
     * So if bitLen is, for example, 6, and the hash value in binary
     * is 101010, since our score is 52 bits we want every element which
     * is in binary: 101010?????????????????????????????????????????????
     * Where ? can be 0 or 1.
     *
     * To get the min score we just use the initial hash value left
     * shifted enough to get the 52 bit value. Later we increment the
     * 6 bit prefis (see the hash.bits++ statement), and get the new
     * prefix: 101011, which we align again to 52 bits to get the maximum
     * value (which is excluded from the search). So we get everything
     * between the two following scores (represented in binary):
     *
     * 1010100000000000000000000000000000000000000000000000 (included)
     * and
     * 1010110000000000000000000000000000000000000000000000 (excluded).
     */
    *min = hash;
    *min <<= targetBitLen - bitLen;
    *max = ++hash;
    *max <<= targetBitLen - bitLen;
}
static uint64_t move_x(uint64_t code, int8_t d, const uint8_t bitLen)
{
    if (d == 0)
        return code;
    uint64_t x = code & 0xaaaaaaaaaaaaaaaaULL; //isolate even bits
    uint64_t y = code & 0x5555555555555555ULL; //isolate odd bits
    /*
     * Create a mask so we could transmit addition/subtraction carry/borrow onto even bits.
     * The idea is to make sure all odd bits are 1s, so if there is a carry by adding one,
     * it's guaranteed there will be a carry over from odd bit to even bit.
     *
     * Given hash: 1001, mask: 0101
     *      Isolated even bits rep: 1000, mask + 1 = 0110
     *       1000
     *     + 0110 = 1110
     * For this case, there is no carry :)
     *
     * Given hash: 1011, mask: 0101
     *       Isolated even bits rep: 1010, mask + 1 = 0110
     *        1010
     *      + 0110 = 10000, see the carry from the second bit is transmitted all the way to the left?
     *
     *      It's like: 011 + 001 = 100, if we only look at the even bits
    */
    uint64_t msk = 0x5555555555555555ULL >> (64 - bitLen);
    if (d > 0)
        x = x + (msk + 1);
    else
    {
        x = x | msk;
        x = x - (msk + 1);
    }
    x &= (0xaaaaaaaaaaaaaaaaULL >> (64 - bitLen));
    return (x | y);
}
