#ifndef MORTON_H
#define MORTON_H
#include <cassert>
#include <cstdint>
#include <iostream>
using namespace std;
#include "envelope.h"
enum class GH_DIRECTION
{
    NORTH,
    SOUTH,
    EAST,
    WEST,
    NORTHEAST,
    SOUTHWEST,
    NORTHWEST,
    SOUTHEAST
};
const Envelope llRange(-180,-90,180,90);
const Envelope wmRange(-20037508.34,-20037508.34,20037508.34,20037508.34); //epsg:3857
const uint64_t EvenBitsMask = 0xAAAAAAAAAAAAAAAAULL;
const uint64_t OddBitsMask = 0x5555555555555555ULL;
uint64_t naiveEncodeGeohash(double lat, double lng, uint8_t bitLen, Envelope bbox = llRange);
Envelope naiveDecodeGeohash(const int64_t code, const uint8_t bitLen, Envelope bbox = llRange);
uint64_t geohash_move_x(int64_t code, int8_t d, const uint8_t bitLen);
uint8_t guessNumberOfBits(double radius);
uint8_t estimate_geohash_steps_by_radius(double range_meters);
uint64_t naiveMoveXByOne(uint64_t hash, const uint8_t bitLen, bool moveEast = true);
uint64_t naiveMoveYByOne(uint64_t hash, const uint8_t bitLen, bool moveNorth = true);
uint64_t naiveNeighbor(uint64_t hash, uint8_t bitLen, GH_DIRECTION dir);
void naiveCellMaxMin(uint64_t hash, uint8_t bitLen, uint64_t *min, uint64_t *max, uint8_t targetBitLen);
#endif // MORTON_H
