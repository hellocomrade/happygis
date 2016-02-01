#ifndef MORTON_H
#define MORTON_H
#include <cassert>
#include <cstdint>
#include <iostream>
using namespace std;
#include "envelope.h"
const Envelope llRange(-180,-90,180,90);
const Envelope wmRange(-20037508.34,-20037508.34,20037508.34,20037508.34); //epsg:3857
const uint64_t EvenBitsMask = 0xAAAAAAAAAAAAAAAALL;
const uint64_t OddBitsMask = 0x5555555555555555LL;
uint64_t naiveEncodeGeohash(double lat, double lng, uint8_t bitLen, Envelope bbox = llRange);
Envelope naiveDecodeGeohash(const int64_t code, const uint8_t bitLen, Envelope bbox = llRange);
uint8_t guessNumberOfBits(double radius);
uint8_t estimate_geohash_steps_by_radius(double range_meters);
uint64_t naiveMoveXByOne(uint64_t hash, const uint8_t bitLen, bool moveEast = true);
uint64_t naiveMoveYByOne(uint64_t hash, const uint8_t bitLen, bool moveNorth = true);
#endif // MORTON_H#ifndef MORTON_H
