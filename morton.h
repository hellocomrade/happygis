#ifndef MORTON_H
#define MORTON_H
#include <cassert>
#include <cstdint>
#include "envelope.h"
const Envelope llRange(-180,-90,180,90);
int64_t naiveEncodeGeohash(double lat, double lng, uint8_t bitLen);
Envelope naiveDecodeGeohash(const int64_t code, const uint8_t bitLen);
#endif // MORTON_H
