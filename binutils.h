#ifndef BINUTILS_H
#define BINUTILS_H
#include <cmath>
#include <limits>
#include "envelope.h"
class BinUtils
{
private:
    long long numCols;
    double extentMin;
    double extentMax;
    const double binSize;
public:
    BinUtils(double binSize);
    long long getId(double x,double y);
    Envelope queryEnvelope(long long id);

};

#endif // BINUTILS_H
