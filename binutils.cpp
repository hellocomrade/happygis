#include "binutils.h"

BinUtils::BinUtils(double binSize):binSize(binSize)
{
    long long maxBinsPerAxis=2;//(long long)sqrt(std::numeric_limits<long long>::max());
    double size=(binSize<1)?maxBinsPerAxis*binSize:maxBinsPerAxis;
    extentMax=size/2;
    extentMin=extentMax-size;
    numCols=(long long)ceil(size/binSize);
}
long long BinUtils::getId(double x, double y)
{
    double down=(this->extentMax-y)/this->binSize;
    double over=(x-extentMin)/this->binSize;
    return ((long long)down*numCols)+(long long)over;
}
Envelope BinUtils::queryEnvelope(long long id)
{
    long long down=id/this->numCols;
    long long over=id%this->numCols;
    double xmin=extentMin+(over*this->binSize);
    double xmax=xmin+this->binSize;
    double ymax=extentMax-down*this->binSize;
    double ymin=ymax-this->binSize;
    return Envelope(xmin,ymin,xmax,ymax);
}
