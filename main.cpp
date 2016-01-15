#include <iostream>
#include "binutils.h"
#include "morton.h"
#include "envelope.h"
using namespace std;
void printEnv(long long id, const Envelope &env)
{
    cout<<id<<": ("<<env.xmin<<","<<env.ymin<<") , ("<<env.xmax<<","<<env.ymax<<")"<<endl;
}
int main()
{
    BinUtils bin(1.0);
    long long id=bin.getId(0,1);
    printEnv(id,bin.queryEnvelope(id));
    id=bin.getId(1,1);
    printEnv(id,bin.queryEnvelope(id));
    id=bin.getId(0,0);
    printEnv(id,bin.queryEnvelope(id));
    id=bin.getId(1,0);
    printEnv(id,bin.queryEnvelope(id));
    id=bin.getId(-1,0);
    printEnv(id,bin.queryEnvelope(id));
    int64_t idd = naiveEncodeGeohash(90,90,64);
    printEnv(idd, naiveDecodeGeohash(idd, 64));
    idd = naiveEncodeGeohash(45.5,45.5,64);
    printEnv(idd, naiveDecodeGeohash(idd, 64));
    idd = naiveEncodeGeohash(-75.8,91.5,64);
    printEnv(idd, naiveDecodeGeohash(idd, 64));
    return 0;
}

