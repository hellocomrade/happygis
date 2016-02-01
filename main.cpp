#include <iostream>
#include <iomanip>
#include "binutils.h"
#include "morton.h"
#include "envelope.h"
#include "latlon2Dist.h"
using namespace std;
void printEnv(uint64_t id, const Envelope &env)
{
    cout<< std::fixed << id <<": ("<<env.xmin<<","<<env.ymin<<") , ("<<env.xmax<<","<<env.ymax<<")"<<endl;
}
void testGeohash()
{
    BinUtils bin(1.0);
    uint64_t id=bin.getId(0,1);
    printEnv(id,bin.queryEnvelope(id));
    id=bin.getId(1,1);
    printEnv(id,bin.queryEnvelope(id));
    id=bin.getId(0,0);
    printEnv(id,bin.queryEnvelope(id));
    id=bin.getId(1,0);
    printEnv(id,bin.queryEnvelope(id));
    id=bin.getId(-1,0);
    printEnv(id,bin.queryEnvelope(id));

    uint64_t idd = naiveEncodeGeohash(90,90,64);
    printEnv(idd, naiveDecodeGeohash(idd, 64));
    idd = naiveEncodeGeohash(45.5,45.5,64);
    printEnv(idd, naiveDecodeGeohash(idd, 64));
    idd = naiveEncodeGeohash(-85.8, 43.5,64);
    printEnv(idd, naiveDecodeGeohash(idd, 64));

    idd = naiveEncodeGeohash(4842397.8495074, -10601981.5353662, 52, wmRange);
    printEnv(idd, naiveDecodeGeohash(idd, 52, wmRange));
    idd = naiveEncodeGeohash(4842397.8495074, -10601981.5353662, 50, wmRange);
    printEnv(idd, naiveDecodeGeohash(idd, 50, wmRange));

    idd = naiveEncodeGeohash(0, 0, 2, wmRange);
    printEnv(idd, naiveDecodeGeohash(idd, 2, wmRange));

    double radius = 50.0;
    uint8_t bits = guessNumberOfBits(radius);
    idd = naiveEncodeGeohash(4842397.8495074, -10601981.5353662, bits, wmRange);
    printEnv(idd, naiveDecodeGeohash(idd, bits, wmRange));
    bits = 2*estimate_geohash_steps_by_radius(radius);
    idd = naiveEncodeGeohash(4842397.8495074, -10601981.5353662, bits, wmRange);
    printEnv(idd, naiveDecodeGeohash(idd, bits, wmRange));

    //epsg:3857
    for(int i=0;i < 65; i+=2)
    {
        //-20037508.34,-20037508.34,20037508.34,20037508.34
        idd = naiveEncodeGeohash(4842397.8495074, -10601981.5353662, i, wmRange);
        Envelope e = naiveDecodeGeohash(idd, i, wmRange);
        cout << i << " : (" << (e.xmax - e.xmin) << ", " << (e.ymax - e.ymin) << ")" << endl;
    }
    //epsg:4326
    for(int i=0;i < 65; i+=2)
    {
        uint64_t idd = naiveEncodeGeohash(0, 0, i);
        Envelope e = naiveDecodeGeohash(idd, i, llRange);
        cout << i << " : " << std::setprecision(20) << lonlatDistHaversine(e.xmax, e.ymax, e.xmin, e.ymin) << endl;
    }
}
void testGeohashNeighbors()
{
    //z-order curve, index 0 is at lower left corner
    //2 bits: 0,1,2,3
    //1,3
    //0,2
    cout << "Index 0: " << naiveEncodeGeohash(-180, -90, 2)<<endl;
    cout << "Index 1: " << naiveEncodeGeohash(-180, 0, 2)<<endl;
    cout << "Index 2: " << naiveEncodeGeohash(180, -90, 2)<<endl;
    cout << "Index 3: " << naiveEncodeGeohash(180, 90, 2)<<endl;

/************************Move East-West***************************************************/
    cout << "Move East-West" << endl;
    cout << "Expect 2: " << naiveMoveXByOne(0, 2) << endl;
    cout << "Expect 0: " << naiveMoveXByOne(2, 2, false) << endl;
    cout << "Expect 3: " << naiveMoveXByOne(1, 2) << endl;
    cout << "Expect 1: " << naiveMoveXByOne(3, 2, false) << endl;

    //outliers
    cout << "Expect 2: " << naiveMoveXByOne(2, 2) << endl;
    cout << "Expect 0: " << naiveMoveXByOne(0, 2, false) << endl;
    cout << "Expect 3: " << naiveMoveXByOne(3, 2) << endl;
    cout << "Expect 1: " << naiveMoveXByOne(1, 2, false) << endl;

/************************Move North-South***************************************************/
    cout << "Move North-South" << endl;
    cout << "Expect 1: " << naiveMoveYByOne(0, 2) << endl;
    cout << "Expect 0: " << naiveMoveYByOne(1, 2, false) << endl;
    cout << "Expect 3: " << naiveMoveYByOne(2, 2) << endl;
    cout << "Expect 2: " << naiveMoveYByOne(3, 2, false) << endl;

    //outliers
    cout << "Expect 1: " << naiveMoveYByOne(1, 2) << endl;
    cout << "Expect 0: " << naiveMoveYByOne(0, 2, false) << endl;
    cout << "Expect 3: " << naiveMoveYByOne(3, 2) << endl;
    cout << "Expect 2: " << naiveMoveYByOne(2, 2, false) << endl;

    //4 bits: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
    //this forms a 16-cell grid that is just sufficient to test all neigbors around the one at the center
    //This 4 by 4 grid's indexes looks like:
    // 5, 7, 13, 15
    // 4, 6, 12, 14
    // 1, 3,  9, 11
    // 0, 2,  8, 10
    //Given center at index 12:
    uint64_t hash12 = naiveEncodeGeohash(0, 0, 4);
    cout << "Index 12: " << hash12 <<endl;//1100
    //Find neighbor at east by adding 1 on even bit set (even is defined by couting from left, which is
    //the starting point when the longitude is hashed. So even:[1,0] is [1,1] now.
    cout << "Index 12's East:" << endl;
    cout << "Expect 14: " << naiveMoveXByOne(hash12, 4) << endl; //1110
    cout << "Index 12's West:" << endl;
    cout << "Expect 6: " << naiveMoveXByOne(hash12, 4, false) << endl; //0110, move west
    cout << "Index 12's North:" << endl;
    cout << "Expect 13: " << naiveMoveYByOne(hash12, 4) << endl; //1101, move north
    cout << "Index 12's West:" << endl;
    cout << "Expect 9: " << naiveMoveYByOne(hash12, 4, false) << endl; //1001, move south


}

void testlonlatDistance()
{
    cout << "Expect 785329.:" << std::setprecision(20) << lonlatDistHaversine(0., 0., 5., 5.) << endl;
    cout << "Expect around 13860000:" << std::setprecision(20) << lonlatDistHaversine(-5.3, 50.2, 70.234, -58.4) << endl;
    cout << "Expect 0:" << std::setprecision(20) << lonlatDistHaversine(0, 0, 0, 0) << endl;
    cout << "Expect 0:" << std::setprecision(20) << lonlatDistHaversine(100, 100, 100, 100) << endl;
    cout << "Expect around 0:" << std::setprecision(20) << lonlatDistHaversine(-180, 0, 180, 0) << endl;
    cout << "Expect around 200020km:" << std::setprecision(20) << lonlatDistHaversine(0, 90, 0, -90) << endl;
    cout << "Expect around 200020km:" << std::setprecision(20) << lonlatDistHaversine(0, 0, 0, 180) << endl;
    cout << "Expect around 200020km:" << std::setprecision(20) << lonlatDistHaversine(180, 0, 0, 0) << endl;
    cout << "Expect around 0:" << std::setprecision(20) << lonlatDistHaversine(360, 0, 0, 0) << endl;
}

int main()
{
    testlonlatDistance();
    testGeohashNeighbors();
    testGeohash();
    return 0;
}



