#include "unittest.h"

UnitTest::UnitTest(QObject *parent) : QObject(parent)
{

}
//http://stackoverflow.com/questions/686353/c-random-float-number-generation
//range: [low, high)
double UnitTest::randomDouble(double low, double high)
{
    return low + random()/(RAND_MAX/(high - low));
}
void UnitTest::printEnvelope(uint64_t id, Envelope &&env)
{
    cout<< std::fixed << id <<": ("<<env.xmin<<","<<env.ymin << ") , (" << env.xmax << ","<< env.ymax<<")" << endl;
}
void UnitTest::preparePointData(vector<pair<double,double> > &pnts, const Envelope& bbox, const uint64_t size)
{
    set<pair<double,double> > container;
    srand(time(NULL));
    while(container.size() < size)
        container.insert(std::make_pair<double,double>(randomDouble(bbox.xmin, bbox.xmax) ,randomDouble(bbox.ymin, bbox.ymax)));
    for(auto itor = container.begin(); itor != container.end(); ++itor)
        pnts.push_back(std::move(*itor));
}
void UnitTest::TestHaversineDistance()
{
    cout << "Expect 785330.:" << std::setprecision(20) << lonlatDistHaversine(0., 0., 5., 5.) << endl;
    cout << "Expect around 13860000:" << std::setprecision(20) << lonlatDistHaversine(-5.3, 50.2, 70.234, -58.4) << endl;
    cout << "Expect 0:" << std::setprecision(20) << lonlatDistHaversine(0, 0, 0, 0) << endl;
    cout << "Expect 0:" << std::setprecision(20) << lonlatDistHaversine(100, 100, 100, 100) << endl;
    cout << "Expect around 0:" << std::setprecision(20) << lonlatDistHaversine(-180, 0, 180, 0) << endl;
    cout << "Expect around 200020km:" << std::setprecision(20) << lonlatDistHaversine(0, 90, 0, -90) << endl;
    cout << "Expect around 200020km:" << std::setprecision(20) << lonlatDistHaversine(0, 0, 0, 180) << endl;
    cout << "Expect around 200020km:" << std::setprecision(20) << lonlatDistHaversine(180, 0, 0, 0) << endl;
    cout << "Expect around 0:" << std::setprecision(20) << lonlatDistHaversine(360, 0, 0, 0) << endl;
}
void UnitTest::TestGeohash()
{
    uint64_t idd = naiveEncodeGeohash(90,90,64);
    printEnvelope(idd, std::move(naiveDecodeGeohash(idd, 64)));
    idd = naiveEncodeGeohash(45.5,45.5,64);
    printEnvelope(idd, naiveDecodeGeohash(idd, 64));
    idd = naiveEncodeGeohash(-85.8, 43.5,64);
    printEnvelope(idd, naiveDecodeGeohash(idd, 64));

    idd = naiveEncodeGeohash(4842397.8495074, -10601981.5353662, 52, wmRange);
    printEnvelope(idd, naiveDecodeGeohash(idd, 52, wmRange));
    idd = naiveEncodeGeohash(4842397.8495074, -10601981.5353662, 50, wmRange);
    printEnvelope(idd, naiveDecodeGeohash(idd, 50, wmRange));

    idd = naiveEncodeGeohash(0, 0, 2, wmRange);
    printEnvelope(idd, naiveDecodeGeohash(idd, 2, wmRange));

    double radius = 50.0;
    uint8_t bits = guessNumberOfBits(radius, 30.0);
    idd = naiveEncodeGeohash(4842397.8495074, -10601981.5353662, bits, wmRange);
    printEnvelope(idd, naiveDecodeGeohash(idd, bits, wmRange));
    bits = 2*estimate_geohash_steps_by_radius(radius);
    idd = naiveEncodeGeohash(4842397.8495074, -10601981.5353662, bits, wmRange);
    printEnvelope(idd, naiveDecodeGeohash(idd, bits, wmRange));

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
void UnitTest::TestGeohash_CellLowerUpperRange()
{
    //4 bits: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
    //this forms a 16-cell grid that is just sufficient to test all neigbors around the one at the center
    //This 4 by 4 grid's indexes looks like:
    // 5, 7, 13, 15
    // 4, 6, 12, 14
    // 1, 3,  9, 11
    // 0, 2,  8, 10
    //Given center at index 3:
    uint64_t hash3 = naiveEncodeGeohash(-1, -1, 4);
    Envelope env = naiveDecodeGeohash(hash3+1, 4);
    cout << "(" << env.xmin << "," << env.ymin << " ; " << env.xmax << "," << env.ymax << ")" << endl;

    uint64_t min, max;
    naiveCellMaxMin(hash3, 4, &min, &max, 8);
    //hash3: 011, min:0110000, max:1000000 after left shift to occupy 8 bits
    //cout << "Expect 48,64: " << min << "," << max << endl;
    QCOMPARE(48ul, min);
    QCOMPARE(64ul, max);
    //demonstrate pyramids from 2 bits to 4 bits
    //          5, 7, 13, 15
    //1,3       4, 6, 12, 14
    //0,2       1, 3,  9, 11
    //          0, 2,  8, 10
    //
    //index 0 at 2-bits level will be transformed to the (0,1,2,3) quadrant, 4 is the max-exclusive
    naiveCellMaxMin(0, 2, &min, &max, 4);
    //cout << "Expect 0,4: " << min << "," << max << endl;
    QCOMPARE(0ul, min);
    QCOMPARE(4ul, max);
    naiveCellMaxMin(1, 2, &min, &max, 4);
    //cout << "Expect 4,8: " << min << "," << max << endl;
    QCOMPARE(4ul, min);
    QCOMPARE(8ul, max);
    naiveCellMaxMin(2, 2, &min, &max, 4);
    //cout << "Expect 8,12: " << min << "," << max << endl;
    QCOMPARE(8ul, min);
    QCOMPARE(12ul, max);
    naiveCellMaxMin(3, 2, &min, &max, 4);
    //cout << "Expect 12,16: " << min << "," << max << endl;
    QCOMPARE(12ul, min);
    QCOMPARE(16ul, max);
}
void UnitTest::TestGeohash_Neighbors()
{
    //z-order curve, index 0 is at lower left corner
    //2 bits: 0,1,2,3
    //1,3
    //0,2
    QCOMPARE(0ul, naiveEncodeGeohash(-180, -90, 2));
    QCOMPARE(1ul, naiveEncodeGeohash(-180, 0, 2));
    QCOMPARE(2ul, naiveEncodeGeohash(180, -90, 2));
    QCOMPARE(3ul, naiveEncodeGeohash(180, 90, 2));

    //To find a range, we just increment the hash by 1 to find the upper range at of the box


/************************Move East-West***************************************************/
    cout << "Move East-West" << endl;
    QCOMPARE(2ul, naiveMoveXByOne(0, 2));
    QCOMPARE(0ul, naiveMoveXByOne(2, 2, false));
    QCOMPARE(3ul, naiveMoveXByOne(1, 2));
    QCOMPARE(1ul, naiveMoveXByOne(3, 2, false));

    //outliers
    QCOMPARE(2ul, naiveMoveXByOne(2, 2));
    QCOMPARE(0ul, naiveMoveXByOne(0, 2, false));
    QCOMPARE(3ul, naiveMoveXByOne(3, 2));
    QCOMPARE(1ul, naiveMoveXByOne(1, 2, false));

/************************Move North-South***************************************************/
    cout << "Move North-South" << endl;
    QCOMPARE(1ul, naiveMoveYByOne(0, 2));
    QCOMPARE(0ul, naiveMoveYByOne(1, 2, false));
    QCOMPARE(3ul, naiveMoveYByOne(2, 2));
    QCOMPARE(2ul, naiveMoveYByOne(3, 2, false));

    //outliers
    QCOMPARE(1ul, naiveMoveYByOne(1, 2));
    QCOMPARE(0ul, naiveMoveYByOne(0, 2, false));
    QCOMPARE(3ul, naiveMoveYByOne(3, 2));
    QCOMPARE(2ul, naiveMoveYByOne(2, 2, false));

    //4 bits: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
    //this forms a 16-cell grid that is just sufficient to test all neigbors around the one at the center
    //This 4 by 4 grid's indexes looks like:
    // 5, 7, 13, 15
    // 4, 6, 12, 14
    // 1, 3,  9, 11
    // 0, 2,  8, 10
    //Given center at index 12:
    uint64_t hash12 = naiveEncodeGeohash(0, 0, 4);
    QCOMPARE(hash12, 12ul);//1100
    //Find neighbor at east by adding 1 on even bit set (even is defined by couting from left, which is
    //the starting point when the longitude is hashed. So even:[1,0] is [1,1] now.
    cout << "Index 12's East:" << endl;
    QCOMPARE(naiveMoveXByOne(hash12, 4), 14ul); //1110
    assert(14 == naiveNeighbor(hash12, 4, GH_DIRECTION::EAST));
    cout << "Index 12's West:" << endl;
    QCOMPARE(naiveMoveXByOne(hash12, 4, false), 6ul); //0110, move west
    //assert(6 == naiveNeighbor(hash12, 4, GH_DIRECTION::WEST));
    cout << "Index 12's North:" << endl;
    QCOMPARE(naiveMoveYByOne(hash12, 4), 13ul); //1101, move north
    //assert(13 == naiveNeighbor(hash12, 4, GH_DIRECTION::NORTH));
    cout << "Index 12's West:" << endl;
    QCOMPARE(naiveMoveYByOne(hash12, 4, false), 9ul); //1001, move south
    //assert(9 == naiveNeighbor(hash12, 4, GH_DIRECTION::SOUTH));
    cout << "Index 12's NorthEast:" << endl;
    QCOMPARE(naiveNeighbor(hash12, 4, GH_DIRECTION::NORTHEAST), 15ul); //1111, move northeast
    cout << "Index 12's NorthWest:" << endl;
    QCOMPARE(naiveNeighbor(hash12, 4, GH_DIRECTION::NORTHWEST), 7ul); //0111, move northwest
    cout << "Index 12's SouthEast:" << endl;
    QCOMPARE(naiveNeighbor(hash12, 4, GH_DIRECTION::SOUTHEAST), 11ul); //1011, move southeast
    cout << "Index 12's SouthWest:" << endl;
    QCOMPARE(naiveNeighbor(hash12, 4, GH_DIRECTION::SOUTHWEST), 3ul); //0011, move southeast

    //outliers:
    uint64_t hash15 = naiveEncodeGeohash(180,90,4);
    QCOMPARE(hash15, 15ul);//1111
    QCOMPARE(naiveNeighbor(hash15, 4, GH_DIRECTION::NORTHEAST), 15ul);
    QCOMPARE(naiveNeighbor(hash15, 4, GH_DIRECTION::NORTHWEST), 13ul);
    QCOMPARE(naiveNeighbor(hash15, 4, GH_DIRECTION::SOUTHEAST), 14ul);
    QCOMPARE(naiveNeighbor(hash15, 4, GH_DIRECTION::SOUTHWEST), 12ul);
    QCOMPARE(naiveNeighbor(hash15, 4, GH_DIRECTION::EAST), 15ul);
    QCOMPARE(naiveNeighbor(hash15, 4, GH_DIRECTION::WEST), 13ul);
}
void UnitTest::TestGeohash_QueryByRadius()
{
    const uint8_t bitDepth = 52;
    const uint32_t pntCount = 1e6;
    vector<pair<double,double> > pnts;
    /*
     * Filling random point data to emulate points on the ground
    */
    preparePointData(pnts, llRange, pntCount);
    cout << pntCount << " points are ready" << endl;
    //my centos 6 equips with g++4.4.7 that doesn't support c++11
    //for(auto p:pnts)centerHash
    //for(auto p = pnts.begin(); p != pnts.end(); ++p)
    //    cout << (*p).first <<"," << p->second <<endl;

    /*
     * First, let's geohash all points and save them in a Set so we can conduct binary search later
     * with guaranteed O(logN)
     * Default bit depth is 52
    */
    set<uint64_t> container;
    for(auto p = pnts.begin(); p != pnts.end(); ++p)
        container.insert(naiveEncodeGeohash(p->first, p->second, bitDepth));
    cout << "Point data have been geohashed" << endl;

    /*
     * Search POI surrounding my location:
     *
     * Set up the lat/lon as the center for search and of course, the radius
     *
    */
    double radius = 1e5;
    double centerLon = 0.0;
    double centerLat = 0.0;
    /*
     * First let's get the number of bits that are necessary to represent this central lon/lat based upon radius.
     * In other word, at this particular bit depth level, this raduis is recognizable (has sufficient resolution)
     * , which means the cell at this bit depth level is greater than the radius.
    */
    uint8_t estBitDepth = guessNumberOfBits(radius, centerLat);
    uint64_t min = 0ULL, max = 0ULL;
    /*
     * Encode central lon/lat
    */
    uint64_t centerHash = naiveEncodeGeohash(centerLon, centerLat, estBitDepth);
    /*
     * According to the hash value at estimated bit depth level, find out what are the corresponding min, max values at level 52.
     * In other words, with this cell size (at estimated bit-depth level), how many cells at bit-depth level 52 should be covered.
    */
    naiveCellMaxMin(centerHash, estBitDepth, &min, &max, bitDepth);
    /*
     * Define search ranges
    */
    vector<pair<uint64_t, uint64_t> > ranges;
    ranges.push_back(std::make_pair(min, max));
    GH_DIRECTION dirs[]{GH_DIRECTION::NORTH, GH_DIRECTION::EAST, GH_DIRECTION::NORTHEAST, GH_DIRECTION::NORTHWEST, GH_DIRECTION::SOUTH, GH_DIRECTION::SOUTHEAST, GH_DIRECTION::SOUTHWEST, GH_DIRECTION::WEST};
    uint64_t neighborHash = 0ULL;
    /*
     * Put all neighbors in, be awared: we didn't try to optimize this by merging overlapping ranges or removing duplicates or invalid,
     * which should be done in a real world app.
     */
    for(int i = 0; i < 8; ++i)
    {
        neighborHash = naiveNeighbor(centerHash, estBitDepth, dirs[i]);
        naiveCellMaxMin(neighborHash, estBitDepth, &min, &max, bitDepth);
        ranges.push_back(std::make_pair(min, max));
    }
    set<uint64_t>::iterator low;
    uint32_t count = 0;
    Envelope dll(0, 0, 0, 0);
    /*
     * Use a set to hold the hits through geohash-based query in order to avoid any duplicates
    */
    set<uint64_t> result;
    for(int i = 0; i < 9; ++i)
    {
        /*
         * This lower bound value is the starting point, if its value is no greater than ranges[i].second,
         * we consider all values in between [low, ranges[i].second) as candidates.
         * Then, we will have to scan all of them to get rid of false candidates.
        */
        low = container.lower_bound(ranges[i].first);
        for(auto itor = low; *itor < ranges[i].second; ++itor)
        {
            dll = naiveDecodeGeohash(*itor, bitDepth);
            /*
             * Haversine distance
             * We decode hash into an envelope and then use the center of this envelop to calcuate the distance to the central point
            */
            if(lonlatDistHaversine(centerLon, centerLat, (dll.xmax + dll.xmin)/2.0, (dll.ymax + dll.ymin)/2.0) <= radius)
                result.insert(*itor);
        }
    }
    cout << "Geohash found " << result.size() << " hits!"<<endl;
    /*
     * Verify the correctness of geohash based query by comparing the number of hits with a brute-force linear scanning against entire
     * dataset.
    */
    for(auto itor = container.begin(); itor != container.end(); ++itor)
    {
        dll = naiveDecodeGeohash(*itor, bitDepth);
        if(lonlatDistHaversine(centerLon, centerLat, (dll.xmax + dll.xmin)/2.0, (dll.ymax + dll.ymin)/2.0) <= radius)
            ++count;
    }
    cout << "Linear scan found " << count << " hits!"<< endl;
    QCOMPARE(count, (uint32_t)result.size());
}

QTEST_MAIN(UnitTest)
