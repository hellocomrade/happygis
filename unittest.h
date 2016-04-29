#ifndef UNITTEST_H
#define UNITTEST_H
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <set>
#include <QtTest/QtTest>
#include "latlon2Dist.h"
#include "envelope.h"
#include "morton.h"

using namespace std;

class UnitTest : public QObject
{
    Q_OBJECT
public:
    explicit UnitTest(QObject *parent = 0);
private:
    inline double randomDouble(double low, double high);
    inline void printEnvelope(uint64_t id, Envelope &&env);
    void preparePointData(vector<pair<double,double> > &pnts, const Envelope& bbox, const uint64_t size);

signals:

private slots:
    void TestHaversineDistance();
    void TestGeohash();
    void TestGeohash_CellLowerUpperRange();
    void TestGeohash_Neighbors();
    void TestGeohash_QueryByRadius();
};

#endif // UNITTEST_H
