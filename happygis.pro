QT += core
QT += testlib
QT -= gui

TARGET = happygis
CONFIG += console
CONFIG -= app_bundle
CONFIG += testcase
#CONFIG -= c++0x
#CONFIG += c++11
#QMAKE_CXX = g++-4.8
#QMAKE_LINK = g++-4.8
#QMAKE_CXXFLAGS -= -std=c++0x
QMAKE_CXXFLAGS += -std=c++11

TEMPLATE = app

SOURCES += \#main.cpp \
    binutils.cpp \
    morton.cpp \
    latlon2Dist.cpp \
    unittest.cpp

HEADERS += \
    binutils.h \
    envelope.h \
    morton.h \
    latlon2Dist.h \
    unittest.h

