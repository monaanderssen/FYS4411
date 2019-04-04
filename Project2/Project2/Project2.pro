TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_MAC_SDK=macosx10.14

INCLUDEPATH += /usr/local/include

LIBS += -L/usr/local/lib

LIBS += -larmadillo -llapack -lblas

SOURCES += \
        main.cpp \
    particle.cpp \
    wavefunction.cpp


QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

HEADERS += \
    particle.h \
    wavefunction.h \
    metropolis.h
