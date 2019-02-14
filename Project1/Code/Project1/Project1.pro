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
    sphericalharmonicoscillator.cpp \
    particle.cpp \
    wavefunction.cpp \
    harmonicoscillator.cpp

HEADERS += \
    sphericalharmonicoscillator.h \
    particle.h \
    wavefunction.h \
    harmonicoscillator.h \
    metropolis.h


QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3
