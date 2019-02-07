TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    sphericalharmonicoscillator.cpp \
    particle.cpp \
    wavefunction.cpp \
    harmonicoscillator.cpp \
    metropolis.cpp

HEADERS += \
    sphericalharmonicoscillator.h \
    particle.h \
    wavefunction.h \
    harmonicoscillator.h \
    metropolis.h
