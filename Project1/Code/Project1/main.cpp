#include <iostream>
#include "sphericalharmonicoscillator.h"
#include "metropolis.h"

using namespace std;
using namespace arma;

int main()
{
    SphericalHarmonicOscillator test;
    test.setDimension(3);
    test.setNumberOfParticles(10);
    test.setAlpha(0.5);
    test.setPsi();
    Metropolis<SphericalHarmonicOscillator> ttt(test);
    vec E = ttt.bruteForceSolve(0.5, 1000);
    cout<< mean(E) << endl;
    return 0;
}
