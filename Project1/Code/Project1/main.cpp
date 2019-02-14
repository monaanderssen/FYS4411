#include <iostream>
#include "sphericalharmonicoscillator.h"
#include "metropolis.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
    SphericalHarmonicOscillator test;
    test.setDimension(3);
    test.setNumberOfParticles(10);
    test.setAlpha(0.1);
    test.setomegaHo(1);
    test.setPsi();
    Metropolis<SphericalHarmonicOscillator> ttt(test);
    vec E = ttt.bruteForceSolve(1, 1000);
    cout<< mean(E) << endl;
    return 0;
}
