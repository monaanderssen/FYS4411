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
    test.setNumberOfParticles(100);
    test.setAlpha(0.1);
    test.setomegaHo(1);
    test.setPsi();
    Metropolis<SphericalHarmonicOscillator> ttt(test);
    int iters = 100000;
    vec En = ttt.bruteForceSolve(0.1, iters);
    cout<< En(0) << " " << En(1) << endl;
    return 0;
}
