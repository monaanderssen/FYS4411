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
    test.setNumberOfParticles(5);
    test.setAlpha(0.45);
    test.setomegaHo(1);
    test.setPsi();
    Metropolis<SphericalHarmonicOscillator> ttt(test);
    int iters = 10000;
    vec En = ttt.bruteForceSolve(0.1, iters);
    cout<< En(0) << " " << En(1) << endl;
    return 0;
}
