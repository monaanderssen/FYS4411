#pragma once
#include "wavefunction.h"
#include "particle.h"
#include <cstdlib>
#include <random>
#include <armadillo>

using namespace arma;
using namespace std;

template <class T>
class Metropolis
{
public:
    Metropolis(T PDF_){
        PDF = PDF_;
        n = PDF.getNumberOfParticles();
        dimension = PDF.getDimension();
    }
    double bruteForceStep(double step);
    vec bruteForceSolve(double step, int iterations);
    T PDF;
    double p1;
    int n;
    int dimension;
};


template <class T>
double Metropolis<T>::bruteForceStep(double step){
    mat r_n(dimension,n);
    r_n.randn();
    r_n *= step;
    double Aij = randu();
    for(int i = 0; i < n; i++){
        PDF.getParticle(i)->changePosition(r_n.col(i));

    }
    double p2 = PDF.PDF();
    if(Aij <= p2/p1){
        p1 = p2;
        return PDF.localEnergy();


    }
    else{
        for(int i = 0; i < n; i++){
            PDF.getParticle(i)->changePosition(-r_n.col(i));
        }
        return PDF.localEnergy();


    }
}

template <class T>
vec Metropolis<T>::bruteForceSolve(double step, int iterations){
    p1 = PDF.PDF();
    vec ret(2);
    double E2 = 0;
    double E = 0;
    double bruteStep;
    for(int i = 0; i < iterations; i++){
        bruteStep = bruteForceStep(step);
        E += bruteStep;
        E2 += bruteStep*bruteStep;
    }
    E /= iterations;
    E2 /= iterations;
    ret(0) = E;
    ret(1) = E2 - E*E;
    return ret;
}
