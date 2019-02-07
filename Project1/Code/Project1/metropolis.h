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
    Metropolis(T PDF_){PDF = PDF_;}
    double bruteForceStep(double step);
    vec bruteForceSolve(double step, int iterations);
    T PDF;
    double p1;
    int n = PDF.getNumberOfParticles();
    int dimension = PDF.getDimension();
};


template <class T>
double Metropolis<T>::bruteForceStep(double step){
    mat r_n(dimension,n);
    r_n.randu();
    r_n = r_n*step;
    double Aij = randu();
    for(int i = 0; i < n; i++){
        PDF.getParticle(i).changePosition(r_n.col(i));
    }
    double p2 = PDF.PDF();
    if(Aij <= p2/p1){
        p1 = p2;
        return p2*PDF.localEnergy();

    }
    else{
        for(int i = 0; i < n; i++){
            r_n = - r_n;
            PDF.getParticle(i).changePosition(r_n.col(i));
        }
        return p1*PDF.localEnergy();

    }
}

template <class T>
vec Metropolis<T>::bruteForceSolve(double step, int iterations){
    p1 = PDF.PDF();
    vec E(iterations);
    for(int i = 0; i < iterations; i++){
        E(i) = bruteForceStep(step);
    }
    return E;
}
