#include "metropolis.h"
#include "wavefunction.h"
#include "particle.h"

template <class T>
double Metropolis<T>::bruteForceStep(double step){
    mat r_n(dimension,n);
    r_n.randu();
    r_n = r_n*step;
    double Aij = randu();
    for(int i = 0; i < n; i++){
        PDF.getParticle(i).changePosition(r_n.colptr(i));
    }
    double p2 = PDF.PDF();
    if(Aij <= p2/p1){
        p1 = p2;
        return p2*PDF.localEnergy();

    }
    else{
        for(int i = 0; i < n; i++){
            r_n = - r_n;
            PDF.getParticle(i).changePosition(r_n.colptr(i));
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
