#pragma once
#include "wavefunction.h"
#include "particle.h"
#include <cstdlib>
#include <random>
#include <armadillo>
#include <iostream>
#include <fstream>

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
	void bruteForceTofile(double step, int iterations); //runs over diferent values of alpha
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

template <class T>
void Metropolis<T>::bruteForceTofile(double step, int iterations) {
	string fileName;
	fileName += "MetropolisSphericalHarmonicsStep" + to_string(step) + "Iterations" + to_string(iterations) + ".txt";
	ofstream myFile;
	myFile.open(fileName);
	double startAlpha = 0.1;
	double alphaStep = 0.01;
	for (int i = 0; i < 100 ; i++) {
		double tempAlpha = i * alphaStep + startAlpha;
		PDF.setAlpha(tempAlpha);
		vec out = bruteForceSolve(step, iterations);
		myFile << tempAlpha << " " << out(0) << " " << out(1) << endl;// " " << E(1) / sqrt(iterations) << endl;
	}
	myFile.close();
}




