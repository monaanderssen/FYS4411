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
	vec importantSamplingStep(double timeStep);
	vec importantSamplingSolve(double timeStep, int iterations); //finds only localEnergy
	void importantSamplingToFile(double timeStep, int inerations);
	vec dEnergy(double timeStep); //returns a vector with energy and derivative of energy using importantSampling
	vec minimice(double timeStep,double startAlpha,double tol, double stepLength); //finding minimum using stepest desent
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
	fileName += "BfMcDim" + to_string(dimension)+ "Npart"+to_string(n) + "Iter" + to_string(iterations) + ".txt";
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


template <class T>
vec Metropolis<T>::importantSamplingStep(double timeStep) {
	mat positionChanges(dimension, n);
	positionChanges = PDF.changeInPosition(timeStep);
	double Aij = randu();
	double G1 = PDF.G(positionChanges, timeStep);
	for (int i = 0; i < n; i++) {
		PDF.getParticle(i)->changePosition(positionChanges.col(i));
	}
	double p2 = PDF.PDF();
	double G2 = PDF.G(-positionChanges, timeStep);
	double divident = G2 * p2;
	double divisor = G1 * p1;
	if (Aij <= divident / divisor) {
		p1 = p2;
		vec out(2);
		out(0) = PDF.localEnergy();
		out(1) = PDF.AlphaDerivative();
		return out;
		//return PDF.localEnergy();
	}
	else {
		for (int i = 0; i < n; i++) {
			PDF.getParticle(i)->changePosition(-positionChanges.col(i));
		}
		vec out(2);
		out(0) = PDF.localEnergy();
		out(1) = PDF.AlphaDerivative();
		return out;
		//return PDF.localEnergy();
	}
}

template<class T> 
vec Metropolis<T>::importantSamplingSolve(double timeStep, int iterations) {
	p1 = PDF.PDF();
	vec ret(2);
	double E2 = 0;
	double E = 0;
	vec importantStep;
	for (int i = 0; i < iterations; i++) {
		importantStep = importantSamplingStep(timeStep);
		E += importantStep(0);
		E2 += importantStep(0) * importantStep(0);
	}
	E /= iterations;
	E2 /= iterations;
	ret(0) = E;
	ret(1) = E2 - E * E;
	return ret;
}

template<class T>
void Metropolis<T>::importantSamplingToFile(double timeStep, int iterations) {
	string fileName;
	fileName += "BfImporttimStep"+ to_string(timeStep)+"Dim" + to_string(dimension) + "Npart" + to_string(n) + "Iter" + to_string(iterations) + ".txt";
	ofstream myFile;
	myFile.open(fileName);
	double startAlpha = 0.1;
	double alphaStep = 0.01;
	for (int i = 0; i < 100; i++) {
		double tempAlpha = i * alphaStep + startAlpha;
		PDF.setAlpha(tempAlpha);
		vec out = importantSamplingSolve(timeStep, iterations);
		myFile << tempAlpha << " " << out(0) << " " << out(1) << endl;// " " << E(1) / sqrt(iterations) << endl;
	}
	myFile.close();
}

template<class T>
vec Metropolis<T>::dEnergy(double timeStep) {
	double iterations = 1000;
	//double timeStep = 0.01;
	p1 = PDF.PDF();
	vec ret(2);
	double dE1 = 0;
	double psiAlpha = 0;
	double E = 0;
	vec importantStep;
	for (int i = 0; i < iterations; i++) {
		importantStep = importantSamplingStep(timeStep);
		E += importantStep(0);
		psiAlpha += importantStep(1);
		dE1 += importantStep(0) * importantStep(1);
	}
	E /= iterations;
	psiAlpha /= iterations;
	dE1 /= iterations;
	vec out(2);
	out(0) = E;
	out(1) = 2 * (dE1 - psiAlpha * E);
	return out;
}

template<class T>
vec Metropolis<T>::minimice(double timeStep, double startAlpha, double tol, double stepLength) {
	int MAXITER = 200; //maximal number of iterations
	int iter = 0;
	double alphaMin = startAlpha;
	for (int i = 0; i < MAXITER; i++) {
		PDF.setAlpha(alphaMin);
		vec energies = dEnergy(timeStep);
		double E = energies(0);
		double dE = energies(1);
		//cout << dE << endl;
		if (abs(dE) < tol) {
			vec out(2);
			out(0) = alphaMin;
			out(1) = iter;
			return out;
		}
		iter++;
		alphaMin -= stepLength * dE;
	}
	cout << "Her gjekk d gæli" << endl;
	vec out(2);
	out(0) = -1;
	out(1) = MAXITER;
	return out;
}