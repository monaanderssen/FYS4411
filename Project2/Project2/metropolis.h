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
    double importantSamplingStep(double timeStep);
	void SGDBruteForce(double step, double stepLength, int iterations, int MAXITER,int miniBachSize);
	void SGDImportance(double timeStep, double stepLength, int iterations, int MAXITER, int miniBachSize); 
	void SGDGibbs(double stepLength,int iterations, int MAXITER, int miniBachSize);
	double gibsStep();
	T PDF;
    double p1;
    int n;
    int dimension;
};


template <class T>
double Metropolis<T>::bruteForceStep(double step){
    vec r_n(dimension*n);
    r_n.randn();
    r_n *= step;
    double Aij = randu();
    //for(int i = 0; i < n; i++){
        //PDF.getParticle(i)->changePosition(r_n.col(i));
    //}
	PDF.x += r_n;
    double p2 = PDF.PDF();
    if(Aij <= p2/p1){
        p1 = p2;
        return PDF.localEnergy();


    }
    else{
        //for(int i = 0; i < n; i++){
          //  PDF.getParticle(i)->changePosition(-r_n.col(i));
        //}
		PDF.x -= r_n;
        return PDF.localEnergy();


    }
}





template <class T>
double Metropolis<T>::importantSamplingStep(double timeStep) {
    vec positionChanges(dimension* n);
    positionChanges = PDF.changeInPosition(timeStep);
    double Aij = randu();
    double G1 = PDF.G(positionChanges, timeStep);
    PDF.x += positionChanges;
    double p2 = PDF.PDF();
    double G2 = PDF.G(-positionChanges, timeStep);
    double divident = G2 * p2;
    double divisor = G1 * p1;
    if (Aij <= divident / divisor) {
        p1 = p2;
        double out;
        out = PDF.localEnergy();
		//out(1) = PDF.AlphaDerivative();
        return out;
        //return PDF.localEnergy();
    }
    else {
		PDF.x -= positionChanges;
        double out;
        out = PDF.localEnergy();
		//out(1) = PDF.AlphaDerivative();
        return out;
        //return PDF.localEnergy();
    }
}









template<class T>
void Metropolis<T>::SGDBruteForce(double step, double stepLength, int iterations, int MAXITER, int miniBachSize) {
	p1 = PDF.PDF();
	int tempIter = iterations;
	string fileName;
	fileName = "BruteForceInteraction" + to_string(PDF.interaction) + "Dim" + to_string(dimension) + "Npart" + to_string(n) + "Hiden" + to_string(PDF.N) + "Iter" + to_string(iterations) + "LearningRate" + to_string(stepLength) + ".txt";
	ofstream file;
	file.open(fileName);
	double E = 0;
	double E2 = 0;
	int numbMiniBaches = iterations / miniBachSize;
	int M = PDF.getM();
	int N = PDF.getN();
	for (int i = 0; i < MAXITER; i++) {
		int minibach = randi(distr_param(0, numbMiniBaches-1));
		int minibachMin = miniBachSize * minibach;
		int minibachMax = (minibach+ 1)*miniBachSize;
		double minibachE = 0;
		if (i <= (MAXITER - 20)) {
			iterations = 10000;
		}
		else
		{
			iterations = tempIter;
		}
		E = 0;
		E2 = 0;
		vec da = zeros(M);
		vec Eda = zeros(M);
		vec db = zeros(N);
		vec Edb = zeros(N);
		mat dw = zeros(M, N);
		mat Edw = zeros(M, N);
		PDF.setParticles();
		PDF.setX();
		for (int k = 0; k < 10000; k++) {
			bruteForceStep(step); 
		}
		for (int j = 0; j < iterations; j++) {
			double tempE = bruteForceStep(step);
			if (i > (MAXITER - 20)) { //only the last steps in the SGD are writen to file
				file << tempE << " "; 
			}
			E += tempE;
			E2 += tempE * tempE;
			if (minibachMin <= j && j <= minibachMax) {
				
				vec tempa = PDF.derivativeLogPsiOverA();
				vec tempb= PDF.derivativeLogPsioverB();
				mat tempw = PDF.derivativeLogPsioverW();
				minibachE += tempE;
				da += tempa;
				Eda += tempE * tempa;
				db += tempb;
				Edb += tempE * tempb;
				dw += tempw;
				Edw += tempE * tempw;
			}

		}
		vec aGrad = 2 * (Eda / miniBachSize - da * minibachE / (miniBachSize*miniBachSize));
		vec bGrad = 2 * (Edb / miniBachSize - db * minibachE / (miniBachSize*miniBachSize));
		mat wGrad = 2 * (Edw / miniBachSize - dw * minibachE / (miniBachSize*miniBachSize));
		PDF.a -= stepLength * aGrad;
		PDF.b -= stepLength * bGrad;
		PDF.w -= stepLength * wGrad;
		E /= iterations;
		E2 /= iterations;
		cout << E << " " << (E2 - E * E) / sqrt(iterations) << endl;
		if (i > (MAXITER - 20)) {
			file << endl;
		}
	}
	file.close();
}



template<class T>
void Metropolis<T>::SGDImportance(double timeStep, double stepLength, int iterations, int MAXITER, int miniBachSize) {
	p1 = PDF.PDF();
	int tempIter = iterations;
	string fileName;
	fileName = "ImportanceInteraction" + to_string(PDF.interaction) + "Dim" + to_string(dimension) + "Npart" + to_string(n) + "Hiden" + to_string(PDF.N) + "Iter" + to_string(iterations) + "LearningRate" + to_string(stepLength) + ".txt";
	ofstream file;
	file.open(fileName);
	double E = 0;
	double E2 = 0;
	int numbMiniBaches = iterations / miniBachSize;
	//cout << numbMiniBaches << endl;
	int M = PDF.getM();
	int N = PDF.getN();
	for (int i = 0; i < MAXITER; i++) {
		int minibach = randi(distr_param(0, numbMiniBaches - 1));
		int minibachMin = miniBachSize * minibach;
		int minibachMax = (minibach + 1)*miniBachSize;
		if (i <= (MAXITER - 20)) {
			iterations = 10000;
		}
		else iterations = tempIter;
		double minibachE = 0;
		E = 0;
		E2 = 0;
		vec da = zeros(M);
		vec Eda = zeros(M);
		vec db = zeros(N);
		vec Edb = zeros(N);
		mat dw = zeros(M, N);
		mat Edw = zeros(M, N);
		PDF.setParticles();
		PDF.setX();
		for (int k = 0; k < 1000; k++) {
			importantSamplingStep(timeStep);
		}
		for (int j = 0; j < iterations; j++) {
			double tempE = importantSamplingStep(timeStep);
			if (i > (MAXITER - 20)) {
				file << tempE << " ";
			}
			E += tempE;
			E2 += tempE * tempE;
			if (minibachMin <= j && j <= minibachMax) {

				vec tempa = PDF.derivativeLogPsiOverA();
				vec tempb = PDF.derivativeLogPsioverB();
				mat tempw = PDF.derivativeLogPsioverW();
				minibachE += tempE;
				da += tempa;
				Eda += tempE * tempa;
				db += tempb;
				Edb += tempE * tempb;
				dw += tempw;
				Edw += tempE * tempw;
			}

		}
		vec aGrad = 2 * (Eda / miniBachSize - da * minibachE / (miniBachSize*miniBachSize));
		vec bGrad = 2 * (Edb / miniBachSize - db * minibachE / (miniBachSize*miniBachSize));
		mat wGrad = 2 * (Edw / miniBachSize - dw * minibachE / (miniBachSize*miniBachSize));
		PDF.a -= stepLength * aGrad;
		PDF.b -= stepLength * bGrad;
		PDF.w -= stepLength * wGrad;
		E /= iterations;
		E2 /= iterations;
		cout << E << " " << (E2 - E * E) / sqrt(iterations) << endl;
		if (i > (MAXITER - 20)) {
			file << endl;
		}
	}
}

template<class T>
double Metropolis<T>::gibsStep() {
	PDF.gibsNewX();
	PDF.gibsNewH();
	return PDF.gibsLocalEnergy();
}

template<class T>
void Metropolis<T>::SGDGibbs(double stepLength, int iterations, int MAXITER, int miniBachSize) {
	p1 = PDF.PDF();
	int tempIter = iterations;
	string fileName;
	fileName = "GibbsInteraction" + to_string(PDF.interaction) + "Dim" + to_string(dimension) + "Npart" + to_string(n) + "Hiden" + to_string(PDF.N) + "Iter" + to_string(iterations) + "LearningRate" + to_string(stepLength) +"Sigma"+to_string(PDF.sigma)+ ".txt";
	ofstream file;
	file.open(fileName);
	double E = 0;
	double E2 = 0;
	int numbMiniBaches = iterations / miniBachSize;
	//cout << numbMiniBaches << endl;
	int M = PDF.getM();
	int N = PDF.getN();
	for (int i = 0; i < MAXITER; i++) {
		int minibach = randi(distr_param(0, numbMiniBaches - 1));
		int minibachMin = miniBachSize * minibach;
		int minibachMax = (minibach + 1)*miniBachSize;
		if (i <= (MAXITER - 20)) {
			iterations = 100000;
		}
		else iterations = tempIter;
		double minibachE = 0;
		E = 0;
		E2 = 0;
		vec da = zeros(M);
		vec Eda = zeros(M);
		vec db = zeros(N);
		vec Edb = zeros(N);
		mat dw = zeros(M, N);
		mat Edw = zeros(M, N);
		PDF.setParticles();
		PDF.setX();
		for (int j = 0; j < iterations; j++) {
			double tempE = gibsStep();
			if (i > (MAXITER - 20)) {
				file << tempE << " ";
			}
			E += tempE;
			E2 += tempE * tempE;
			if (minibachMin <= j && j <= minibachMax) {

				vec tempa = PDF.gibsDerivativeLogPsiOverA();
				vec tempb = PDF.gibsDerivativeLogPsioverB();
				mat tempw = PDF.gibsDerivativeLogPsioverW();
				minibachE += tempE;
				da += tempa;
				Eda += tempE * tempa;
				db += tempb;
				Edb += tempE * tempb;
				dw += tempw;
				Edw += tempE * tempw;
			}

		}
		vec aGrad = 2 * (Eda / miniBachSize - da * minibachE / (miniBachSize*miniBachSize));
		vec bGrad = 2 * (Edb / miniBachSize - db * minibachE / (miniBachSize*miniBachSize));
		mat wGrad = 2 * (Edw / miniBachSize - dw * minibachE / (miniBachSize*miniBachSize));
		PDF.a -= stepLength * aGrad;
		PDF.b -= stepLength * bGrad;
		PDF.w -= stepLength * wGrad;
		E /= iterations;
		E2 /= iterations;
		cout << E << " " << (E2 - E * E) / sqrt(iterations) << endl;
		if (i > (MAXITER - 20)) {
			file << endl;
		}
	}
}