#pragma once
#include "wavefunction.h"
#include "particle.h"
#include <cstdlib>
#include <random>
#include <armadillo>
#include <iostream>
#include <fstream>
//#include <mpi.h> //har problemer med � kompilere denne

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
    void bruteForceTofile(double step, int iterations, int my_rank); //runs over diferent values of alpha
    vec importantSamplingStep(double timeStep);
    vec importantSamplingSolve(double timeStep, int iterations); //finds only localEnergy
    void importantSamplingToFile(double timeStep, int inerations);
    vec dEnergy(double timeStep,int iterations); //returns a vector with energy and derivative of energy using importantSampling
    vec minimize(double timeStep,double startAlpha,double tol, double stepLength, int iterations=10000, int MAXITER=200); //finding minimum using stepest desent
	vec SGDBruteForce(double step, double tol, double stepLength, int iterations, int MAXITER,int miniBachSize);
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
void Metropolis<T>::bruteForceTofile(double step, int iterations,int my_rank) {
    double E;
    p1 = PDF.PDF();
    for(int i = 0; i < 50000; i++){
        E = bruteForceStep(step);
    }
    string fileName;
    fileName += "BfMcDim" + to_string(dimension)+ "step" + to_string(step) + "Npart"+to_string(n) + "Iter" + to_string(iterations) + "alpha" + to_string(PDF.alpha) + ".txt";
    ofstream myFile;
    myFile.open(fileName);
    for (int i = 0; i < iterations; i++) {
        E = bruteForceStep(step);
        myFile << E << endl;
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
    for(int i = 0; i < 1000; i++){
        importantStep = importantSamplingStep(timeStep);
    }
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
    p1 = PDF.PDF();
    vec importantStep;
    for(int i = 0; i < 1000; i++){ //Let method converge
        importantStep = importantSamplingStep(timeStep);
    }
    string fileName;
    fileName += "BfImporttimStep"+ to_string(timeStep)+"Dim" + to_string(dimension) + "Npart" + to_string(n) + "Iter" + to_string(iterations) + "alpha" + to_string(PDF.alpha) + ".txt";
    ofstream myFile;
    myFile.open(fileName);
    for (int i = 0; i < iterations; i++) {
        importantStep = importantSamplingStep(timeStep);
        myFile << importantStep(0) << endl;
    }
    myFile.close();
}

template<class T>
vec Metropolis<T>::dEnergy(double timeStep, int iterations) {
    //double timeStep = 0.01;
    vec ret(2);
    double dE1 = 0;
    double psiAlpha = 0;
    double E = 0;
    double Etot = 0;
    double psiAlphatot = 0;
    double dE1tot = 0;
    int my_rank, numprocs, idum;

    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    int no_intervalls = iterations/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;
    if ( (my_rank == numprocs-1) &&( myloop_end < iterations) ) myloop_end = iterations;
    idum = -1-my_rank;
    arma::arma_rng::set_seed(idum);

    p1 = PDF.PDF();
    vec importantStep;

    for(int cycles = myloop_begin; cycles <= myloop_end; cycles++){
        importantStep = importantSamplingStep(timeStep);
        E += importantStep(0);
        psiAlpha += importantStep(1);
        dE1 += importantStep(0)*importantStep(1);
    }

    MPI_Allreduce(&E, &Etot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&psiAlpha, &psiAlphatot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&dE1, &dE1tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    Etot /= iterations;
    psiAlphatot /= iterations;
    dE1tot /= iterations;

    vec out(2);
    out(0) = Etot;
    out(1) = 2 * (dE1tot - psiAlphatot * Etot);
    return out;
}

template<class T>
vec Metropolis<T>::minimize(double timeStep, double startAlpha, double tol, double stepLength,int iterations, int MAXITER) {
    int iter = 0;
    double alphaMin = startAlpha;
    double dE;
    vec importantStep;
    PDF.setAlpha(alphaMin);
    p1 = PDF.PDF();
    for(int i = 0; i < 10000; i++){
        importantStep = importantSamplingStep(timeStep);
    }
    for (int i = 0; i < MAXITER; i++) {
        PDF.setAlpha(alphaMin);
        vec energies = dEnergy(timeStep, iterations);
        double E = energies(0);
        dE = energies(1);
        //cout << dE << endl;
        if (abs(dE) < tol) {

            vec out(2);
            out(0) = alphaMin;
            out(1) = iter;
            return out;
        }
        iter++;
        alphaMin -= stepLength * dE;
        cout << alphaMin << endl;
    }
    cout << "WARNING: did not converge: dE = " << dE <<  endl;
    vec out(2);
    out(0) = -1;
    out(1) = MAXITER;
    return out;
}

template<class T>
vec Metropolis<T>::SGDBruteForce(double step, double tol, double stepLength, int iterations, int MAXITER, int miniBachSize) {
	p1 = PDF.PDF();
	double E = 0;
	double E2 = 0;
	int numbMiniBaches = iterations / miniBachSize;
	int M = PDF.getM();
	int N = PDF.getN();
	for (int i = 0; i < MAXITER; i++) {
		int minibach = randi(distr_param(0, numbMiniBaches-1));
		int minibachMin = miniBachSize * minibach;
		int minibachMax = (minibachMin + 1)*minibach;
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
			bruteForceStep(step);
		}
		for (int j = 0; j < iterations; j++) {
			double tempE = bruteForceStep(step);
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
	}
}
