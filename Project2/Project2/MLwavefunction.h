#pragma once
#include <armadillo>
#include "wavefunction.h"

using namespace arma;

class MLWavefunction : public Wavefunction 
{
public:
	MLWavefunction();
	~MLWavefunction();
	MLWavefunction(int dimension, int numberOfParticles, int numberOfHidenModes, double sigma);
	void setWeightsAndBiases();
	void setX();
	double MLWave();
	double PDF();
	double localEnergy();
	double interactionTerm();
	int getM() { return M; }
	int getN() { return N; }
	vec derivativeLogPsiOverA();
	vec derivativeLogPsioverB();
	mat derivativeLogPsioverW();
	int M;  //number of visible/input variables numberOfParticles * dimension
	int N;  //number of hiden variables
	vec a;
	vec b;
	mat w;
	vec x;
	vec h;
	double omega=1.0;
	double sigma;
	bool interaction;
	vec driftForce();
	vec changeInPosition(double timeStep);
	double G(vec positionChanges, double timeStep);
	double gibsWavefunction();
	double gibsPDF();
	vec gibsDerivativeLogPsiOverA();
	vec gibsDerivativeLogPsioverB();
	mat gibsDerivativeLogPsioverW();
	double gibsLocalEnergy();
	void gibsNewX();
	void gibsNewH();
};

