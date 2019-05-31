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
	void setInteraction() { interaction = true; }//only run this if you want interaction
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
	vec a; //Visible biases
	vec b; //hiden biases
	mat w; //weights
	vec x; //Visible nodes
	vec h; //Hiden nodes
	double omega=1.0; //in this project we are only considering omega=1
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

