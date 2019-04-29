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
	double MLWave();
	double PDF();
	double localEnergy();
	double interactionTerm();
	vec derivativeLogPsiOverA();
	vec derivativeLogPsioverB();
	mat derivativeLogPsioverW();
private:
	int M;  //number of visible/input variables numberOfParticles * dimension
	int N;  //number of hiden variables
	vec a;
	vec b;
	mat w;
	vec x;
	vec h;
	double omega;
	double sigma;
	bool interaction;
};

