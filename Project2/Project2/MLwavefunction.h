#pragma once
#include <armadillo>
#include "wavefunction.h"

using namespace arma;

class MLWavefunction : public Wavefunction 
{
public:
	MLWavefunction();
	~MLWavefunction();
	MLWavefunction(int dimension, int numberOfParticles, int numberOfHidenModes);
	void setWeightsAndBiases();
private:
	int M;  //number of visible/input variables numberOfParticles * dimension
	int N;  //number of hiden variables
	vec a;
	vec b;
	mat w;
	double omega;
};

