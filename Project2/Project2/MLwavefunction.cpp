#include "MLwavefunction.h"

MLWavefunction::MLWavefunction()
{
}

MLWavefunction::~MLWavefunction()
{
}

MLWavefunction::MLWavefunction(int dimension, int numberOfParticles, int numberOfHidenModes) {
	M = dimension * numberOfParticles;
	N = numberOfHidenModes;
	this->dimension = dimension;
	this->numberOfParticles = numberOfParticles;
	setParticles();
	setWeightsAndBiases();
}

void MLWavefunction::setWeightsAndBiases() {
	a = randu(M);
	b = randu(N);
	w = randu(M, N);
}