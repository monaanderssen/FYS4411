#include "MLwavefunction.h"

MLWavefunction::MLWavefunction()
{
}

MLWavefunction::~MLWavefunction()
{
}

MLWavefunction::MLWavefunction(int dimension, int numberOfParticles, int numberOfHidenModes, double sigma) {
	M = dimension * numberOfParticles;
	N = numberOfHidenModes;
	this->sigma = sigma;
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


double MLWavefunction::MLWave() { //with out the partitionfunction 
	double sum1 = norm(x - a);
	double product = 1;
	for (int j = 0; j < N; j++) {
		double exponent = b(j) + dot(x, w.col(j))/(sigma*sigma);
		product *= (1 + exp(exponent));
	}
	return product * exp(sum1*sum1 / (2 * sigma*sigma));
}

double MLWavefunction::PDF() {
	double temp = MLWave();
	return temp * temp;
}

