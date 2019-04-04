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


double MLWavefunction::MLWave() { 
	double sum1 = 0.0;
	double product1 = 1.0;
	for (int j = 0; j < N; j++) {
		double poversum = 0;
		int acounter = 0;
		for (int i = 0; i < numberOfParticles; i++) {
			vec temp = this->getParticle(i)->getPosition();
			for (int l = 0; l < dimension; l++) {
				poversum += temp[l] * w(acounter, j);
				acounter++;
			}
		}
		product1 *= (1 + exp(b[j] + poversum / sigma));
	}
	int acounter = 0;
	for (int i = 0; i < numberOfParticles; i++) {
		vec temp = this->getParticle(i)->getPosition();
		for (int l = 0; l < dimension; l++) {
			sum1 += (temp[l]  - a[acounter]) *(temp[l] -a[acounter]);
			acounter++;
		}
	}
	return product1 * exp(sum1);
}