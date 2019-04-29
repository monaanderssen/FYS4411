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

void MLWavefunction::setX() {
	int counter = 0;
	for (int i = 0; i < dimension; i++) {
		vec temp = this->getParticle(i)->getPosition();
		for (int j = 0; j < numberOfParticles; j++) {
			x(counter) = temp(j);
			counter++;
		}
	}
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

double MLWavefunction::localEnergy() {
	double sum1 = 0;
	double sum2 = 0;
	for (int i = 0; i < M; i++) {
		double delsum1 = 0;
		for (int j = 0; j < N; j++) {
			double exponent = b(j) + dot(x, w.col(j)) / (sigma*sigma);
			double exponential = exp(exponent);
			delsum1 += w(i, j)*(1 + exp(-exponent));
			sum2 += w(i, j)*w(i, j)*exponential / ((1 + exponential)*(1 + exponential));
		}
		double temp=(x(i) - a(i) + delsum1);
		sum1 += temp * temp;
	}
	double output = 0.5*(M - (sum1 + sum2) / (sigma*sigma)) / (sigma*sigma) + omega * dot(x, x);
	if (interaction) {
		output += interactionTerm();
	}
	return output;
}

double MLWavefunction::interactionTerm() { //not finished
	return 0;
}

double MLWavefunction::derivativeLogPsiOverA(int i) {
	return (x(i) - a(i)) / (sigma*sigma);
}

double MLWavefunction::derivativeLogPsioverB(int j) {
	double exponent = -b(j) - dot(x, w.col(j)) / (sigma*sigma);
	return 1 / (1 + exp(exponent));
}

double MLWavefunction::derivativeLogPsioverW(int i, int j) {
	double exponent = -b(j) - dot(x, w.col(j)) / (sigma*sigma);
	return x(i) / ((1 + exp(exponent))*sigma*sigma);
}

double MLWavefunction::derivativeLogPsiOverAlpha(int i) {
	if (i < M) {
		return derivativeLogPsiOverA(i);
	}
	i -= M;
	if (i < N) {
		return derivativeLogPsioverB(i);
	}
	i -= N;
	//finn resten
}