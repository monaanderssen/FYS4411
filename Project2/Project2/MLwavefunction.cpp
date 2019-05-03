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
	setX();
}

void MLWavefunction::setWeightsAndBiases() {
	a = randn(M)*1;
	b = randn(N)*1;
	w = randn(M, N)*1;
}

void MLWavefunction::setX() {
	int counter = 0;
	x = zeros(M);
	//cout << numberOfParticles << endl;
	for (int i = 0; i < numberOfParticles; i++) {
		vec temp = this->getParticle(i)->getPosition();
		for (int j = 0; j < dimension; j++) {
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
	return product * exp(-sum1*sum1 / (2 * sigma*sigma));
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
			delsum1 += w(i, j)/(1 + exp(-exponent));
			sum2 += w(i, j)*w(i, j)*exponential / ((1 + exponential)*(1 + exponential));
		}
		//cout <<"hei" <<x(i)<<endl;
		double temp=(-(x(i) - a(i)) + delsum1);
		sum1 += temp * temp;
	}
	double output = 0.5*((M - (sum1 + sum2) / (sigma*sigma)) / (sigma*sigma) + omega*omega * dot(x, x));
	if (interaction) {
		output += interactionTerm();
	}
	return output;
}

double MLWavefunction::interactionTerm() { //not finished
	return 0;
}

vec MLWavefunction::derivativeLogPsiOverA() {
	return (x - a) / (sigma*sigma);
}

vec MLWavefunction::derivativeLogPsioverB() {
	vec out(N);
	for (int j = 0; j < N; j++) {
		double exponent = -b(j) - dot(x, w.col(j)) / (sigma*sigma);
		out(j)=1 / (1 + exp(exponent));
	}
	return out;
}
mat MLWavefunction::derivativeLogPsioverW() {
	mat out(M, N);
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			double exponent = -b(j) - dot(x, w.col(j)) / (sigma*sigma);
			out(i,j)=x(i) / ((1 + exp(exponent))*sigma*sigma);
		}
	}
	return out;
}
