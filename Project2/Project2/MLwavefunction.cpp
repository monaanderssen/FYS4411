#include "MLwavefunction.h"
#include <math.h>
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
	a = randn(M)*0.3;
	b = randn(N)*0.3;
	w = randn(M, N)*0.3;
	h = zeros(N);
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


vec MLWavefunction::driftForce() {
	vec out = -(x - a);
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			double exponent = b(j) + dot(x, w.col(j)) / (sigma*sigma);
			out(i) += w(i, j) / (1 + exp(-exponent));
		}
	}
	return 2*out /(sigma*sigma);
}

vec MLWavefunction::changeInPosition(double timeStep) {
	return 0.5*driftForce() * timeStep + randn(M)*sqrt(timeStep);
}

double MLWavefunction::G(vec positionChanges, double timeStep) {
	double pi = 3.14159265; //fell free to extend
	int n = getNumberOfParticles();
	double divisor = pow((4 * pi * 0.5 * timeStep), 3.0*n / 2);
	double exponentDivisor = 4 * 0.5*timeStep;
	vec temp = positionChanges - 0.5*timeStep*driftForce();
	double out = exp(-dot(temp, temp) / exponentDivisor) / divisor;
	return out;
}

double MLWavefunction::gibsWavefunction() {
	return sqrt(MLWave());
}

double MLWavefunction::gibsPDF() {
	return MLWave();
}

vec MLWavefunction::gibsDerivativeLogPsiOverA() {
	return 0.5 *derivativeLogPsiOverA();
}

vec MLWavefunction::gibsDerivativeLogPsioverB() {
	return 0.5*derivativeLogPsioverB();
}

mat MLWavefunction::gibsDerivativeLogPsioverW() {
	return 0.5 *derivativeLogPsioverW();
}

double MLWavefunction::gibsLocalEnergy() {
	double sum1 = 0;
	double sum2 = 0;
	for (int i = 0; i < M; i++) {
		double delsum1 = 0;
		for (int j = 0; j < N; j++) {
			double exponent = b(j) + dot(x, w.col(j)) / (sigma*sigma);
			double exponential = exp(exponent);
			delsum1 += w(i, j) / (1 + exp(-exponent));
			sum2 += w(i, j)*w(i, j)*exponential / ((1 + exponential)*(1 + exponential));
		}
		//cout <<"hei" <<x(i)<<endl;
		double temp = (-(x(i) - a(i)) + delsum1);
		sum1 += temp * temp;
	}
	double output = 0.5*((M*0.5 - (0.25*sum1 + 0.5*sum2) / (sigma*sigma)) / (sigma*sigma) + omega * omega * dot(x, x));
	if (interaction) {
		output += interactionTerm();
	}
	return output;
}

void MLWavefunction::gibsNewX() {
	x = mvnrnd(a + w * h, eye(M,M)*sigma*sigma);
}

void MLWavefunction::gibsNewH() {
	vec prob = randu(N);
	for (int i = 0; i < N; i++) {
		double p = 1 / (1 + exp(b(i) + dot(x, w.col(i))/(sigma*sigma)));
		if (prob(i) < p) {
			h(i) = 0;
		}
		else h(i) = 1;
	}
}