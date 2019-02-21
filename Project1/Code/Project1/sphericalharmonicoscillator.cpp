#include "sphericalharmonicoscillator.h"
#include "wavefunction.h"
#include "particle.h"
#include <vector>
#include <armadillo>

using std::vector;
using namespace arma;


SphericalHarmonicOscillator::SphericalHarmonicOscillator()
{
}


void SphericalHarmonicOscillator::setPsi() {
    this->setParticles();
}



double SphericalHarmonicOscillator::harmonicOscillatorWavefunction() {
	double product = 1.0;
	double totalSum = 0.0;
    int N = this->getNumberOfParticles();
	for (int i = 0; i < N; i++) {
        vec temp = this->getParticle(i)->getPosition();
        double sum = dot(temp,temp);
        //product *= exp(-alpha * sum);
		totalSum += sum;
	}
	product *= exp(-alpha * totalSum);
	return product;
}


double SphericalHarmonicOscillator::localEnergyAnalytical() {
	double sum = 0.0;
    int dimensions = this->getDimension();
    int N = this->getNumberOfParticles();
	for (int i = 0; i < N; i++) {
        vec temp = this->getParticle(i)->getPosition();
        sum += dot(temp,temp) * (omegaHo*omegaHo/2 - 2*alpha*alpha);
	}
    sum += dimensions * numberOfParticles *alpha;
	return sum;
}

double SphericalHarmonicOscillator::PDF(){
    double x = harmonicOscillatorWavefunction();
    return x*x;
}


double SphericalHarmonicOscillator::localEnergyNumericalDerivative() {
	double dx = 0.001;
	double dx2 = dx * dx;
	double sum = 0.0;
	int dimentions = this->getDimension();
	int N = this->getNumberOfParticles();
	for (int i = 0; i < N; i++) {
		vec temp = this->getParticle(i)->getPosition();
		double r2i = dot(temp, temp);
		double tempSum = 0.0;
		for (int j = 0; j < dimentions; j++) {
			double positionElement = 2 * temp[j] * dx;
			double diferential = exp(-alpha * (positionElement + dx2)) + exp(-alpha * (dx2 - positionElement)) - 2;
			tempSum += diferential;
		}
		sum += r2i * omegaHo*omegaHo / 2;
		sum -= tempSum / (2*dx2);
	}
	return sum;
}

vec SphericalHarmonicOscillator::driftForce(int i) {
	vec temp = getParticle(i)->getPosition();
	return -4 * alpha*temp;
}

double SphericalHarmonicOscillator::localEnergy() {
	if (analytical) {
		return localEnergyAnalytical();
	}
	else {
		return localEnergyNumericalDerivative();
	}
}