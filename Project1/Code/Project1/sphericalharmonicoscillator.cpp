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
    int N = this->getNumberOfParticles();
	for (int i = 0; i < N; i++) {
        vec temp = this->getParticle(i)->getPosition();
        double sum = dot(temp,temp);
        product *= exp(-alpha * sum);
	}
	return product;
}


double SphericalHarmonicOscillator::localEnergy() {
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
