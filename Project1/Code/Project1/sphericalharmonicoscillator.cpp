#include "sphericalharmonicoscillator.h"
#include "wavefunction.h"
#include "particle.h"
#include <vector>
#include <cmath>

using std::vector;


SphericalHarmonicOscillator::SphericalHarmonicOscillator()
{
}


void SphericalHarmonicOscillator::setPsi() {
	psi.setDimensions(dimension);
    psi.setNumberOfParticles(numberOfParticles);
    psi.setParticles();
}



double SphericalHarmonicOscillator::harmonicOscillatorWavefunction() {
	double product = 1.0;
	int dimensions = psi.getDimensions();
	int N = psi.getNumberOfParticles();
	for (int i = 0; i < N; i++) {
		double sum = 0.0;
		vector<double> temp = psi.getParticle(i).getPosition();
		for (int j = 0; j < dimensions; j++) {
			sum += temp[j] * temp[j];
		}
		product *= exp(-alpha * sum);
	}
	return product;
}


double SphericalHarmonicOscillator::localEnergy() {
	double sum = 0.0;
	int dimensions = psi.getDimensions();
	int N = psi.getNumberOfParticles();
	for (int i = 0; i < N; i++) {
		vector<double> temp = psi.getParticle(i).getPosition();
		double tempSum = 0.0;
		for (int j = 0; j < dimensions; j++) {
			tempSum += temp[j] * temp[j];
		}
		sum += tempSum * (omegaHo / 2 -2*alpha*alpha);
	}
	sum += dimensions * numberOfParticles *alpha;
	return sum;
}
