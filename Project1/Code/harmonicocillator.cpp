#include "harmonicocillator.h"
#include "wavefunction.h"
#include "particle.h"
#include <vector>

using std::vector;


Harmonicocillator::Harmonicocillator()
{
}


void Harmonicocillator::setPsi() {
	psi.setDimensions(dimension);
	psi.setNumberOfParticels(numberOfParticles);
	psi.setParticels();
}



double Harmonicocillator::harmonicOcillatorWavefunction() {
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


double Harmonicocillator::localEnergy() {
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