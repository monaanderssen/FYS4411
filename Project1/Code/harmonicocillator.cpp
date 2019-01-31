#include "harmonicocillator.h"
#include "wavefunction.h"
#include <vector>

using std::vector;

void harmonicocillator::setPsi() {
	psi.setDimensions(dimension);
	psi.setNumberOfParticels(numberOfParticles);
	psi.setParticels();
}



double harmonicocillator::harmonicOcillatorWavefunction() {
	double product = 1.0;
	int dimension = psi.getDimensions();
	int N = psi.getNumberOfParticles();
	for (int i = 0; i < N; i++) {
		double sum = 0.0;
		vector<double> temp = psi.getParticle(i).getPosition();
		for (int j = 0; j < dimension; j++) {
			sum += temp[j] * temp[j];
		}
		product *= exp(-alpha * sum);
	}
	return product;
}


double harmonicocillator::localEnergy() {
	double sum = 0.0;
	int dimension = psi.getDimensions();
	int N = psi.getNumberOfParticles();
	for (int i = 0; i < N; i++) {
		vector<double> temp = psi.getParticle(i).getPosition();
		double tempSum = 0.0;
		for (int j = 0; j < dimension; j++) {
			tempSum += temp[j] * temp[j];
		}
		sum += tempSum * (omegaHo / 2 -2*alpha*alpha);
	}
	sum += dimension * numberOfParticles *alpha;
	return sum;
}