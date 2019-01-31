#include "harmonicocillator.h"
#include "wavefunction.h"

double harmonicocillator::harmonicOcillatorWavefunction() {
	double product = 1;
	for (int i = 0; i < psi.getNumberOfParticles(); i++) {
		double sum = 0;
		vector<double> temp = psi.getParticle(i).getPosition();
		for (int j = 0; j < psi.getDimensions(); j++) {
			sum += temp[j] * temp[j];
		}
		product *= exp(-alpha * sum);
	}
	return product;
}