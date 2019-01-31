#include "particle.h"
#include "wavefunction.h"
#include <vector>
#include <random>

using std::vector;


Wavefunction::Wavefunction() {

}

Particle Wavefunction::getParticle(int particleNumber) {
	return particels.at(particleNumber);
}

void Wavefunction::setParticels() {
	vector<Particle> temp(numberOfParticels);
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(-1.0, 1.0);
	for (int i = 0; i < numberOfParticels; i++) {
		vector<double> position(dimensions);
		for (int j = 0; j < dimensions; j++) {
			position[j] = distribution(generator);
		}
		temp[i].setPosition(position);
	}
	particels = temp;
}