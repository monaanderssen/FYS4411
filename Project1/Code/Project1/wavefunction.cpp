#include "particle.h"
#include "wavefunction.h"
#include <vector>
#include <random>
#include <cmath>
using std::vector;


Wavefunction::Wavefunction() {

}

Particle Wavefunction::getParticle(int particleNumber) {
    return particles.at(particleNumber);
}

void Wavefunction::setParticles() {
    vector<Particle> temp(numberOfParticles);
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    for (int i = 0; i < numberOfParticles; i++) {
		vector<double> position(dimensions);
		for (int j = 0; j < dimensions; j++) {
			position[j] = distribution(generator);
		}
		temp[i].setPosition(position);
	}
    particles = temp;
}


