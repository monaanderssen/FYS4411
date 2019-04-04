#include "particle.h"
#include "wavefunction.h"
#include <vector>
#include <random>
#include <cmath>
using std::vector;


Wavefunction::Wavefunction() {

}

Particle * Wavefunction::getParticle(int particleNumber) {
    Particle * p = &particles.at(particleNumber);
    return p;
    //return particles.at(particleNumber);
}

void Wavefunction::setParticles() {
    vector<Particle> temp(numberOfParticles);
    for (int i = 0; i < numberOfParticles; i++) {
        vec position(dimension);
        position.randn();
        position*= 0.1;
        temp[i].setPosition(position);
    }
    particles = temp;
}


