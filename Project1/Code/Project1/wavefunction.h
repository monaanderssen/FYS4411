#pragma once
#include "particle.h"
using std::vector;

class Wavefunction
{
public:
	Wavefunction();
	void setDimensions(int dimensions) { this->dimensions = dimensions; }
    void setNumberOfParticles(int numberOfParticels) { this->numberOfParticles = numberOfParticles; }
    void setParticles();
    int getNumberOfParticles() { return numberOfParticles; }
	int getDimensions() { return dimensions; }
	Particle getParticle(int particleNumber);
    vector<Particle> getParticels() { return particles; }


private:
	int dimensions;
    int numberOfParticles;
    vector<Particle> particles;

};


