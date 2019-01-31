#pragma once
#include "particle.h"
using std::vector;

class Wavefunction
{
public:
	Wavefunction();
	void setDimensions(int dimensions) { this->dimensions = dimensions; }
	void setNumberOfParticels(int numberOfParticels) { this->numberOfParticels = numberOfParticels; }

	void setParticels();
	int getNumberOfParticles() { return numberOfParticels; }
	int getDimensions() { return dimensions; }
	Particle getParticle(int particleNumber);
	vector<Particle> getParticels() { return particels; }


private:
	int dimensions;
	int numberOfParticels;
	vector<Particle> particels;

};


