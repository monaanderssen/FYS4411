#pragma once
#include "wavefunction.h"
#include "particle.h"

class SphericalHarmonicOscillator
{
public:
    SphericalHarmonicOscillator();
	void setomegaHo(double omegaHo) { this->omegaHo = omegaHo; }
	void setDimension(double dimension) { this->dimension = dimension; }
	void setPsi();
	void setAlpha(double alpha) { this->alpha = alpha; }
	void setNumberOfParticles(int numberOfParticles) { this->numberOfParticles = numberOfParticles; }
    double harmonicOscillatorWavefunction();
	double localEnergy();
private:
	int numberOfParticles;
	double alpha;
	double omegaHo;
	double dimension;
	Wavefunction psi;

};


