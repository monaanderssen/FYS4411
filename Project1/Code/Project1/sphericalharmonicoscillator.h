#pragma once
#include "wavefunction.h"
#include "particle.h"

class SphericalHarmonicOscillator : public Wavefunction
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
    double PDF();
private:
	int numberOfParticles;
	double alpha;
    double omegaHo = 1.;
    int dimension;
	Wavefunction psi;

};


