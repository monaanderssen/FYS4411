#pragma once
#include "wavefunction.h"

class harmonicocillator
{
public:
	harmonicocillator();
	void setomegaHo(double omegaHo) { this->omegaHo = omegaHo; };
	void setDimension(double dimension) { this->dimension = dimension; }
	void setPsi();
	void setAlpha(double alpha) { this->alpha = alpha; };
	void setNumberOfParticles(int numberOfParticles) { this->numberOfParticles; }
	double harmonicOcillatorWavefunction();
	double localEnergy();
private:
	int numberOfParticles;
	double alpha;
	double omegaHo;
	double dimension;
	Wavefunction psi;

};

harmonicocillator::harmonicocillator()
{
}

