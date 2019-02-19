#pragma once
#include "wavefunction.h"
#include "particle.h"

class SphericalHarmonicOscillator : public Wavefunction
{
public:
    SphericalHarmonicOscillator();
	void setomegaHo(double omegaHo) { this->omegaHo = omegaHo; }
	void setPsi();
	void setAlpha(double alpha) { this->alpha = alpha; }
    double harmonicOscillatorWavefunction();
	double localEnergy();
    double PDF();
	double localEnergyNumericalDerivative(); 
private:
	double alpha;
    double omegaHo;

};


