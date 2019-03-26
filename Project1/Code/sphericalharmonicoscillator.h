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
	void setAnalytical(bool analytical) { this->analytical = analytical; }
    void setParticleSpacing(double variance_ = 0.1){this->particleSpacing(variance_); }
    double harmonicOscillatorWavefunction();
	double localEnergyAnalytical();
    double PDF();
	double localEnergyNumericalDerivative(); 
	vec driftForce(int i); //return drift force for particle i
	double localEnergy();
	mat changeInPosition(double timeStep); //gives the change of position, don't change the wavefunction.
	double G(mat positionChange, double timeStep);
	double AlphaDerivative();
    double alpha;
private:

    double omegaHo;
	bool analytical = true;
};


