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
    double harmonicOscillatorWavefunction();
	double localEnergyAnalytical();
    double PDF();
	double localEnergyNumericalDerivative(); 
	vec driftForce(int i); //return drift force for particle i
	double localEnergy();
	mat changeInPosition(double timeStep); //gives the change of position, don't change the wavefunction.
	double G(mat positionChange, double timeStep);
	double AlphaDerivative();
private:
	double alpha;
    double omegaHo;
	bool analytical = true;
};


