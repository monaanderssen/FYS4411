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
	double localEnergyAnalytical(); //local energy with analytical derivatives
    double PDF();
	double localEnergyNumericalDerivative(); //local energy with numerical derivatives
	vec driftForce(int i); //return drift force for particle i
	double localEnergy();
	mat changeInPosition(double timeStep); //gives the change of position, don't change the wavefunction.
	double G(mat positionChange, double timeStep);
	double AlphaDerivative(); //derivative of the wavefunction with respect to alpha 
    double alpha;
private:

    double omegaHo;
	bool analytical = true; //desides if localEnergy uses numerical or analytical derivatives
};


