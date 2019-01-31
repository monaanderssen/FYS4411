#pragma once
#include "wavefunction.h"

class harmonicocillator
{
public:
	harmonicocillator();
	void setomegaHo(double omegaHo) { this->omegaHo = omegaHo; };
	void setPsi();
	void setAlpha(double alpha) { this->alpha = alpha; };
private:
	double alpha;
	double omegaHo;
	Wavefunction psi;
	double harmonicOcillatorWavefunction();

};

harmonicocillator::harmonicocillator()
{
}

