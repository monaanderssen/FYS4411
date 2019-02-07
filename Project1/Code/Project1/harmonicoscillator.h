#pragma once
#include "wavefunction.h"
#include "particle.h"


class HarmonicOscillator : public Wavefunction
{
public:
    HarmonicOscillator(int dimensions, int numberOfParticles, double a_, double beta_, double gamma_);
    void set_alpha(double alpha_);
    double g(int particleNumber);
    double f(int particleOne, int particleTwo);
    double harmonicWavefunction();
    double PDF();


private:
    double a, alpha, beta, gamma;
};

