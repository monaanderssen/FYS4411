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
    double L1norm(int i, int j);
    double V_int(double L1norm);
    double f(double L1norm);
    double uDerivative(double L1norm);
    double uDoubleDerivative(double L1norm);
    double gradient(int i);
    double doubleGradient(int i);


private:
    double a, alpha, beta, gamma;
};
