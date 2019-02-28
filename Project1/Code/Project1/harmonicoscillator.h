#pragma once
#include "wavefunction.h"
#include "particle.h"


class HarmonicOscillator : public Wavefunction
{
public:
    HarmonicOscillator(int dimensions, int numberOfParticles, double a_, double beta_, double gamma_);
    void set_alpha(double alpha_);
    double g(int particleNumber);
    double harmonicWavefunction();
    double PDF();
    double norm(int i, int j);
    double V_int(double norm);
    double f(double norm);
    double uDerivative(double norm);
    double uDoubleDerivative(double norm);
    vec gradientPhi(int i);
    double doubleGradientPhi(int i);
    double laplacianPsiOverPsi(int k);
    vec driftForce(int k);
    double psi();
    double localEnergy();

private:
    double a, alpha, beta, gamma;
};
