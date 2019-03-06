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
    double L1norm(int i, int j);
    double V_int(double L1norm);
    double f(double L1norm);
    double uDerivative(double L1norm);
    double uDoubleDerivative(double L1norm);
    vec gradientPhi(int i);
    double doubleGradientPhi(int i);
    double laplacianPsiOverPsi(int k);
    vec driftForce(int k);
    double psi();
    double localEnergy();

private:
    double a, alpha, beta, gamma;
};
