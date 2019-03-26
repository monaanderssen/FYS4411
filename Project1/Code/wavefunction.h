#pragma once
#include "particle.h"
using std::vector;

class Wavefunction : public Particle
{
public:
    Wavefunction();
    void setNumberOfParticles(int numberOfParticles) { this->numberOfParticles = numberOfParticles; }
    void setParticles();
    int getNumberOfParticles() { return numberOfParticles; }
    void particleSpacing(double variance_);
    Particle * getParticle(int particleNumber);
    vector<Particle> getParticles() { return particles; }
    int numberOfParticles;
    double variance;
private:
    vector<Particle> particles;

};




