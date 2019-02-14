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
    Particle * getParticle(int particleNumber);
    vector<Particle> getParticles() { return particles; }
    int numberOfParticles;
private:
    vector<Particle> particles;

};




