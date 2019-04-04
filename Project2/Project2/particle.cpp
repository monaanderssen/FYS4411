#include"particle.h"



Particle::Particle(){
}

void Particle::setPosition(vec &position) {
    this->position = position;
}

void Particle::changePosition(vec change) {
    position += change;
}
