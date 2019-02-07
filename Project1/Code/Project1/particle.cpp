#include"particle.h"



Particle::Particle(){
}

void Particle::setPosition(vec &position) {
	this->position = position;
	dimension = position.size();
}

void Particle::changePosition(vec change) {
    position += change;
}
