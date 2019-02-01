#include"particle.h"

using std::vector;

Particle::Particle(){
}

void Particle::setPosition(vector<double> &position) {
	this->position = position;
	dimension = position.size();
}

void Particle::changePosition(double change, int cordinate) {
	position.at(cordinate) += change; 
}