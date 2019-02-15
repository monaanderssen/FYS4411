#pragma once
#include <vector>
#include <cmath>
#include <armadillo>

using namespace arma;
using std::vector;

class Particle
{
public:
    Particle();
    void setDimension(int dimension) { this->dimension = dimension;}
    int getDimension() { return dimension; }
    void setPosition(vec &position);
    vec getPosition() { return position; }
    void changePosition(vec change);
    int dimension;
private:
    vec position;
};
