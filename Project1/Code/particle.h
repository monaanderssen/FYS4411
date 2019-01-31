#pragma once
#include <vector>
using std::vector;

class Particle
{
public:
	Particle();
	void setDimension(int dimension) { this->dimension = dimension;}
	int getDimension() { return dimension; }
	void setPosition(vector<double> &position);
	vector<double> getPosition() { return position; }
	void changePosition(double change, int cordinate); 
private:
	int dimension;
	vector<double> position;
};



