#include "harmonicoscillator.h"
#include "wavefunction.h"
#include "particle.h"

HarmonicOscillator::HarmonicOscillator() {

}

HarmonicOscillator::HarmonicOscillator(int dimensions, int numberOfParticles, double a_, double beta_, double gamma_)
{
	a = a_;
	beta = beta_;
	gamma = gamma_;
	this->setDimension(dimensions);
	this->setNumberOfParticles(numberOfParticles);
	this->setParticles();

}
void HarmonicOscillator::setAlpha(double alpha_) {
	alpha = alpha_;
}
double HarmonicOscillator::g(int i) {
	vec x = this->getParticle(i)->getPosition();
	if (dimension == 3) {
		return exp(-alpha * (x[0] * x[0] + x[1] * x[1] + beta * x[2] * x[2]));
	}
	else if (dimension == 2) {
		return exp(-alpha * (x[0] * x[0] + x[1] * x[1]));
	}
	else {
		return exp(-alpha * x[0] * x[0]);
	}
}

double HarmonicOscillator::harmonicWavefunction() {
	double sum1 = 1;
	double sum2 = 0;
	for (int i = 0; i < numberOfParticles; i++) {
		sum1 *= this->g(i);
		for (int j = i + 1; j < numberOfParticles; j++) {
			sum2 += log(f(Norm(i, j)));
		}
	}
	return sum1 * exp(sum2);
}

double HarmonicOscillator::PDF() {
	double x = harmonicWavefunction();
	return x * x;
}

double HarmonicOscillator::Norm(int i, int j) {
	return norm(this->getParticle(i)->getPosition() - this->getParticle(j)->getPosition(), 2);
}


double HarmonicOscillator::V_int(double Norm) {
	if (Norm > a) {
		return 0;
	}
	else {
		return pow(10, 20);
	}
}

double HarmonicOscillator::f(double Norm) {
	if (Norm > a) {
		return 1 - (a / Norm);
	}
	else {
		return 0;
	}
}

double HarmonicOscillator::uDerivative(double Norm) {
	if (Norm > a) {
		return a / (Norm*Norm - a * Norm);
	}
	else {
		return 0;
	}
}

double HarmonicOscillator::uDoubleDerivative(double Norm) {
	if (Norm > a) {
		return (a*a - 2 * a*Norm) / ((Norm*Norm - a * Norm)*(Norm*Norm - a * Norm));
	}
	else {
		return 0;
	}
}

vec HarmonicOscillator::gradientPhi(int i) {
	vec r = this->getParticle(i)->getPosition();
	vec temp(3);
	temp(0) = -2 * alpha*(r(0));
	temp(1) = -2 * alpha*(r(1));
	temp(2) = -2 * alpha*(beta*r(2));
	return temp;
}

double HarmonicOscillator::doubleGradientPhi(int i) {
	vec r = this->getParticle(i)->getPosition();
	return -2 * alpha*((2 + beta) - 2 * alpha*(r(0)*r(0) + r(1)*r(1) + beta * r(2)*r(2)));
}



double HarmonicOscillator::laplacianPsiOverPsi(int k) {
	double sum = doubleGradientPhi(k);
	vec r_k = this->getParticle(k)->getPosition();
	for (int j = 0; j <numberOfParticles; j++) {
		if (j != k) {
			vec r_j = this->getParticle(j)->getPosition();
			sum += 2 * dot(gradientPhi(k), r_k - r_j)*uDerivative(Norm(k, j)) / Norm(k, j);
			sum += (uDoubleDerivative(Norm(k, j)) + 2 * uDerivative(Norm(k, j)) / Norm(k, j));
			for (int i = 0; i<numberOfParticles; i++) {
				if (i != k) {
					vec r_i = this->getParticle(i)->getPosition();
					sum += dot(r_k - r_j, r_k - r_i)*uDerivative(Norm(k, j))*uDerivative(Norm(i, k)) / (Norm(k, j)*Norm(i, k));
				}
			}
		}
	}
	return sum;
}


vec HarmonicOscillator::driftForce(int k) {
	vec sum = zeros(3);
	vec r_k = this->getParticle(k)->getPosition();
	for (int j = 0; j<numberOfParticles; j++) {
		if (j != k) {
			vec r_j = this->getParticle(j)->getPosition();
			sum += (r_k - r_j) / Norm(k, j) *uDerivative(Norm(k, j));
		}
	}
	return 2 * (sum + gradientPhi(k));
}

double HarmonicOscillator::psi() {
	double sum = 0;
	double sum2 = 1;

	for (int i = 0; i <numberOfParticles; i++) {
		vec r_i = this->getParticle(i)->getPosition();
		sum += (r_i(0)*r_i(0) + r_i(1)*r_i(1) + r_i(2)*r_i(2)*beta);
		for (int j = 0; j<i; j++) {
			sum2 *= f(Norm(i, j));
		}
	}
	return exp(-alpha * sum)*sum2;
}

double HarmonicOscillator::localEnergy() {
	double sum = 0;
	double sum2 = 0;

	for (int i = 0; i<numberOfParticles; i++) {
		vec r_i = this->getParticle(i)->getPosition();
		sum += r_i(0)*r_i(0) + r_i(1)*r_i(1) + r_i(2)*r_i(2)*gamma - laplacianPsiOverPsi(i);
		for (int j = 0; j<i; j++) {
			sum2 += V_int(Norm(j, i));
		}
	}
	return 0.5*sum + sum2;
}

double HarmonicOscillator::AlphaDerivative() {
	double sum = 0.0;
	for (int i = 0; i < numberOfParticles; i++) {
		vec temp = getParticle(i)->getPosition();
		sum += (temp(0)*temp(0) + temp(1)*temp(1) + beta * temp(2)*temp(2));
	}
	return -1 * sum;
}

mat HarmonicOscillator::changeInPosition(double timeStep) { //this is aqtualy exactly the samme funktion as for sphericalHarm
	int n = getNumberOfParticles();
	int D = getDimension();
	mat y(D, n); //realy y-x in equations
	y.randn();
	y *= sqrt(timeStep);
	for (int i = 0; i < n; i++) {
		y.col(i) += driftForce(i)*timeStep / 2;
	}
	return y;
}

double HarmonicOscillator::G(mat positionChange, double timeStep) {
	double pi = 3.14159265; //fell free to extend
	int n = getNumberOfParticles();
	double divisor = pow((4 * pi * 0.5 * timeStep), 3.0*n / 2);
	double exponentDivisor = 4 * 0.5*timeStep;
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		vec temp = positionChange.col(i) - driftForce(i)*timeStep / 2;
		sum += dot(temp, temp);
	}
	double outPut = exp(-sum / exponentDivisor) / divisor;
	return outPut;
}







