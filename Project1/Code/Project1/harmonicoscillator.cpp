#include "harmonicoscillator.h"
#include "wavefunction.h"
#include "particle.h"

HarmonicOscillator::HarmonicOscillator(int dimensions,int numberOfParticles, double a_, double beta_, double gamma_)
{
    a = a_;
    beta = beta_;
    gamma = gamma_;
    this->setDimension(dimensions);
    this->setNumberOfParticles(numberOfParticles);
    this->setParticles();

}
void HarmonicOscillator::set_alpha(double alpha_){
    alpha = alpha_;
}
double HarmonicOscillator::g(int i){
    vec x = this->getParticle(i)->getPosition();
    if(dimension == 3){
        return exp(-alpha*(x[0]*x[0] + x[1]*x[1] + beta*x[2]*x[2]));
    }
    else if(dimension == 2){
        return exp(-alpha*(x[0]*x[0] + x[1]*x[1]));
    }
    else{
        return exp(-alpha*x[0]*x[0]);
    }
}

double HarmonicOscillator::harmonicWavefunction(){
    double sum1 = 1;
    double sum2 = 0;
    for(int i = 0; i < numberOfParticles; i++){
        sum1 *= this->g(i);
        for(int j = i + 1; j < numberOfParticles; j++){
            sum2 += log(this->f(i,j));
        }
    }
    return sum1*exp(sum2);
}

double HarmonicOscillator::PDF(){
    double x = harmonicWavefunction();
    return x*x;
}

double HarmonicOscillator::f(int particleOne, int particleTwo){
    return 0.5;
}

double HarmonicOscillator::L1norm(int i, int j){
    return norm(this->getParticle(i)->getPosition() - this->getParticle(j)->getPosition(), 2);
}


double HarmonicOscillator::V_int(double L1norm){
    if(L1norm > a){
        return 0;
    }
}

double HarmonicOscillator::f(double L1norm){
    if(L1norm > a){
        return 1-(a/L1norm);
    }
    else {
        return 0;
    }
}

double HarmonicOscillator::uDerivative(double L1norm){
    if(L1norm > a){
        return a/(L1norm*L1norm - a*L1norm);
    }
    else {
        return 0;
    }
}

double HarmonicOscillator::uDoubleDerivative(double L1norm){
    if(L1norm > a){
        return (a*a - 2*a*L1norm)/((L1norm*L1norm - a*L1norm)*(L1norm*L1norm - a*L1norm));
    }
    else {
        return 0;
    }
}

vec HarmonicOscillator::gradientPhi(int i){
    vec r = this->getParticle(i)->getPosition();
    vec temp(3);
    temp(0) = -2*alpha*(r(0));
    temp(1) = -2*alpha*(r(1));
    temp(2) = -2*alpha*(beta*r(2));
    return temp;
}

double HarmonicOscillator::doubleGradientPhi(int i){
    vec r = this->getParticle(i)->getPosition();
    return -2*alpha*((2 + beta) -2*alpha*(r(0)*r(0) + r(1)*r(1) + beta*r(2)*r(2)));
}



double HarmonicOscillator::laplacianPsiOverPsi(int k){
    double sum = doubleGradientPhi(k);
    vec r_k = this->getParticle(k)->getPosition();
    for (int j=0; j <numberOfParticles;j++) {
        if(j!=k){
        vec r_j = this->getParticle(j)->getPosition();
        double norm = L1norm(k,j);
        sum += 2*dot(gradientPhi(k),r_k-r_j)*uDerivative(norm)/norm;
        sum += (uDoubleDerivative(norm)+2*uDerivative(norm)/norm);
        for (int i=0;i<numberOfParticles;i++) {
            if(i !=k){
            vec r_i = this->getParticle(i)->getPosition();
            double norm2 = L1norm(i,k);
            sum += dot(r_k-r_j,r_k-r_i)*uDerivative(norm)*uDerivative(norm2)/(norm*norm2);
            }
        }
        }
    }
}


vec HarmonicOscillator::driftForce(int k){
    vec sum = zeros(3);
    vec r_k = this->getParticle(k)->getPosition();
    for (int j=0;j<numberOfParticles;j++) {
        if(j != k){
            vec r_j = this->getParticle(j)->getPosition();
            double norm = L1norm(k,j);
            sum+= (r_k-r_j)/norm *uDerivative(norm);
        }
    }
    return 2*(sum + gradientPhi(k));
}

double HarmonicOscillator::psi(){
    double sum = 0;
    double sum2 = 1;

    for (int i = 0;i <numberOfParticles;i++) {
        vec r_i = this->getParticle(i)->getPosition();
        sum += (r_i(0)*r_i(0) + r_i(1)*r_i(1) + r_i(2)*r_i(2)*beta);
        for(int j;j<i; j++){
            sum2 *= f(i,j);
        }
    }
    return exp(-alpha*sum)*sum2;
}












