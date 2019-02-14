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
