#include <iostream>
#include "metropolis.h"
#include "MLwavefunction.h"
using namespace std;

int main()
{
	MLWavefunction test(2, 2, 2, 1);
	//test.setInteraction();
	Metropolis<MLWavefunction> ttt(test);
	//cout << test.localEnergy();
	//ttt.SGDBruteForce(1, 0.000001, 0.1, 1000000, 200, 10000);
	ttt.SGDImportance(1, 0.001, 0.05, 10000, 400, 1000);
	//test.gibsNewX();
	//ttt.SGDGibbs(.1,10000,400,1000);

}
