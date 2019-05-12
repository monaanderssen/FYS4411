#include <iostream>
#include "metropolis.h"
#include "MLwavefunction.h"
using namespace std;

int main()
{
	MLWavefunction test(1, 2, 4, 1);
	Metropolis<MLWavefunction> ttt(test);
	//cout << test.localEnergy();
	//ttt.SGDBruteForce(1, 0.000001, 0.01, 1000000, 200, 10000);
	//ttt.SGDImportance(0.001, 0.001, 0.01, 1000000, 50, 10000);
	test.gibsNewX();
}
