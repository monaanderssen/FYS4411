#include <iostream>
#include "metropolis.h"
#include "MLwavefunction.h"
using namespace std;

int main()
{
	MLWavefunction test(1, 2, 2, 1);
	Metropolis<MLWavefunction> ttt(test);
	cout << test.localEnergy();
	//ttt.SGDBruteForce(0.0001, 0.000001, 0.001, 100000, 200, 10000);

}
