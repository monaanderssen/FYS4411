#include <iostream>
#include "metropolis.h"
#include "MLwavefunction.h"
using namespace std;

int main()
{
	MLWavefunction test(2, 2, 2, sqrt(0.5)); //Dimension, number of particles, hidden nodes, sigma.
	//test.setInteraction();  //Sett in tis to include interactions.
	Metropolis<MLWavefunction> ttt(test);
	ttt.SGDBruteForce(1, 0.1, 1<<20, 200, 10000); //the parameters are steplength, learningrate, iterations, lerningcycles, bachSize
	ttt.SGDImportance(1, 0.1, 1<<20, 200, 10000); // ------------------------""-------------------
	ttt.SGDGibbs(.1,1<<20,1000,10000); //The parameters are learningrate, iterations, lerningcycles, bachSize

}
