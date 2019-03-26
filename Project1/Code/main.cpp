#include <iostream>
#include "sphericalharmonicoscillator.h"
#include "metropolis.h"
#include <armadillo>
#include <mpi.h>
#include <vector>
#include "harmonicoscillator.h"

using namespace std;
using namespace arma;

int main()
{

    /*

    //Write energies and variance to file for different alphas using openMPI

    int num_alphas = 100;
    int iterations = 131072;
    int dimension = 1;
    int npart = 100;
    double stepsize = 0.1;

    vec alpha = linspace(0.1,1,num_alphas);


    int my_rank, numprocs, idum;
    MPI_Init (NULL,NULL);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    ofstream myFile;
    string filename;
    filename += "ISSPD" + to_string(dimension) + "N" + to_string(npart) + "iters" + to_string(iterations) + "proc" + to_string(my_rank) + "step" + to_string(stepsize) + ".txt";
    myFile.open(filename);

    int no_intervalls = num_alphas/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;
    if ( (my_rank == numprocs-1) &&( myloop_end < num_alphas) ) myloop_end = num_alphas;
    idum = -1-my_rank;
    arma::arma_rng::set_seed(idum);

    for(int cycles = myloop_begin; cycles <= myloop_end; cycles++){
        SphericalHarmonicOscillator spHa;
        spHa.setAlpha(alpha(cycles-1));
        spHa.setDimension(dimension);
        spHa.setomegaHo(1);
        spHa.setParticleSpacing(0.3);
        spHa.setNumberOfParticles(npart);
        spHa.setPsi();
        Metropolis<SphericalHarmonicOscillator> bruteF(spHa);
        vec En = bruteF.importantSamplingSolve(stepsize,iterations);
        myFile << to_string(En(0)) + " " + to_string(En(1)) + " " + to_string(alpha(cycles - 1)) << endl;

    }

    myFile.close();

    MPI_Finalize();

    */

    /*
    //Write local energy from every step to file (for grid of alphas) in order to use blocking.

    int num_alphas = 100;
    int iterations = 2048;
    int dimension = 3;
    int npart = 100;
    double stepsize = 0.004;

    vec alpha = linspace(0.1,1,num_alphas);

    int my_rank, numprocs, idum;
    MPI_Init (NULL,NULL);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    int no_intervalls = num_alphas/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;
    if ( (my_rank == numprocs-1) &&( myloop_end < num_alphas) ) myloop_end = num_alphas;
    idum = -1-my_rank;
    arma::arma_rng::set_seed(idum);

    for(int cycles = myloop_begin; cycles <= myloop_end; cycles++){
        HarmonicOscillator spHa(dimension, npart, 0.0043, 2.82843, 2.82843, 0.5);
        spHa.setAlpha(alpha(cycles-1));
        Metropolis<HarmonicOscillator> bruteF(spHa);
        bruteF.importantSamplingToFile(stepsize,iterations);
    }

    MPI_Finalize();
    */


    /*
    //Do gradient descent to find optimal alpha
    int dimensions = 3;
    int npart = 100;
    double stepsize = 0.000001;
    double timestep = 0.004;
    MPI_Init(NULL,NULL);
    HarmonicOscillator HarOsc(dimensions, npart, 0.0043, 2.82843, 2.82843,0.5);
    Metropolis<HarmonicOscillator> gradDesc(HarOsc);
    vec E = gradDesc.minimize(timestep,0.4835465,1e-8,stepsize,10000,200);
    MPI_Finalize();
    cout << E(0) << " " << E(1) << endl;
    */

    /*
    //Write the local energy from every step to file for single alpha.
    int iterations = 1048576;
    int dimension = 3;
    int npart = 100;
    double stepsize = 0.004;

    HarmonicOscillator spHa(dimension, npart, 0.0043, 2.82843, 2.82843, 0.5);
    spHa.setAlpha(0.4835465);
    Metropolis<HarmonicOscillator> bruteF(spHa);
    bruteF.importantSamplingToFile(stepsize,iterations);
    */
    return 0;
}
