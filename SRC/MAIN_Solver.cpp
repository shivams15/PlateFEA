#include "Grid.h"
#include "FEMSolver.h"
#include "petsc.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{
	PetscErrorCode ierr;
	ierr = PetscInitialize(&argc, &argv, (char *)0, "Initializing Program");

	int nr, na;
	double w, E, nu;
	int err = 0;

	if(argc < 11){
		printf("Insufficient data provided!\n");
		err++;
	}
	else{
		for(int i=1; i<argc; i++){
			if(!strcmp(argv[i], "-nr"))
				nr = atoi(argv[++i]);
			else if(!strcmp(argv[i], "-na"))
				na = atoi(argv[++i]);
			else if(!strcmp(argv[i], "-w"))
				w = atof(argv[++i]);
			else if(!strcmp(argv[i], "-E"))
				E = atof(argv[++i]);
			else if(!strcmp(argv[i], "-nu"))
				nu = atof(argv[++i]);
			else err++;
		}
	}

	if(nr<1 || na<1 || w<=1 || E<=0 || nu<=-1 || nu>=0.5)
		err++;

	if(err){
		printf("Error processing command line arguments!\nPlease use the following guidelines to provide required data.\n\n");
		printf("-nr	<value>	(positive integer that specifies the number of elements in the radial direction\n");
		printf("-na	<value>	(positive integer that specifies half the number of elements in the angular direction\n");
		printf("-w	<value> (value greater than 1 that specifies the ratio of plate width to the hole radius\n");
		printf("-E	<value> (Young's modulus of the plate material. It must be positive.\n");
		printf("-nu	<value> (Poisson's ratio for the plate material. It must be between -1 and 0.5\n");
		printf("Note: The plate length is automatically set to twice the width\n");
		printf("      All lengths are normalized by the hole radius, and all stresses are normalized by the applied normal stress in x-direction.\n");
	}
	else{
		Grid *grid = new Grid(nr,na, w);
		FEMSolver *fSolver = new FEMSolver(nu, E, grid);
		fSolver->Solve();
	}
	return 0;
}