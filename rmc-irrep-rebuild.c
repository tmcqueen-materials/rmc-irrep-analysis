#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdint.h>
#include <math.h>
#include "defs.h"
#include "globals.h"
#define USAGE_MESSAGE "Script Use: rmc-irrep-rebuild [xyzinput].cfg [irrepinput].cfg [ircontrib].irrep\n"
#include "io.h"
#include "memlapack.h"
#include "supercell.h"
#include "xyzinput.h"
#include "irrepinput.h"
#include "ircontrib.h"
#include "eyefile.h"
#include "irreprebuild.h"

int main(int argc, char *argv[]) { // argv[1] = xyzinput, argv[2] = irrepinput, argv[3] = ircontrib
	char *opfName = NULL;

	// results
	LAPACK *atomMatrix = NULL; // numAtoms x 3 matrix (real)
	int *atomMatrixLabels = NULL; // numAtoms matrix (real)

	// check number of input parameters
	if (argc != 4) printUsageAndExit(1);

        // Disable output buffering of status and error messages (helps debugging)
        setvbuf(stdout, NULL, _IONBF, 0);
        setvbuf(stderr, NULL, _IONBF, 0);

	/*
	 * Read ideal unit cell, the "xyzinput"
	 */
	readXYZFile(FTA_FILE, argv[1]);

	/*
	 * Read in irrep coefficients.
	 */
	readIRContribFile(FTA_FILE,argv[3]);

        /*
         * Read in irreducible representations.
         */
        readIrrepFile(FTA_FILE,argv[2]);

	// compute number of atoms
	numAtoms = numSites*numCells;

	// Allocate result matrices
	printStatus("Allocating result matrices...\n");
        atomMatrix = lp_malloc(numAtoms,3,false);
        atomMatrixLabels = safe_malloc((UINT64)numAtoms*(UINT64)sizeof(int));

	/*
	 * Compute supercell atom positions
	 */
        printStatus("Generating supercell atom positions...\n");
	irrepRebuild(1.0, irrepCoeffs, &atomMatrix, &atomMatrixLabels);

	/*
	 * Write result to file
	 */
	opfName = safe_malloc((UINT64)strlen(argv[3])+(UINT64)5);
	strncpy(opfName, argv[3], strlen(argv[3])+5);
	strncat(opfName, "_eye", 4);
	writeEyeFile(opfName, atomMatrix, atomMatrixLabels);
	safe_free(opfName); opfName = NULL;

	printStatus("Memory Used (excluding stack and overhead): %i MB\n", maxMemAllocated/1048576);
	
	// Free memory 
	safe_free(atomMatrixLabels);
	lp_free(atomMatrix);
	freeGlobals();

	if (memAllocated > (UINT64)0)
		printStatus("%i bytes of memory left unfreed.\n", memAllocated);

	// Done
	return 0;
}

