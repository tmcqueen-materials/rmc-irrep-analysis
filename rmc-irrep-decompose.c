#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdint.h>
#include <math.h>
#include "defs.h"
#include "globals.h"
#define USAGE_MESSAGE "Script Use: rmc-irrep-decompose [xyzinput].cfg [irrepinput].cfg [eyeinput].eye\n"
#include "io.h"
#include "memlapack.h"
#include "supercell.h"
#include "xyzinput.h"
#include "irrepinput.h"
#include "ircontrib.h"
#include "eyefile.h"
#include "irrepdecompose.h"

int main(int argc, char *argv[]) { // argv[1] = xyzinput, argv[2] = irrepinput, argv[3] = eyeinput
	char *opfName = NULL;

	// results
	LAPACK *irMatrix = NULL; // numCells x numIrreps (complex)

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
	 * Read in atom positions.
	 */
	readEyeFile(FTA_FILE, argv[3]);

        /*
         * Read in irreducible representations.
         */
        readIrrepFile(FTA_FILE,argv[2]);

	// Allocate result matrix
	printStatus("Allocating result matrix...\n");
	irMatrix = lp_malloc(numCells,numIrreps,true);

        printStatus("Generating irrep contributions...\n");
	irrepDecompose(atomMatrix,atomMatrixLabels,&irMatrix,false,false,NULL);
	
	/*
	 * Write result to file
	 */
	opfName = safe_malloc((UINT64)strlen(argv[3])+(UINT64)7);
	strncpy(opfName, argv[3], strlen(argv[3])+7);
	strncat(opfName, "_irrep", 7);
	writeIRContribFile(opfName, irMatrix, false, false, NULL);
	safe_free(opfName); opfName = NULL;

	printStatus("Memory Used (excluding stack and overhead): %i MB\n", maxMemAllocated/1048576);
	
	// Free memory 
	lp_free(irMatrix);
	freeGlobals();

	if (memAllocated > (UINT64)0)
		printStatus("%i bytes of memory left unfreed.\n", memAllocated);

	// Done
	return 0;
}

