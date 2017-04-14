#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdint.h>
#include <math.h>
#include "defs.h"
#include "globals.h"
#define USAGE_MESSAGE "Script Use: rmc-breakup [xyzinput].cfg [eyeinput].eye\n"
#include "io.h"
#include "memlapack.h"
#include "supercell.h"
#include "xyzinput.h"
#include "irrepinput.h"
#include "ircontrib.h"
#include "eyefile.h"
#include "breakup.h"

int main(int argc, char *argv[]) { // argv[1] = xyzinput, argv[2] = eyeinput
	char *opfName = NULL;
	char toapp[7] = { 0 };
	int i,j,k,l;

	// results
	LAPACK **atomMatrixOutput = NULL;
	int **atomMatrixLabelsOutput = NULL;
	int bx = 0, by = 0, bz = 0;

	// check number of input parameters
	if (argc != 3) printUsageAndExit(1);

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
	readEyeFile(FTA_FILE, argv[2]);

        // query user for number of subboxes along each dimension
        printStatus("Supercell is %i x %i x %i sub-cells, with %i atoms.\n", cellsPerSupercell[0], cellsPerSupercell[1], cellsPerSupercell[2], numAtoms);
        printStatus("Enter number of sub-boxes desired along x, y, and z, separated by spaces:");
        i = fscanf(stdin, "%i %i %i", &j, &k, &l);
        if (i != 3 || j < 1 || k < 1 || l < 1) printErrorAndExit(8, "One or more number of sub-boxes per supercell missing or invalid!\n");
        bx = j;
        by = k;
        bz = l;

	// Allocate result matrices
	printStatus("Allocating result matrices...\n");
	atomMatrixOutput = safe_malloc((UINT64)bx*(UINT64)by*(UINT64)bz*(UINT64)sizeof(LAPACK *));
	atomMatrixLabelsOutput = safe_malloc((UINT64)bx*(UINT64)by*(UINT64)bz*(UINT64)sizeof(int *));
	for (i = 0; i < bx*by*bz; i++) {
		atomMatrixOutput[i] = lp_malloc(numAtoms/bx/by/bz,3,false);
		atomMatrixLabelsOutput[i] = safe_malloc((UINT64)(numAtoms/bx/by/bz)*(UINT64)sizeof(int));
	}

        printStatus("Generating sub-boxes...\n");
	breakup(bx,by,bz,atomMatrix,atomMatrixLabels,&atomMatrixOutput,&atomMatrixLabelsOutput);

	/*
	 * Adjust global variables to subbox size
	 */
	numAtoms = numAtoms/bx/by/bz;
	cellsPerSupercell[0] /= bx;
        cellsPerSupercell[1] /= by;
        cellsPerSupercell[2] /= bz;
        supercellEdge[0] /= (FTYPE)bx;
        supercellEdge[1] /= (FTYPE)by;
        supercellEdge[2] /= (FTYPE)bz;
		
	/*
	 * Write result to files
	 */
	printStatus("Writing Results to file...\n");
	for (i = 0; i < bx*by*bz; i++) {
	    snprintf(toapp,7,"_%i",i);
	    opfName = safe_malloc((UINT64)strlen(argv[2])+(UINT64)7);
	    strncpy(opfName, argv[2], strlen(argv[2])+7);
	    strncat(opfName, toapp, 7);
            writeEyeFile(opfName, atomMatrixOutput[i], atomMatrixLabelsOutput[i]);
	    safe_free(opfName); opfName = NULL;
	}

	printStatus("Memory Used (excluding stack and overhead): %i MB\n", maxMemAllocated/1048576);
	
	// Free memory 
        for (i = 0; i < bx*by*bz; i++) {
		lp_free(atomMatrixOutput[i]);
                safe_free(atomMatrixLabelsOutput[i]);
        }
        safe_free(atomMatrixOutput);
        safe_free(atomMatrixLabelsOutput);
	freeGlobals();

	if (memAllocated > (UINT64)0)
		printStatus("%i bytes of memory left unfreed.\n", memAllocated);

	// Done
	return 0;
}

