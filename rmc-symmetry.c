#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdint.h>
#include <math.h>
#include "defs.h"
#include "globals.h"
#define USAGE_MESSAGE "Script Use: rmc-symmetry [xyzinput].cfg [irrep].cfg [orientations].cfg [eyeFilesToAnalyze].txt [perIrrepResults].txt\n"
#include "io.h"
#include "memlapack.h"
#include "supercell.h"
#include "xyzinput.h"
#include "irrepinput.h"
#include "ircontrib.h"
#include "eyefile.h"
#include "orientation.h"
#include "irreprebuild.h"
#include "irrepdecompose.h"

int main(int argc, char *argv[]) { // argv[1] = xyz input, argv[2] = irreps, argv[3] = list of orientations, argv[4] = file listing eyeFilesToAnalyze
	int j,k,m,n;
	FILE *opf = NULL;
	COO mir;
	int kx,ky,kz;
	FTYPE fi,fj,fk,fl;
	FTYPE runAvg[MAX_ORIENTATIONS][2] = { { 0.0 } }; // mean and sum of squares of differences (for variance)
	FTYPE curRun[MAX_ORIENTATIONS] = { 0.0 };
	int numAvg = 0;
	char *curLine;
	LAPACK *irrepCoeffsOriented = NULL; // numCells x numIrreps (complex)
	LAPACK *irrepOffsetsOriented = NULL; // 1 x 3 (real)
	LAPACK *irrepCoeffsUnoriented = NULL; // numCells x numIrreps (complex)
	LAPACK *irrepOffsetsUnoriented = NULL; // 1 x 3 (real)
	LAPACK *atomsOriented = NULL; // numAtoms x 3 (real)

	// check number of input parameters
	if (argc != 6) printUsageAndExit(1);

	// Disable output buffering of status and error messages (helps debugging)
	setvbuf(stdout, NULL, _IONBF, 0);
	setvbuf(stderr, NULL, _IONBF, 0);

	/*
	 * Read xyz file
	 */
	readXYZFile(TMP_FILE, argv[1]);

        /*
         * Open list of files to analyze
         */
        openFile(FTA_FILE, argv[4]);

	/*
	 * Open perIrrep output file
	 */
	opf = fopen(argv[5], "wb");
	if (!opf) printErrorAndExit(2, "Could not open %s for outputting per-irrep info!\n", argv[5]);

        /*
         * There must be at least one file to open and analyze
         */
        curLine = getNextLine(FTA_FILE);
        if (curLine == NULL || strlen(curLine) < 1) printErrorAndExit(2, "No files to open in eyeFilesToAverage.txt!\n");

        /*
         * Read in atom positions from first file (to fill in necessary globals for next file reading steps)
         */
        readEyeFile(TMP_FILE, curLine);

        /*
         * Read irreps file
         */
        readIrrepFile(TMP_FILE, argv[2]);

        /*
         * Read orientation file
         */
        readOrientationFile(TMP_FILE, argv[3]);

	/*
 	 * Allocate memory to hold the [un]oriented irrep coefficients and oriented atom positions
	 */
	irrepCoeffsOriented = lp_malloc(numCells,numIrreps,true);
	irrepOffsetsOriented = lp_malloc(1,3,false); 
	irrepCoeffsUnoriented = lp_malloc(numCells,numIrreps,true);
	irrepOffsetsUnoriented = lp_malloc(1,3,false);
	atomsOriented = lp_malloc(numAtoms, 3, false);

	/*
	 * Print running table of results
	 */
	printStatus("  Run  ");
	for (j = 0; j < numOrientations; j++)
	    printStatus(" %10i", j);
	printStatus("\n");
	printStatus("-------");
	for (j = 0; j < numOrientations; j++)
	    printStatus(" ----------");
	printStatus("\n");

	while (curLine != NULL && strlen(curLine) > 0) {
	    numAvg++;

	    /*
	     * Read in atom positions
	     */
	    readEyeFile(TMP_FILE, curLine);

	    /*
	     * Compute baseline (unoriented) irrep map
	     */
	    irrepDecompose(atomMatrix, atomMatrixLabels, &irrepCoeffsUnoriented, true, false, &irrepOffsetsUnoriented);

	    printStatus("%7i", numAvg);

            /*
             * Sequentially apply every possible orientation/order and compute correlation with input
             */
            for (j = 0; j < numOrientations; j++) { // orientations
                // compute oriented irrep coefficients
                applyOrientation(j, atomMatrix, &atomsOriented);
                irrepDecompose(atomsOriented, atomMatrixLabels, &irrepCoeffsOriented, true, false, &irrepOffsetsOriented);

		fprintf(opf, "%7i %7i", numAvg, j);

		// Compute the least squares metric of how close to "identical" the two are. We do that by computing the following for each unique irrep set,
		// on a per-wavevector basis:
		// sum over all degenerate irreps of individual |orig|^2 (=(a^2+b^2) where a+bi is the coefficient) 
		// sum over all degenerate irreps of individual |new|^2 (=(a^2+b^2) where a+bi is the coefficient)
		// sum of (|orig|^2-|new|^2)^2 over all wavevectors
		// This produces a least-squares type number indicating quality of the match (0 = perfect matching) on that irrep set
		//
		// These can then be summed over all irrep sets to give a number indicating quality of the match for that symmetry operation
		fl = 0.0; // sum over all wavevectors and all irrep sets
                for (n = 0; n < numIrrepByDegen; n++) {
		  fk = 0.0; // sum over all wavevectors
                  for (kx = 0; kx < cellsPerSupercell[0]; kx++) {
                    for (ky = 0; ky < cellsPerSupercell[1]; ky++) {
                      for (kz = 0; kz < cellsPerSupercell[2]; kz++) {
			fi = 0.0; // orig
			fj = 0.0; // new
                        for (k = 0; k < irrepByDegen[n][0]; k++) {
                          m = irrepByDegen[n][k+1];
			  fi += lp_idx(irrepCoeffsUnoriented,FtofullF(kx,ky,kz),m,LP_REAL)*lp_idx(irrepCoeffsUnoriented,FtofullF(kx,ky,kz),m,LP_REAL)+lp_idx(irrepCoeffsUnoriented,FtofullF(kx,ky,kz),m,LP_IMAG)*lp_idx(irrepCoeffsUnoriented,FtofullF(kx,ky,kz),m,LP_IMAG);
 			  fj += lp_idx(irrepCoeffsOriented,FtofullF(kx,ky,kz),m,LP_REAL)*lp_idx(irrepCoeffsOriented,FtofullF(kx,ky,kz),m,LP_REAL)+lp_idx(irrepCoeffsOriented,FtofullF(kx,ky,kz),m,LP_IMAG)*lp_idx(irrepCoeffsOriented,FtofullF(kx,ky,kz),m,LP_IMAG);
			}
			fk += (fi-fj)*(fi-fj);
                      }
                    }
                  }
		  fprintf(opf, " %10.5e", fk);
		  fl += fk;
                }
		curRun[j] = fl;
		printStatus(" %10.5e", fl);

		fprintf(opf, "\n");
            }
	    printStatus("\n");

	    /*
	     * Update averages
	     */
	    for (j = 0; j < numOrientations; j++) {
		fl = curRun[j]-runAvg[j][0];
		runAvg[j][0] = runAvg[j][0] + fl/(FTYPE)numAvg;
		runAvg[j][1] = runAvg[j][1] + fl*(curRun[j]-runAvg[j][0]);
	    }

            /*
	     * Get name of next file to analyze
	     */
	    curLine = getNextLine(FTA_FILE);
	}

	/*
	 * Print out averages
	 */
        printStatus("OVERALL");
	for (j = 0; j < numOrientations; j++)
	    printStatus(" %10.5e", runAvg[j][0]);
	printStatus("\n");
	printStatus(" STDDEV");
	for (j = 0; j < numOrientations; j++)
	    printStatus(" %10.5e", runAvg[j][1]/(FTYPE)numAvg);
	printStatus("\n\n");

	/*
	 * Close per-irrep file
	 */
	fclose(opf); opf = NULL;

	/*
	 * Close list of files to analyze
	 */
	closeFile(FTA_FILE);

	// Free memory
	lp_free(irrepCoeffsOriented);
	lp_free(irrepOffsetsOriented);
	lp_free(irrepCoeffsUnoriented);
	lp_free(irrepOffsetsUnoriented);
        lp_free(atomsOriented);
	freeGlobals();

	printStatus("Memory Used (excluding stack and overhead): %i MB\n", maxMemAllocated/1048576);

        if (memAllocated > (UINT64)0)
                printStatus("%i bytes of memory left unfreed.\n", memAllocated);

	// Done
	return 0;
}

