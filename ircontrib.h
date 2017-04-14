/*
 * Read irrep contributions from a single file. Stored in a numCells x numIrreps array (complex)
 */
void readIRContribFile(int f, char *file) {
	char *curLine;
	int k,kx,ky,kz,i,j;
	FTYPE fi,fj,fk;
	int tmp[6];
	FTYPE ftmp[6];

        #pragma omp single
        {
	 openFile(f,file);
         // first line contains number of irreps, supercell multiplication factors, supercell edge lengths and angles,
         // and whether global offsets are removed
         curLine = getNextLine(f); 
	 if (!curLine)
		printErrorAndExit(6, "There is no first line in IRcontrib file %s!\n", file);
         i = sscanf(curLine, "%i %i %i %i %lf %lf %lf %lf %lf %lf %i %i", &(tmp[0]), &(tmp[1]), &(tmp[2]), &(tmp[3]),
                &(ftmp[0]), &(ftmp[1]), &(ftmp[2]), &(ftmp[3]), &(ftmp[4]), &(ftmp[5]), &(tmp[4]), &(tmp[5]));
         if (i != 12) 
		printErrorAndExit(6, "Error: Only read %i items from line 1 of %s: %s (expected 12)\n", i, file, curLine);
         if (tmp[4] != 0) {
		printStatus("A global cell shift was removed before constructing the irrep file %s.\n", file);
		irrepOffset = true;
	 } else { irrepOffset = false; }
         if (tmp[4] != 0 && tmp[5] != 0) { 
		printWarning("Warning: The aforementioned cell shift was applied on a per-site basis (and will likely give incorrect results).\n");
		irrepOffsetPerSite = true;
	 } else { irrepOffsetPerSite = false; }

	 if ( (numIrreps <= 0 && tmp[0] <= 0) || (tmp[0] != numIrreps && numIrreps > 0))
		printErrorAndExit(6, "Number of irreps = %i is invalid (%i)\n", tmp[0], numIrreps);
	 if ( (cellsPerSupercell[0] <= 0 && tmp[1] <= 0) || (tmp[1] != cellsPerSupercell[0] && cellsPerSupercell[0] > 0))
		printErrorAndExit(6, "Cells per supercell x = %i is invalid (%i)\n", tmp[1], cellsPerSupercell[0]);
         if ( (cellsPerSupercell[1] <= 0 && tmp[2] <= 0) || (tmp[2] != cellsPerSupercell[1] && cellsPerSupercell[1] > 0))
                printErrorAndExit(6, "Cells per supercell y = %i is invalid (%i)\n", tmp[2], cellsPerSupercell[1]);
         if ( (cellsPerSupercell[2] <= 0 && tmp[3] <= 0) || (tmp[3] != cellsPerSupercell[2] && cellsPerSupercell[2] > 0))
                printErrorAndExit(6, "Cells per supercell z = %i is invalid (%i)\n", tmp[3], cellsPerSupercell[2]);
	 if ( (supercellEdge[0] < ZERO_TOL && ftmp[0] < ZERO_TOL) || (fabs(supercellEdge[0]-ftmp[0]) > ZERO_TOL && supercellEdge[0] >= ZERO_TOL))
		printErrorAndExit(6, "Supercell edge x = %.8f is invalid (%.8f)\n", ftmp[0], supercellEdge[0]);
         if ( (supercellEdge[1] < ZERO_TOL && ftmp[1] < ZERO_TOL) || (fabs(supercellEdge[1]-ftmp[1]) > ZERO_TOL && supercellEdge[1] >= ZERO_TOL))
                printErrorAndExit(6, "Supercell edge y = %.8f is invalid (%.8f)\n", ftmp[1], supercellEdge[1]);
         if ( (supercellEdge[2] < ZERO_TOL && ftmp[2] < ZERO_TOL) || (fabs(supercellEdge[2]-ftmp[2]) > ZERO_TOL && supercellEdge[2] >= ZERO_TOL))
                printErrorAndExit(6, "Supercell edge z = %.8f is invalid (%.8f)\n", ftmp[2], supercellEdge[2]);
         if ( (cellAngles[0] < ZERO_TOL && ftmp[3] < ZERO_TOL) || (fabs(cellAngles[0]-ftmp[3]*PI/180.0) > ZERO_TOL && cellAngles[0] >= ZERO_TOL))
                printErrorAndExit(6, "Supercell angle alpha = %.8f is invalid (%.8f)\n", ftmp[3]*PI/180.0, cellAngles[0]);
         if ( (cellAngles[1] < ZERO_TOL && ftmp[4] < ZERO_TOL) || (fabs(cellAngles[1]-ftmp[4]*PI/180.0) > ZERO_TOL && cellAngles[1] >= ZERO_TOL))
                printErrorAndExit(6, "Supercell angle beta = %.8f is invalid (%.8f)\n", ftmp[4]*PI/180.0, cellAngles[1]);
         if ( (cellAngles[2] < ZERO_TOL && ftmp[5] < ZERO_TOL) || (fabs(cellAngles[2]-ftmp[5]*PI/180.0) > ZERO_TOL && cellAngles[2] >= ZERO_TOL))
                printErrorAndExit(6, "Supercell angle gamma = %.8f is invalid (%.8f)\n", ftmp[5]*PI/180.0, cellAngles[2]);

	 numIrreps = tmp[0];
	 for (i = 0; i < 3; i++) {
	    cellsPerSupercell[i] = tmp[i+1];
	    supercellEdge[i] = ftmp[i];
	    cellAngles[i] = ftmp[i+3]*PI/180.0;
	    cellEdge[i] = supercellEdge[i] / (FTYPE)cellsPerSupercell[i];
	 }
	 numCells = cellsPerSupercell[0]*cellsPerSupercell[1]*cellsPerSupercell[2];
         numUnique = (cellsPerSupercell[0]/2)+1 + (cellsPerSupercell[1]/2)*cellsPerSupercell[2] + (cellsPerSupercell[0]/2)*cellsPerSupercell[1]*cellsPerSupercell[2];
         if (numUnique > MAX_WAVEVECTORS || numCells > MAX_WAVEVECTORS)
                printErrorAndExit(8, "Too many wavevectors: numIrreps=%i, cellsPerSupercell=(%i,%i,%i), numUnique=%i, numCells=%i\n", numIrreps, cellsPerSupercell[0],
                        cellsPerSupercell[1], cellsPerSupercell[2], numUnique, numCells);

	 // Construct LtoC and CtoL matrices (from cellAngles, supercellEdge, cellEdge)
	 BuildLtoCandCtoL();

	 // Read in offset vectors (cartesian, x,y,z) after allocating memory
	 // 1 if global, numSites if done per-site
	 if (tmp[4] != 0) {

	    // memory allocation
	    if (tmp[5] != 0 && numSites < 1) printErrorAndExit(7, "Offsets removed, but numSites not yet defined! Read in XYZ file first.\n");
	    if (irrepOffsets != NULL && (lp_rows(irrepOffsets) != ((tmp[5] != 0) ? (numSites) : (1)) || lp_cols(irrepOffsets) != 3)) {
		printStatus("Reallocating irrepOffsets array from %ix%i to %ix%i\n", lp_rows(irrepOffsets), lp_cols(irrepOffsets), ((tmp[5] != 0) ? (numSites) : (1)), 3);
		lp_free(irrepOffsets);
		irrepOffsets = NULL;
	    }

	    if (irrepOffsets == NULL)
		irrepOffsets = lp_malloc(((tmp[5] != 0) ? (numSites) : (1)),3,false);

	    // reading of vectors
	    for (i = 0; i < ( (tmp[5] != 0) ? (numSites) : (1) ); i++) {
                curLine = getNextLine(f); if (!curLine) { printErrorAndExit(9, "Premature end of file at global offset %i!\n", i); }
                j = sscanf(curLine, "%lf %lf %lf", &fi, &fj, &fk);
                if (j != 3) printWarning("Warning: Only read %i items on line offset=%i: %s (expected 3)\n", j, i, curLine);
                lp_idx(irrepOffsets, i, 0, LP_REAL) = fi;
                lp_idx(irrepOffsets, i, 1, LP_REAL) = fj;
		lp_idx(irrepOffsets, i, 2, LP_REAL) = fk;
	    }
	 }

	 // (re)allocate memory for irrepCoeffs as needed
	 if (irrepCoeffs != NULL && (lp_rows(irrepCoeffs) != numCells || lp_cols(irrepCoeffs) != numIrreps)) {
		printStatus("Reallocating irrepCoeffs array from %ix%i to %ix%i\n", lp_rows(irrepCoeffs), lp_cols(irrepCoeffs), numCells, numIrreps);
		lp_free(irrepCoeffs);
		irrepCoeffs = NULL;
	 }
	 if (irrepCoeffs == NULL)
		irrepCoeffs = lp_malloc(numCells, numIrreps, true);

	 /*
	  * Read in the individual coefficients as a function of
	  * irrep and wavevector. Order is:
	  * irrep1 kx=0 ky=0 kz=0 (real+imag)
	  * irrep1 kx=0 ky=0 kz=1 (real+imag)
	  * ...
	  * with the irrep order the same as that in the irrep file
	  */

         // successive lines contain all IR contribs. Read directly into array
	 k = 0;
       	 for (i = 0; i < numIrreps; i++) {
           for (kx = 0; kx < cellsPerSupercell[0]; kx++) {
             for (ky = 0; ky < cellsPerSupercell[1]; ky++) {
               for (kz = 0; kz < cellsPerSupercell[2]; kz++) {
		curLine = getNextLine(f); if (!curLine) { printErrorAndExit(10, "Premature end of file at irrep=%i k=(%i,%i,%i)!\n", i, kx, ky, kz); }
		j = sscanf(curLine, "%lf %lf", &fi, &fj);
		if (j != 2) printWarning("Warning: Only read %i items on line irrep=%i k=(%i,%i,%i): %s (expected 2)\n", j, i, kx, ky, kz, curLine);
		lp_idx(irrepCoeffs, cellIndex(kx,ky,kz), i, LP_REAL) = fi;
		lp_idx(irrepCoeffs, cellIndex(kx,ky,kz), i, LP_IMAG) = fj;
		k++;
	       }
	     }
	   }
	 }
	 closeFile(f);
         printStatus("Total of %i IRcontribs read in (not printing).\n", k);
	} // end single
}

void writeIRContribFile(char *opfName, LAPACK *irMatrix, bool irOffset, bool irOffsetPerSite, LAPACK *irOffsetMatrix) {
        FILE *opf;
	int i,kx,ky,kz;
        opf = fopen(opfName, "wb");
        if (!opf) printErrorAndExit(7, "Error: Cannot open output file %s!\n", opfName);
        printStatus("Outputting irrep coefficients to irrep file %s...", opfName);
        // first line contains number of irreps, supercell multiplication factors, supercell edge lengths and angles
        // and whether global offsets (overall, or per site) are removed
	fprintf(opf, "%i %i %i %i %.8f %.8f %.8f %.5f %.5f %.5f %i %i\n", numIrreps, cellsPerSupercell[0], cellsPerSupercell[1], cellsPerSupercell[2],
		supercellEdge[0], supercellEdge[1], supercellEdge[2], cellAngles[0]*180.0/PI, cellAngles[1]*180.0/PI, cellAngles[2]*180.0/PI, irOffset, irOffsetPerSite);

	// if global offsets are used, next 1 (or numSites) lines contain offset vectors
	if (irOffset) {
	    for (i = 0; i < lp_rows(irOffsetMatrix); i++)
	        fprintf(opf, "%.9g %.9g %.9g\n", lp_idx(irOffsetMatrix,i,0,LP_REAL), lp_idx(irOffsetMatrix,i,1,LP_REAL), lp_idx(irOffsetMatrix,i,2,LP_REAL));
	}

	// successive lines contain real and imaginary components of each contribution
        for (i = 0; i < numIrreps; i++) {
           for (kx = 0; kx < cellsPerSupercell[0]; kx++) {
             for (ky = 0; ky < cellsPerSupercell[1]; ky++) {
               for (kz = 0; kz < cellsPerSupercell[2]; kz++) {
		fprintf(opf, "%.9g %.9g\n", lp_idx(irMatrix,cellIndex(kx,ky,kz),i,LP_REAL), lp_idx(irMatrix,cellIndex(kx,ky,kz),i,LP_IMAG));
               }
             }
           }
        }

        fclose(opf);
        printStatus("Done.\n");
}

