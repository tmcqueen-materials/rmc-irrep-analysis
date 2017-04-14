void readEyeFileHeader(int f) {
	char *curLine = NULL;
	char matchStr[16] = { 0 };
	int i = 0, j = 0, k = 0, l = 0;
	FTYPE fi;
	LAPACK *H = NULL;

	// First line is number of particles
	curLine = getNextLine(f);
	if (!curLine) printErrorAndExit(7, "Error: Cannot read line 1 from eye file!\n");
	i = sscanf(curLine, "Number of particles = %i", &j);
	if (i < 1 || j < 1) printErrorAndExit(7, "Error: Missing or invalid number of atoms (%i)!\n", j);
	numAtoms = j;
	curLine = getNextLine(f); // "A = 1.0 Angstrom (basic length-scale)", ignore [stronger file checking could be done]
        if (!curLine) printErrorAndExit(7, "Error: Cannot read line 2 from eye file!\n");

	/*
	 * The matrix H is stored in column-major format, and H is the matrix such that C = H*L, where L is the 
	 * atom positions as given in the eyefile in supercell relative lattice units, and C is a cartesian coordinate
	 * system in angstrom units.
	 *
	 * To be consistent, we actually extract the supercell edge lengths and angles from this matrix, and then
	 * setup our own LtoC matrix that is always in upper triangular form.
	 */
	H = lp_malloc(3,3,false);
	for (i = 1; i < 4; i++) {
	    for (j = 1; j < 4; j++) {
		curLine = getNextLine(f);
		if (!curLine) printErrorAndExit(7, "Error: Cannot read line %i from eye file!\n", 3+3*(i-1)+(j-1));
		sprintf(matchStr, "H0(%i,%i) = %%lf A",i,j);
		k = sscanf(curLine, matchStr, &fi);
		if (k < 1) printErrorAndExit(7, "Error, missing or invalid H(%i,%i) value: %.6f\n", i, j, fi);
		lp_idx(H,j-1,i-1,LP_REAL) = fi;
	    }
	}

	for (i = 0; i < 3; i++)
	    supercellEdge[i] = magCol(H,i);
	cellAngles[0] = angCol(H,1,2);
	cellAngles[1] = angCol(H,0,2);
	cellAngles[2] = angCol(H,0,1);

	lp_free(H);

        // query user for number of subcells along each dimension, if needed
        if (cellsPerSupercell[0] == 0 && cellsPerSupercell[1] == 0 && cellsPerSupercell[2] == 0) {
            printStatus("Supercell is %.4f x %.4f x %.4f A, with %i atoms.\n", supercellEdge[0], supercellEdge[1], supercellEdge[2], numAtoms);
            printStatus("Enter number of subcells along x, y, and z, separated by spaces:");
            i = fscanf(stdin, "%i %i %i", &j, &k, &l);
            if (i != 3 || j < 1 || k < 1 || l < 1) printErrorAndExit(8, "One or more number of subcells per supercell missing or invalid!\n");
            cellsPerSupercell[0] = j;
            cellsPerSupercell[1] = k;
            cellsPerSupercell[2] = l;
            numCells = j*k*l;
            for (i = 0; i < 3; i++) cellEdge[i] = supercellEdge[i] / (FTYPE)cellsPerSupercell[i];
        }

	// Fill in conversion matrices
	BuildLtoCandCtoL();

	// make sure eta matrix is 0
        curLine = getNextLine(f);
        if (!curLine) printErrorAndExit(7, "Error: Cannot read line 12 from eye file!\n");
        i = sscanf(curLine, "eta(1,1) = %lf A", &fi);
        if (i < 1 || fabs(fi) > ZERO_TOL) printErrorAndExit(7, "Error on line 12, missing or non-zero eta value: %.6f\n", fi);
        curLine = getNextLine(f);
        if (!curLine) printErrorAndExit(7, "Error: Cannot read line 13 from eye file!\n");
        i = sscanf(curLine, "eta(1,2) = %lf A", &fi);
        if (i < 1 || fabs(fi) > ZERO_TOL) printErrorAndExit(7, "Error on line 13, missing or non-zero eta value: %.6f\n", fi);
        curLine = getNextLine(f);
        if (!curLine) printErrorAndExit(7, "Error: Cannot read line 14 from eye file!\n");
        i = sscanf(curLine, "eta(1,3) = %lf A", &fi);
        if (i < 1 || fabs(fi) > ZERO_TOL) printErrorAndExit(7, "Error on line 14, missing or non-zero eta value: %.6f\n", fi);
        curLine = getNextLine(f);
        if (!curLine) printErrorAndExit(7, "Error: Cannot read line 15 from eye file!\n");
        i = sscanf(curLine, "eta(2,2) = %lf A", &fi);
        if (i < 1 || fabs(fi) > ZERO_TOL) printErrorAndExit(7, "Error on line 15, missing or non-zero eta value: %.6f\n", fi);
        curLine = getNextLine(f);
        if (!curLine) printErrorAndExit(7, "Error: Cannot read line 16 from eye file!\n");
        i = sscanf(curLine, "eta(2,3) = %lf A", &fi);
        if (i < 1 || fabs(fi) > ZERO_TOL) printErrorAndExit(7, "Error on line 16, missing or non-zero eta value: %.6f\n", fi);
        curLine = getNextLine(f);
        if (!curLine) printErrorAndExit(7, "Error: Cannot read line 17 from eye file!\n");
        i = sscanf(curLine, "eta(3,3) = %lf A", &fi);
        if (i < 1 || fabs(fi) > ZERO_TOL) printErrorAndExit(7, "Error on line 17, missing or non-zero eta value: %.6f\n", fi);

}

void readEyeFile(int f, char *file) {
	char *curLine;
	int i,j,k,l,curAtom,tmp[3];
	FTYPE ftmp[4];
        char atomLabel[MAX_IDEAL_SITES_LABEL_LEN] = { 0 };

       #pragma omp single
       {
	openFile(f, file);
	readEyeFileHeader(f); // fills in supercellEdge, numAtoms

        // (re)allocate memory for atomMatrix and atomMatrixLabels as needed
        if (atomMatrix != NULL && (lp_rows(atomMatrix) != numAtoms || lp_cols(atomMatrix) != 3)) {
                printStatus("Reallocating atomMatrix array from %ix%i to %ix%i\n", lp_rows(atomMatrix), lp_cols(atomMatrix), numAtoms, 3);
                lp_free(atomMatrix);
                atomMatrix = NULL;
		safe_free(atomMatrixLabels);
		atomMatrixLabels = NULL;
        }
        if (atomMatrix == NULL)
                atomMatrix = lp_malloc(numAtoms,3,false);
	if (atomMatrixLabels == NULL)
		atomMatrixLabels = safe_malloc((UINT64)numAtoms*(UINT64)sizeof(int));

        /*
         * Read in atom positions and labels.
         */

        for (curAtom = 0; curAtom < numAtoms; curAtom++) {
            curLine = getNextLine(f);
            i = sscanf(curLine, " %lf%n", &(ftmp[3]), &j);
            if (i < 1 || fabs(ftmp[3]) > ZERO_TOL) printErrorAndExit(3, "Error reading input atom float from eye file %s, atom %i!\n", file, curAtom);
            curLine += j;
            i = sscanf_str(curLine, atomLabel, MAX_IDEAL_SITES_LABEL_LEN);
            if (i < 1 || strlen(atomLabel) < 1) printErrorAndExit(3, "Error reading input atom label from eye file %s, atom %i!\n", file, curAtom);
            i = sscanf(curLine, " %*s %lf %lf %lf %i %i %i", &(ftmp[0]), &(ftmp[1]), &(ftmp[2]), &(tmp[0]), &(tmp[1]), &(tmp[2]));
            if (i < 6 || tmp[0] != 0 || tmp[1] != 0 || tmp[2] != 0) 
		printErrorAndExit(3, "Error reading input atom coordinates from eye file %s, atom %i!\n", file, curAtom);

	    // find matching atom label
	    (atomMatrixLabels)[curAtom] = -1;
	    for (j = 0; j < numSites; j++)
		if (strcmp(atomLabel, idealPositionLabel[j]) == 0) (atomMatrixLabels)[curAtom] = j;
	    if ((atomMatrixLabels)[curAtom] == -1)
		printErrorAndExit(8, "readEyeFile: Could not find atom label %s from supercell in ideal atom labels!\n", atomLabel);

            // fill in atom matrix
	    for (j = 0; j < 3; j++)
		lp_idx((atomMatrix),curAtom,j,LP_REAL) = (FTYPE)cellsPerSupercell[j]*mod1(ftmp[j]);

        }

	closeFile(f);
       } // end omp single
}

// BEGIN MINI EYE FILE OUTPUT FUNCTION

void writeEyeFile(char *opfName, LAPACK *atomMatrix, int *atomMatrixLabels) {
	FILE *opf;
	int i,j,k,l,m = 0;
        opf = fopen(opfName, "wb");
        if (!opf) printErrorAndExit(7, "Error: Cannot open output file %s!\n", opfName);
        printStatus("Outputting atom positions to eye file %s...", opfName);
        // header first
        fprintf(opf, "Number of particles = %i\n", numAtoms);
        fprintf(opf, "A = 1.0 Angstrom (basic length-scale)\n");
	for (i = 0; i < 3; i++) {
	    for (j = 0; j < 3; j++) {
		fprintf(opf, "H0(%i,%i) = %.6f A\n", i+1, j+1, lp_idx(ssLtoC,j,i,LP_REAL));
	    }
	}
        fprintf(opf, "eta(1,1) = 0.000000\n");
        fprintf(opf, "eta(1,2) = 0.000000\n");
        fprintf(opf, "eta(1,3) = 0.000000\n");
        fprintf(opf, "eta(2,2) = 0.000000\n");
        fprintf(opf, "eta(2,3) = 0.000000\n");
        fprintf(opf, "eta(3,3) = 0.000000\n");
        for (l = 0; l < numSites; l++) {
          for (i = 0; i < cellsPerSupercell[0]; i++) {
            for (j = 0; j < cellsPerSupercell[1]; j++) {
              for (k = 0; k < cellsPerSupercell[2]; k++) {
                m = i*cellsPerSupercell[1]*cellsPerSupercell[2]*numSites+j*cellsPerSupercell[2]*numSites+k*numSites+l;
                fprintf(opf,"  0.000 %s  %.7f %.7f %.7f 0 0 0\n", idealPositionLabel[atomMatrixLabels[m]],
                        lp_idx(atomMatrix,m,0,LP_REAL)/(FTYPE)cellsPerSupercell[0],
                        lp_idx(atomMatrix,m,1,LP_REAL)/(FTYPE)cellsPerSupercell[1],
                        lp_idx(atomMatrix,m,2,LP_REAL)/(FTYPE)cellsPerSupercell[2]);
              }
            }
          }
        }
        fclose(opf);
        printStatus("Done.\n");
}

// END MINI EYE FILE OUTPUT FUNCTION

