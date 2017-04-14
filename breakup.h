/*
 * Breaks a big box into a series of smaller boxes. 
 * Assumes atom positions have been read in to global variables.
 */
void breakup(int bx, int by, int bz, LAPACK *atomMatrix, int *atomMatrixLabels, LAPACK ***atomMatrixOutput, int ***atomMatrixLabelsOutput) {
	int i,j,tmp[3],k,l,kx,ky,kz,tmp2[3];
	FTYPE fj,ftmp;

	// matrices
	int curAtom = 0;
	LAPACK *amo = NULL;
	bool *siteUsed = NULL; // numCells x numSites (false if no atom in that position yet)

	// results
	// atomMatrixOutput: bx*by*bz pointers to matrices, each of size numAtoms/(bx*by*bz) x 3 (real)
	// atomMatrixLabelsOutput: bx*by*bz pointers to arrays, each of size numAtoms/(bx*by*bz)

        if (atomMatrix == NULL || atomMatrixLabels == NULL)
                printErrorAndExit(16, "breakup: One or more input matrices not allocated!\n");

	if (atomMatrixOutput == NULL || (*atomMatrixOutput) == NULL || atomMatrixLabelsOutput == NULL || (*atomMatrixLabelsOutput) == NULL)
		printErrorAndExit(16, "breakup: One or more output matrices not allocated!\n");

	if (mod(cellsPerSupercell[0],bx) != 0 || bx < 1 || mod(cellsPerSupercell[1],by) != 0 || by < 1 ||
	    mod(cellsPerSupercell[2],bz) != 0 || bz < 1)
		printErrorAndExit(16, "breakup: bx=%i, by=%i, and/or bz=%i invalid (cellsPerSupercell: %i,%i,%i)!\n",
			bx, by, bz, cellsPerSupercell[0], cellsPerSupercell[1], cellsPerSupercell[2]);

	tmp2[0] = cellsPerSupercell[0]/bx;
	tmp2[1] = cellsPerSupercell[1]/by;
	tmp2[2] = cellsPerSupercell[2]/bz;

	/*
	 * Assign input atom positions to the closest corresponding "ideal"
	 * position, making sure there are no duplicates.
	 */

	// allocate siteUsed
	siteUsed = safe_malloc((UINT64)sizeof(bool)*(UINT64)numSites*(UINT64)numCells);

	// match atomMatrix atoms to atoms in 'ideal' xyz input
	for (curAtom = 0; curAtom < numAtoms; curAtom++) {
	    // find closest matching site
	    fj = rluMinDistance(lp_idx((atomMatrix),curAtom,0,LP_REAL), lp_idx((atomMatrix),curAtom,1,LP_REAL), lp_idx((atomMatrix),curAtom,2,LP_REAL),
			lp_idx(idealPosition,0,0,LP_REAL), lp_idx(idealPosition,0,1,LP_REAL), lp_idx(idealPosition,0,2,LP_REAL), LtoC);
	    j = 0;
	    for (i = 1; i < numSites; i++) {
            	ftmp = rluMinDistance(lp_idx((atomMatrix),curAtom,0,LP_REAL), lp_idx((atomMatrix),curAtom,1,LP_REAL), lp_idx((atomMatrix),curAtom,2,LP_REAL),
                        lp_idx(idealPosition,i,0,LP_REAL), lp_idx(idealPosition,i,1,LP_REAL), lp_idx(idealPosition,i,2,LP_REAL), LtoC);
		if (ftmp < fj) { fj = ftmp; j = i; }
	    }
            // find subcell indices
            for (i = 0; i < 3; i++)
                tmp[i] = mod(prip5(lp_idx((atomMatrix),curAtom,i,LP_REAL)-lp_idx(idealPosition,j,i,LP_REAL)),cellsPerSupercell[i]);

	    // sanity checks
	    if (fj > MAX_IDEAL_ACTUAL_DISTANCE)
		printWarning("Warning: Atom %i is %.6f > %.6f A away from the nearest ideal position in subcell (%i,%i,%i).\n", 
			curAtom, fj, MAX_IDEAL_ACTUAL_DISTANCE, tmp[0], tmp[1], tmp[2]);
	    if (strcmp(idealPositionLabel[(atomMatrixLabels)[curAtom]], idealPositionLabel[j]) != 0)
		printWarning("Warning: Atom %i label %s does not match ideal position label %s in subcell (%i,%i,%i).\n",
			curAtom, idealPositionLabel[(atomMatrixLabels)[curAtom]], idealPositionLabel[j], tmp[0], tmp[1], tmp[2]);

	    // mark spot as used
	    i = cellIndex(tmp[0],tmp[1],tmp[2]);
	    if (siteUsed[i*numSites+j] != false)
		printErrorAndExit(8, "Error: Atom %i is closest to atom %i in subcell (%i,%i,%i), but another already matches!\n",
			curAtom, j, tmp[0], tmp[1], tmp[2]);
	    siteUsed[i*numSites+j] = true;
	    // copy info to relevant sub-box
	    amo = (*atomMatrixOutput)[(tmp[0]/tmp2[0])*bz*by+(tmp[1]/tmp2[1])*bz+(tmp[2]/tmp2[2])];
	    i = (mod(tmp[0],tmp2[0])*tmp2[1]*tmp2[2] +
		mod(tmp[1],tmp2[1])*tmp2[2] +
		mod(tmp[2],tmp2[2]))*numSites+j;

	    lp_idx(amo,i,0,LP_REAL) = lp_idx((atomMatrix),curAtom,0,LP_REAL) - ((int)(tmp[0]/tmp2[0])) * tmp2[0];

	    if (lp_idx(amo,i,0,LP_REAL) < 0.0) lp_idx(amo,i,0,LP_REAL) += (FTYPE)tmp2[0];
	    else if (lp_idx(amo,i,0,LP_REAL) >= (FTYPE)tmp2[0]) lp_idx(amo,i,0,LP_REAL) -= (FTYPE)tmp2[0];

            lp_idx(amo,i,1,LP_REAL) = lp_idx((atomMatrix),curAtom,1,LP_REAL) - ((int)(tmp[1]/tmp2[1])) * tmp2[1];

            if (lp_idx(amo,i,1,LP_REAL) < 0.0) lp_idx(amo,i,1,LP_REAL) += (FTYPE)tmp2[1];
            else if (lp_idx(amo,i,1,LP_REAL) >= (FTYPE)tmp2[1]) lp_idx(amo,i,1,LP_REAL) -= (FTYPE)tmp2[1];

            lp_idx(amo,i,2,LP_REAL) = lp_idx((atomMatrix),curAtom,2,LP_REAL) - ((int)(tmp[2]/tmp2[2])) * tmp2[2];

            if (lp_idx(amo,i,2,LP_REAL) < 0.0) lp_idx(amo,i,2,LP_REAL) += (FTYPE)tmp2[2];
            else if (lp_idx(amo,i,2,LP_REAL) >= (FTYPE)tmp2[2]) lp_idx(amo,i,2,LP_REAL) -= (FTYPE)tmp2[2];

	    (*atomMatrixLabelsOutput)[(tmp[0]/tmp2[0])*bz*by+(tmp[1]/tmp2[1])*bz+(tmp[2]/tmp2[2])][i] =
		(atomMatrixLabels)[curAtom];
	}

	safe_free(siteUsed);
}

