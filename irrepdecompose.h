/*
 * Convert a set of supercell atom positions into irrep form.
 * Assumes xyz atom positions and irrep atom info has been read in to global variables.
 */
void irrepDecompose(LAPACK *atomMatrix, int *atomMatrixLabels, LAPACK **irCoeffs, bool irOffset, bool irOffsetPerSite, LAPACK **irOffsetMatrix) {
	int i,j,tmp[3],k,l,kx,ky,kz;
	FTYPE fi,fj,fk,ftmp;

	// matrices
	LAPACK *AtomstoIR = NULL; // numIrreps x numCells matrix (real) [= irrepTmp * deltaMatrix^T]
	int curAtom = 0;
	bool *deltaMatrixUsed = NULL; // numCells x numSites (false if no atom in that position yet)
	LAPACK *deltaMatrix = NULL; // numCells x 3*numSites matrix (real)
	LAPACK *tmpDM = NULL; // 3 x 1 (for of conversion of delta vector to cartesian, angstrom coordinates)
	LAPACK *irrepTmp = NULL; // numIrreps x 3*numSites matrix (real) [= LtoCnorm*irrepTmp]
        LAPACK *irOM = NULL; // numSitesx3

	// compute matrices
	LAPACK *phaseMatrix = NULL; // numCells x 1 (=e^ikx for all x for a single k) matrix (complex, or numCells x 2 if treated as real for BLAS/LAPACK)
	LAPACK *irMatrixSmall = NULL; // numIrreps x 1 matrix (complex, or numIrreps x 2 if treated as real for BLAS/LAPACK) [all coefficients for a single subcell (k value)]

	// results
	// irCoeffs: numCells x numIrreps (complex)
	// irOffsetMatrix: NULL, or 1x3 (irOffset defined), or numSitesx3 (irOffset and irOffsetPerSite defined)

	/*
	 * Assign input atom positions to the closest corresponding "ideal"
	 * position, making sure there are no duplicates. We directly store the 
	 * "delta" vectors from ideal for easy conversion to irreps.
	 */

	// allocate deltaMatrix and deltaMatrixUsed
	deltaMatrixUsed = safe_malloc((UINT64)sizeof(bool)*(UINT64)numSites*(UINT64)numCells);
	deltaMatrix = lp_malloc(numCells,3*numSites,false);
	tmpDM = lp_malloc(3,1,false);

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

	    // fill in appropriate location
	    i = cellIndex(tmp[0],tmp[1],tmp[2]);
	    if (deltaMatrixUsed[i*numSites+j] != false)
		printErrorAndExit(8, "Error: Atom %i is closest to atom %i in subcell (%i,%i,%i), but another already matches!\n",
			curAtom, j, tmp[0], tmp[1], tmp[2]);
	    deltaMatrixUsed[i*numSites+j] = true;
	    for (k = 0; k < 3; k++)
		lp_idx(tmpDM,k,0,LP_REAL) = modp5(lp_idx((atomMatrix),curAtom,k,LP_REAL)-lp_idx(idealPosition,j,k,LP_REAL));

	    cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1.0,lp_ptr(LtoC),lp_cols(LtoC),lp_ptr(tmpDM),1,0.0,&(lp_idx(deltaMatrix,i,3*j,LP_REAL)),1);
	}

	/*
	 * Remove global cell shifts, if desired.
         * Note that we *always* calculate per-Site cell shifts, in order to provide helpful
         * warnings if something looks horribly wrong.
	 */
	if (irOffset) {
	    // Allocate result matrix if needed
            if (!(*irOffsetMatrix))
                *irOffsetMatrix = lp_malloc(((irOffsetPerSite) ? (numSites) : (1)),3,false);
            else if (lp_rows((*irOffsetMatrix)) != ((irOffsetPerSite) ? (numSites) : (1)) || lp_cols((*irOffsetMatrix)) != 3)
                printErrorAndExit(12, "irrepDecompose: irOffsetMatrix already allocated, but wrong size!\n");
        }
        if (irOffset && irOffsetPerSite) irOM = *irOffsetMatrix;
        else irOM = lp_malloc(numSites,3,false);

	// In the deltaMatrix, the global offsets are simply the sum of the delta vectors over all atoms (over all atoms, and per site)
	fi = 0.0; fj = 0.0; fk = 0.0;
	for (i = 0; i < numSites; i++) {
	  for (j = 0; j < numCells; j++) {
	    fi += lp_idx(deltaMatrix,j,3*i+0,LP_REAL);
	    fj += lp_idx(deltaMatrix,j,3*i+1,LP_REAL);
	    fk += lp_idx(deltaMatrix,j,3*i+2,LP_REAL);
	  }
	  lp_idx(irOM,i,0,LP_REAL) = fi / (FTYPE)numCells;
	  lp_idx(irOM,i,1,LP_REAL) = fj / (FTYPE)numCells;
	  lp_idx(irOM,i,2,LP_REAL) = fk / (FTYPE)numCells;
	  fi = 0.0; fj = 0.0; fk = 0.0;
	}

        // single global cell shift if desired
	if (irOffset && !irOffsetPerSite) {
          fi = 0.0; fj = 0.0; fk = 0.0;
          for (i = 0; i < numSites; i++) {
            fi += lp_idx(irOM,i,0,LP_REAL);
            fj += lp_idx(irOM,i,1,LP_REAL);
            fk += lp_idx(irOM,i,2,LP_REAL);
          }
          lp_idx(*irOffsetMatrix,0,0,LP_REAL) = fi / (FTYPE)numSites;
          lp_idx(*irOffsetMatrix,0,1,LP_REAL) = fj / (FTYPE)numSites;
          lp_idx(*irOffsetMatrix,0,2,LP_REAL) = fk / (FTYPE)numSites;
	}

        // Make sure per-site shifts are "small" to check for errors in input xyz specification
        for (i = 0; i < numSites; i++) {
          fi = lp_idx(irOM,i,0,LP_REAL)*lp_idx(irOM,i,0,LP_REAL)+lp_idx(irOM,i,1,LP_REAL)*lp_idx(irOM,i,1,LP_REAL)+lp_idx(irOM,i,2,LP_REAL)*lp_idx(irOM,i,2,LP_REAL);
          if (sqrt(fi) > MAX_SYSTEMATIC_SHIFT)
            printWarning("Warning: Large per-site shift on site %6i: [ %10.5E %10.5E %10.5E ]\n", i, lp_idx(irOM,i,0,LP_REAL), lp_idx(irOM,i,1,LP_REAL), lp_idx(irOM,i,2,LP_REAL));
        }

        // Apply shift if desired
        if (irOffset) {
            for (i = 0; i < numSites; i++) {
              for (j = 0; j < numCells; j++) {
                lp_idx(deltaMatrix,j,3*i+0,LP_REAL) -= lp_idx(((irOffsetPerSite)?(irOM):(*irOffsetMatrix)),((irOffsetPerSite)?(i):(0)),0,LP_REAL);
                lp_idx(deltaMatrix,j,3*i+1,LP_REAL) -= lp_idx(((irOffsetPerSite)?(irOM):(*irOffsetMatrix)),((irOffsetPerSite)?(i):(0)),1,LP_REAL);
                lp_idx(deltaMatrix,j,3*i+2,LP_REAL) -= lp_idx(((irOffsetPerSite)?(irOM):(*irOffsetMatrix)),((irOffsetPerSite)?(i):(0)),2,LP_REAL);
              }
	    }
	}

        if (!(irOffset && irOffsetPerSite)) lp_free(irOM);

	/*
	 * Convert irrep vectors to cartesian, normalized coordinates
	 */
	irrepTmp = lp_malloc(numIrreps,3*numSites,false);
	for (i = 0; i < numIrreps; i++) {
	    for (j = 0; j < numSites; j++) {
		cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1.0,lp_ptr(LtoCnorm),lp_cols(LtoCnorm),&(lp_idx(irrep,i,3*j,LP_REAL)),1,0.0,
			&(lp_idx(irrepTmp,i,3*j,LP_REAL)),1);
	    }
	}
	
	/*
	 * At this point, we have everything needed calculate irrep contributions. The expensive matrix
	 * construction step is the phase matrix, which varies by wavevector. So we iterate over
	 * wavevectors (kx,ky,kz), construct the phase matrix, then compute the irCoeffs for all irreps at that
	 * wavevector.
	 *
	 * A clean implementation, written so that OpenMP can parallelize it. If needed in the future
	 * this will also allow for easy porting to MPI or similar.
	 */
        // Generate AtomstoIR matrix
        AtomstoIR = lp_malloc(numIrreps,numCells,false);

	// AtomstoIR = irrep * deltaMatrix^T
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,numIrreps,numCells,3*numSites,1.0,lp_ptr(irrepTmp),lp_cols(irrepTmp),
                lp_ptr(deltaMatrix),lp_cols(deltaMatrix),0.0,lp_ptr(AtomstoIR),lp_cols(AtomstoIR));

	// free unneeded matrices
	lp_free(tmpDM);
        lp_free(irrepTmp);
        lp_free(deltaMatrix);

	// Allocate result matrix if needed
	if (!(*irCoeffs))
		*irCoeffs = lp_malloc(numCells,numIrreps,true);
	else if (lp_rows((*irCoeffs)) != numCells || lp_cols((*irCoeffs)) != numIrreps)
		printErrorAndExit(12, "irrepDecompose: irCoeffs already allocated, but wrong size!\n");
	
	#pragma omp parallel default(shared) private(phaseMatrix,irMatrixSmall,k,kx,ky,kz,l,fj,tmp)
	{

	 // per-thread/instance arrays
	 phaseMatrix = lp_malloc(numCells,1,true);
	 irMatrixSmall = lp_malloc(3*numSites,1,true);

	 #pragma omp for
	 for (k = 0; k < numCells; k++) {
		kx = k / (cellsPerSupercell[1]*cellsPerSupercell[2]);
		ky = (k % (cellsPerSupercell[1]*cellsPerSupercell[2])) / cellsPerSupercell[2];
		kz = (k % (cellsPerSupercell[1]*cellsPerSupercell[2])) % cellsPerSupercell[2];

		// construct phase matrix
		l = 0;
		for (tmp[0] = 0; tmp[0] < cellsPerSupercell[0]; tmp[0]++) {
		   for (tmp[1] = 0; tmp[1] < cellsPerSupercell[1]; tmp[1]++) {
		      for (tmp[2] = 0; tmp[2] < cellsPerSupercell[2]; tmp[2]++) {
			fj = 2.0*PI*(FTYPE)kx*(FTYPE)tmp[0]/(FTYPE)cellsPerSupercell[0] + 2.0*PI*(FTYPE)ky*(FTYPE)tmp[1]/(FTYPE)cellsPerSupercell[1] + 2.0*PI*(FTYPE)kz*(FTYPE)tmp[2]/(FTYPE)cellsPerSupercell[2];                
			lp_idx(phaseMatrix,l,0,LP_REAL) = cos(fj);
                        lp_idx(phaseMatrix,l,0,LP_IMAG) = -1.0*sin(fj);
			l++;
		      }
		   }
		}

		// compute irMatrixSmall = AtomstoIR * phaseMatrix [result is numIrreps x 2 when considered as real]
		cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,numIrreps,2,numCells,1.0,lp_ptr(AtomstoIR),lp_cols(AtomstoIR),
			lp_ptr(phaseMatrix),2*lp_cols(phaseMatrix),0.0,lp_ptr(irMatrixSmall),2*lp_cols(irMatrixSmall));

		// copy to final output matrix
		for (l = 0; l < numIrreps; l++) {
		  lp_idx((*irCoeffs),k,l,LP_REAL) = lp_idx(irMatrixSmall,l,0,LP_REAL) / (FTYPE)numCells;
                  lp_idx((*irCoeffs),k,l,LP_IMAG) = lp_idx(irMatrixSmall,l,0,LP_IMAG) / (FTYPE)numCells;
		}
	 }

         // free per-thread/instance arrays
         lp_free(phaseMatrix);
         lp_free(irMatrixSmall);

	} // end omp parallel

	// Free memory 
	lp_free(AtomstoIR);
	safe_free(deltaMatrixUsed);
}

