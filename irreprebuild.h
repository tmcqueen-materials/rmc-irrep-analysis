/*
 * Coverts the passed set of irrep coefficients into an atomMatrix with a scale factor.
 * Assumes ideal unit cell xyz's, and irrep atom specifiers, have been read into global variables.
 */

void irrepRebuild(FTYPE scale, LAPACK *irCoeffs, LAPACK **atomMatrix, int **atomMatrixLabels) {
	int i,j,k,l,m,n,o = 0;
	int subcell = 0;
	int kx,ky,kz = 0;
	FTYPE fi = 0.0;
	
	LAPACK *IRtoAtoms = NULL; // 3*numSites x numCells matrix (complex) [= irrepComplex^T * irCoeffs^T]
	LAPACK *irrepComplex = NULL; // numIrreps x 3*numSites matrix (complex)

	// compute matrices
	int curAtom = 0;
	LAPACK *phaseMatrix = NULL; // numCells x 1 (=e^ikx for all k for a single x) matrix (complex)
	LAPACK *deltaMatrixSmall = NULL; // 3*numSites x 1 matrix (complex) [all site coordinates for a single subcell (x value)]

	// results
        LAPACK *deltaMatrix = NULL; // numAtoms x 3 matrix (real)
        LAPACK *idealMatrix = NULL; // numAtoms x 3 matrix (real)
	// atomMatrix: numAtoms x 3 matrix (real)
	// atomMatrixLabels: numAtoms matrix (real)

	/*
	 * At this point, we have read everything needed to rebuild from irreps. The expensive matrix
	 * construction step is the phase matrix, which varies by subcell index. So we iterate over
	 * subcell indices (i,j,k), construct the phase matrix, then compute the deltaMatrix for all sites in the
	 * subcell.
	 *
	 * A clean implementation, written so that OpenMP can parallelize it. If needed in the future
	 * this will also allow for easy porting to MPI or similar.
	 */

        // Generate IRtoAtoms matrix
        IRtoAtoms = lp_malloc(3*numSites,numCells,true);
	irrepComplex = lp_malloc(numIrreps,3*numSites,true);

	// copy irrep to irrepComplex (all imaginary components = 0), simultaneously converting to cartesian, normalized, coordinates
	for (i = 0; i < numIrreps; i++) {
	    for (j = 0; j < numSites; j++) {
		cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, lp_ptr(LtoCnorm), lp_cols(LtoCnorm), &(lp_idx(irrep,i,3*j,LP_REAL)), 1,
		    0.0, &(lp_idx(irrepComplex,i,3*j,LP_REAL)), 2);
	    }
	}

	// IRtoAtoms = irrepComplex^T * irCoeffs^T
        cblas_zgemm(CblasRowMajor,CblasTrans,CblasTrans,3*numSites,numCells,numIrreps,CMPLX1p0i,lp_ptr(irrepComplex),lp_cols(irrepComplex),
                lp_ptr(irCoeffs),lp_cols(irCoeffs),CMPLX0p0i,lp_ptr(IRtoAtoms),lp_cols(IRtoAtoms));
	// free irrepComplex
	lp_free(irrepComplex);

	// Allocate other shared matrices
        deltaMatrix = lp_malloc(numAtoms,3,false);
        idealMatrix = lp_malloc(numAtoms,3,false);
        if (!(*atomMatrix)) 
		*atomMatrix = lp_malloc(numAtoms,3,false);
	else if (lp_rows((*atomMatrix)) != numAtoms || lp_cols((*atomMatrix)) != 3) 
		printErrorAndExit(12, "irrepRebuild: Atom matrix passed in as allocated, but is incorrect size!\n");
	if (!(*atomMatrixLabels))
        	*atomMatrixLabels = safe_malloc((UINT64)numAtoms*(UINT64)sizeof(int));

        // generate ideal atom positions
        curAtom = 0;
        for (i = 0; i < cellsPerSupercell[0]; i++) {
          for (j = 0; j < cellsPerSupercell[1]; j++) {
            for (k = 0; k < cellsPerSupercell[2]; k++) {
                for (l = 0; l < numSites; l++) {
                    (*atomMatrixLabels)[curAtom] = l;
                    lp_idx(idealMatrix, curAtom, 0,LP_REAL) = 1.0*lp_idx(idealPosition,l,0,LP_REAL)+(FTYPE)i;
                    lp_idx(idealMatrix, curAtom, 1,LP_REAL) = 1.0*lp_idx(idealPosition,l,1,LP_REAL)+(FTYPE)j;
                    lp_idx(idealMatrix, curAtom, 2,LP_REAL) = 1.0*lp_idx(idealPosition,l,2,LP_REAL)+(FTYPE)k;
                    curAtom++;
                }
            }
          }
        }

	#pragma omp parallel default(shared) private(phaseMatrix,deltaMatrixSmall,subcell,i,j,k,m,l,kx,ky,kz,fi,o,n,curAtom)
	{

	 // per-thread/instance arrays
	 phaseMatrix = lp_malloc(numCells,1,true);
	 deltaMatrixSmall = lp_malloc(3*numSites,1,true);

	 #pragma omp for
	 for (subcell = 0; subcell < numCells; subcell++) {
		i = subcell / (cellsPerSupercell[1]*cellsPerSupercell[2]);
		j = (subcell % (cellsPerSupercell[1]*cellsPerSupercell[2])) / cellsPerSupercell[2];
		k = (subcell % (cellsPerSupercell[1]*cellsPerSupercell[2])) % cellsPerSupercell[2];

		// construct phase matrix
		m = 0;
		for (kx = 0; kx < cellsPerSupercell[0]; kx++) {
		   for (ky = 0; ky < cellsPerSupercell[1]; ky++) {
		      for (kz = 0; kz < cellsPerSupercell[2]; kz++) {
			fi = 2.0*PI*(FTYPE)kx*(FTYPE)i/(FTYPE)cellsPerSupercell[0] + 2.0*PI*(FTYPE)ky*(FTYPE)j/(FTYPE)cellsPerSupercell[1] + 2.0*PI*(FTYPE)kz*(FTYPE)k/(FTYPE)cellsPerSupercell[2];                
			lp_idx(phaseMatrix,m,0,LP_REAL) = cos(fi);
                        lp_idx(phaseMatrix,m,0,LP_IMAG) = sin(fi);
			m++;
		      }
		   }
		}

		// compute deltaMatrxiSmall = IRtoAtoms * phaseMatrix [result is 3*numSites x 1]
		cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3*numSites,1,numCells,CMPLX1p0i,lp_ptr(IRtoAtoms),lp_cols(IRtoAtoms),
			lp_ptr(phaseMatrix),lp_cols(phaseMatrix),CMPLX0p0i,lp_ptr(deltaMatrixSmall),lp_cols(deltaMatrixSmall));

		// check that all imaginary contributions vanished, and copy to final output matrix after converting back to lattice units
		for (l = 0; l < numSites; l++) {
		  curAtom = subcell*numSites+l;
		  for (o = 0; o < 3; o++) {
		    if (fabs(lp_idx(deltaMatrixSmall,3*l+o,0,LP_IMAG)) > ZERO_TOL) {
                        printWarning("irrepRebuild Warning: imaginary component for (%i,%i,%i,%i) not zero (component %i): %f > %f\n",i,j,k,l,o,
                                fabs(lp_idx(deltaMatrixSmall,3*l+o,0,LP_IMAG)), ZERO_TOL);
		    }
		  }
		  cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, lp_ptr(CtoL), lp_cols(CtoL), &(lp_idx(deltaMatrixSmall,3*l+0,0,LP_REAL)), 2,
			0.0, &(lp_idx(deltaMatrix,curAtom,0,LP_REAL)), 1);
		}
	 }

         // free per-thread/instance arrays
         lp_free(phaseMatrix);
         lp_free(deltaMatrixSmall);

	} // end omp parallel

	/*
	 * At this point, idealMatrix contains ideal atom positions, atomMatrixLabels contain indices
	 * to idealPositionLabel entries, and deltaMatrix contains the displacement from ideal from the
	 * sum of irrep contributions. 
	 *
	 * now compute atomMatrix = idealMatrix+scale*deltaMatrix, where scale is an arbitrary factor
	 * of how much to scale the displacement. If you want to compute multiple results with different
	 * scale factors, save lots of CPU time by just computing atomMatrix many times!
	 *
	 * in the daxpy call, the second parameter is the scale factor (1.00)
	 */
	cblas_dcopy(3*numAtoms, lp_ptr(idealMatrix), 1, lp_ptr((*atomMatrix)), 1);
	cblas_daxpy(3*numAtoms, scale, lp_ptr(deltaMatrix), 1, lp_ptr((*atomMatrix)), 1);
	// cleanup to make sure the positions stay inside the big box
	for (i = 0; i < numAtoms; i++) {
	    for (j = 0; j < 3; j++) {
		if (lp_idx((*atomMatrix),i,j,LP_REAL) < 0.0) {
		    lp_idx((*atomMatrix),i,j,LP_REAL) += (FTYPE)cellsPerSupercell[j];
		}
		if (lp_idx((*atomMatrix),i,j,LP_REAL) >= (FTYPE)cellsPerSupercell[j]) {
		    lp_idx((*atomMatrix),i,j,LP_REAL) -= (FTYPE)cellsPerSupercell[j];
		}
	    }
	}

	// Free memory 
        lp_free(idealMatrix);
	lp_free(deltaMatrix);
	lp_free(IRtoAtoms);
}

