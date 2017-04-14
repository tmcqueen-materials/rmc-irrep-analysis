bool IsIdentityMatrix(LAPACK *rot) {
        int i, j;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                if (i != j && fabs(lp_idx(rot,i,j,LP_REAL)) > ZERO_TOL) return false;
                if (i == j && fabs(lp_idx(rot,i,j,LP_REAL)-1.0) > ZERO_TOL) return false;
            }
        }
	return true;
}

bool IsIdentityTranslate(LAPACK *trans) {
        int i;
        for (i = 0; i < 3; i++) {
            if (fabs(mod1(lp_idx(trans,i,0,LP_REAL))) > ZERO_TOL) return false;
        }
	return true;
}

bool IsIdentityTransform(LAPACK *rot, LAPACK *trans) {
	if (!IsIdentityMatrix(rot)) return false;
	if (!IsIdentityTranslate(trans)) return false;
	return true;
}

void GenerateOrientations(int level, int *tmp) {
	int i;
	if (level >= numOrientationMatrix) {
		printStatus("%8i [", tmp[numOrientationMatrix]);
		for (i = 0; i < numOrientationMatrix; i++) {
			printStatus(" %4i", tmp[i]);
			orientations[tmp[numOrientationMatrix]][i] = tmp[i];
		}
		printStatus("]\n");
		tmp[numOrientationMatrix] += 1;
		return;
	}
	for (i = 0; i < orientToIdentity[level]; i++) {
		tmp[level] = i;
		GenerateOrientations(level+1,tmp);
	}
	return;
}

/*
 * Read in orientation matrices.
 */
void readOrientationFile(int f, char *file) {
	char *curLine;
	int i,j;
	LAPACK *t3x3 = NULL, *t3x32 = NULL, *t3x1 = NULL;
	int *tmp = NULL;

       #pragma omp single
       {
	openFile(f, file);
	/*
	 * Get rid of old matrices, etc, if they exist
	 */
	freeOrientationGlobals();

	// first line is number of matrices
	curLine = getNextLine(f);
	if (!curLine) printErrorAndExit(5, "readOrientationFile: Could not read first line from file %s\n", file);
	i = sscanf(curLine, " %i", &numOrientationMatrix);
	if (i < 1 || numOrientationMatrix < 1)
		printErrorAndExit(5, "readOrientationFile: Invalid or missing number of matrices in line %s: %i (expected > 0)\n", file,
			numOrientationMatrix);

	// allocate matrices
	orientMatrix = safe_malloc((UINT64)numOrientationMatrix*(UINT64)sizeof(LAPACK *));
	for (i = 0; i < numOrientationMatrix; i++) orientMatrix[i] = lp_malloc(3,3,false);
	orientShiftMatrix = safe_malloc((UINT64)numOrientationMatrix*(UINT64)sizeof(LAPACK *));	
	for (i = 0; i < numOrientationMatrix; i++) orientShiftMatrix[i] = lp_malloc(3,1,false);
        orientToIdentity = safe_malloc((UINT64)numOrientationMatrix*(UINT64)sizeof(int));

	// read in matrices
	for (i = 0; i < numOrientationMatrix; i++) {
	    for (j = 0; j < 3; j++) {
		curLine = getNextLine(f);
		if (!curLine) printErrorAndExit(5, "readOrientationFile: Premature end of file at %i,%i\n", i, j);
		if (sscanf(curLine, " %lf %lf %lf %lf", &(lp_idx(orientMatrix[i],j,0,LP_REAL)),
			&(lp_idx(orientMatrix[i],j,1,LP_REAL)), &(lp_idx(orientMatrix[i],j,2,LP_REAL)),
			&(lp_idx(orientShiftMatrix[i],j,0,LP_REAL))) < 4)
			printErrorAndExit(5, "readOrientationFile: Could not read four numbers from line %i,%i: %s\n", i, j, curLine);
	    }
	}

	// Now we have to figure out the number of times a given operation must be applied to generate the identity matrix.
	// We apply the following special rules:
	// 1) If the rotation matrix is the identity matrix, we go until the translate vector returns to identity
	// 2) If the rotation matrix is not identity, we go until the rotation matrix arrives at identity,
	//    EVEN IF the translate vector portion did not return to identity.
	// These two rules allow the orientations file to be constructed easily as just a list of the generators of a 
	// space group.
	t3x3 = lp_malloc(3,3,false);
	t3x32 = lp_malloc(3,3,false);
	t3x1 = lp_malloc(3,1,false);
	for (i = 0; i < numOrientationMatrix; i++) {
	    cblas_dcopy(9, lp_ptr(orientMatrix[i]), 1, lp_ptr(t3x3), 1);
	    cblas_dcopy(3, lp_ptr(orientShiftMatrix[i]), 1, lp_ptr(t3x1), 1);
	    orientToIdentity[i] = 1;
	    while ((IsIdentityMatrix(orientMatrix[i]) && !IsIdentityTranslate(t3x1)) ||
		   (!IsIdentityMatrix(orientMatrix[i]) && !IsIdentityMatrix(t3x3))) {
		if (orientToIdentity[i] > MAX_ORIENTATIONS_TO_IDENTITY)
			printErrorAndExit(5, "readOrientationFile: Orientation %i takes more than %i operations to return to identity.\n",
				i+1, MAX_ORIENTATIONS_TO_IDENTITY);
		cblas_daxpy(3, 1.0, lp_ptr(orientShiftMatrix[i]), 1, lp_ptr(t3x1), 1);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, lp_ptr(orientMatrix[i]),
				lp_cols(orientMatrix[i]), lp_ptr(t3x3), lp_cols(t3x3), 0.0,
				lp_ptr(t3x32), lp_cols(t3x32));
		cblas_dcopy(9, lp_ptr(t3x32), 1, lp_ptr(t3x3), 1);
		orientToIdentity[i] += 1;
	    }
	    if (i == 0 && (!IsIdentityTransform(t3x3,t3x1) || orientToIdentity[0] > 1))
		printErrorAndExit(5, "readOrientationFile: First orientation matrix is not identity!\n");
	    if (!IsIdentityMatrix(t3x3))
		printErrorAndExit(5, "readOrientationFile: After %i operations, orientation %i does not have an identity rotation matrix.\n",
			orientToIdentity[i], i+1);
	    if (!IsIdentityTranslate(t3x1))
		printStatus("readOrientationFile: After %i operations, the rotation matrix of orientation %i returns to identity, but the translate vector is not identity. Make sure this is the intended behavior (usually ok for centered space groups).\n", orientToIdentity[i], i+1);
	}
	lp_free(t3x1);
	lp_free(t3x32);
	lp_free(t3x3);

	// Now figure out number of orientations and allocate memory to store them
	numOrientations = orientToIdentity[0];
        for (i = 1; i < numOrientationMatrix; i++) numOrientations *= orientToIdentity[i];

	if (numOrientations > MAX_ORIENTATIONS)
	    printErrorAndExit(5, "numOrientations=%i > MAX_ORIENTATIONS=%i. Are you sure orientations.cfg only contains generators?", numOrientations, MAX_ORIENTATIONS);

	printStatus("Number of Orientations: %i. This should equal the multiplicity of the general position of the space group used.\n", numOrientations);

	orientations = safe_malloc((UINT64)numOrientations*(UINT64)sizeof(int *));
	for (i = 0; i < numOrientations; i++) orientations[i] = safe_malloc((UINT64)numOrientationMatrix*(UINT64)sizeof(int));

	// Recursively generate all possible orientation configurations
	printStatus("ORIENTATIONS\n");
	printStatus("-------- ---------------------------------------------------------------------------------\n");
	tmp = safe_malloc((UINT64)(numOrientationMatrix+1)*(UINT64)sizeof(int));
	tmp[numOrientationMatrix] = 0;
	GenerateOrientations(0,tmp);
	safe_free(tmp);
	
	closeFile(f);
       } // end omp single
}

void applyOrientation(int orient, LAPACK *atomMatrix, LAPACK **atomMatrixOutput) {
	int i,j,k,l;
	LAPACK *atomMatrixTemp = NULL;

	// check inputs/outputs
	if (*(atomMatrixOutput) == NULL)
		*atomMatrixOutput = lp_malloc(numAtoms,3,false);
	if (lp_rows((*atomMatrixOutput)) != numAtoms || lp_cols((*atomMatrixOutput)) != 3)
		printErrorAndExit(6, "applyOrientation: Invalid sized atomMatrixOutput array! %ix%i (expected %ix%i)\n",
			lp_rows((*atomMatrixOutput)), lp_cols((*atomMatrixOutput)), numAtoms, 3);
        if (lp_rows(atomMatrix) != numAtoms || lp_cols(atomMatrix) != 3)
                printErrorAndExit(6, "applyOrientation: Invalid sized atomMatrix array! %ix%i (expected %ix%i)\n",
                        lp_rows(atomMatrix), lp_cols(atomMatrix), numAtoms, 3);

	// allocate temporary storage
	atomMatrixTemp = lp_malloc(numAtoms,3,false);

	// Copy input matrix to output
	cblas_dcopy(lp_rows(atomMatrix)*lp_cols(atomMatrix), lp_ptr(atomMatrix), 1, lp_ptr((*atomMatrixOutput)), 1);

	// Apply all relevant transformations one at a time
	for (i = 0; i < numOrientationMatrix; i++) {
	  for (l = 0; l < orientations[orient][i]; l++) {

	    // Copy translation to all atoms (applied after rotation)	
	    for (j = 0; j < numAtoms; j++) 
	      for (k = 0; k < 3; k++)
		lp_idx(atomMatrixTemp,j,k,LP_REAL) = lp_idx(orientShiftMatrix[i],0,k,LP_REAL);

            // Apply rotation and add to translation already stored in temp matrix
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, numAtoms, 3, 3, 1.0, lp_ptr((*atomMatrixOutput)), lp_cols((*atomMatrixOutput)),
                lp_ptr(orientMatrix[i]), lp_cols(orientMatrix[i]), 1.0, lp_ptr(atomMatrixTemp), lp_cols(atomMatrixTemp));

	    // Copy temp matrix to output
	    cblas_dcopy(lp_rows(atomMatrixTemp)*lp_cols(atomMatrixTemp), lp_ptr(atomMatrixTemp), 1, lp_ptr((*atomMatrixOutput)), 1);
	  }
	}

	// free temporary storage
	lp_free(atomMatrixTemp);
}

int orderOrd(int ord) {
	return ord/NUM_SIGN_MAP;
}

int signOrd(int ord) {
	return mod(ord,NUM_SIGN_MAP);
}

COO MapIrrep(int ord, int irrep) {
        COO rv;
	int order = orderOrd(ord);
	int sidx = signOrd(ord);

	rv.ksign[0] = signMap[sidx][0];
	rv.ksign[1] = signMap[sidx][1];
	rv.ksign[2] = signMap[sidx][2];
        // apply order
        rv.ir = abs(irrepOrder[order][irrep])-1;
	if (irrepOrder[order][irrep]*signMap[sidx][irrepXYZ[irrep]-1] < 0) rv.ph = PI;
	else rv.ph = 0.0;
        rv.tph[0] = 0;
        rv.tph[1] = 0;
        rv.tph[2] = 0;
        return rv;
}

// offsetPhases and deltaPhasesPh are modified!
void ComputeDeltaPhasesPh(FTYPE offsets[3], int order, LAPACK *transPhases, LAPACK *deltaPhases, LAPACK **osf, LAPACK **dpp) {
	int m,i,k;
	COO mir;
	FTYPE fi;
	LAPACK *offsetPhases = *(osf);
	LAPACK *deltaPhasesPh = *(dpp);

	for (m = 0; m < numIrreps; m++) {
	    mir = MapIrrep(order,m);
	    for (i = 0; i < 3; i++)
		lp_idx(offsetPhases,i,m,LP_REAL) = (2.0*(offsets[i]-(FTYPE)mir.tph[i]))*PI/(FTYPE)cellsPerSupercell[i];
	}
	if (deltaPhases != deltaPhasesPh) cblas_dcopy(numUnique*numIrreps,lp_ptr(deltaPhases),1,lp_ptr(deltaPhasesPh),1);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, numUnique, numIrreps, 3, -1.0,
	    lp_ptr(transPhases), lp_cols(transPhases), lp_ptr(offsetPhases), lp_cols(offsetPhases),
	    1.0, lp_ptr(deltaPhasesPh), lp_cols(deltaPhasesPh));

	// WE DO NOT round numbers in deltaPhasesPh to the nearest multiple of 2*PI.
	return;
}

// deltaPhasesPh modified
FTYPE DeltaPhasesPhMod2PiWeight(LAPACK **dpp, bool weight, LAPACK *oWeight) {
	int m,k;
	LAPACK *deltaPhasesPh = *(dpp);

	// Make modulo 2*PI
	for (m = 0; m < numIrreps; m++) {
	    for (k = 0; k < numUnique; k++) {
	      if (weight) {
		lp_idx(deltaPhasesPh,k,m,LP_REAL) = lp_idx(oWeight,k,m,LP_REAL)*mod2Pi(lp_idx(deltaPhasesPh,k,m,LP_REAL));
	      } else { // no weight
		lp_idx(deltaPhasesPh,k,m,LP_REAL) = mod2Pi(lp_idx(deltaPhasesPh,k,m,LP_REAL));
	      }
	    }
	}

	// Compute sum of squares of differences
	return cblas_ddot(numUnique*numIrreps,lp_ptr(deltaPhasesPh),1,lp_ptr(deltaPhasesPh),1);
}
