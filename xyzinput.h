/*
 * Read ideal unit cell positions, the "xyzinput", into global variables
 */
void readXYZFile(int f, char *file) {
	int i;
	char *curLine;

	#pragma omp single
	{
	 // allocate memory
	 if (idealPosition != NULL) lp_free(idealPosition);
	 idealPosition = lp_malloc(MAX_IDEAL_SITES, 3, false);
	 openFile(f,file);
	 // we do not know how many ideal sites there are, read lines until end of file
	 numSites = 0;
	 for (curLine = getNextLine(f); curLine != NULL && numSites < MAX_IDEAL_SITES; numSites++, curLine = getNextLine(f)) {
		i = sscanf_str(curLine, idealPositionLabel[numSites], MAX_IDEAL_SITES_LABEL_LEN);
		if (i < 1) printWarning("Warning: Could not read label of site from line %i: %s\n", numSites, curLine);
		i = sscanf(curLine, " %*s %lf %lf %lf", &(lp_idx(idealPosition,numSites,0,LP_REAL)), 
			&(lp_idx(idealPosition,numSites,1,LP_REAL)), &(lp_idx(idealPosition,numSites,2,LP_REAL)));
		if (i != 3) printWarning("Warning: Only read %i items from line %i: %s (expected 4)\n", i+1, numSites, curLine);
	 } 
	 closeFile(f);
	 printStatus("Total of %i ideal sites read in:\n", numSites);
	 if (numSites >= MAX_IDEAL_SITES) printWarning("More ideal sites in xyz file than allowed by MAX_IDEAL_SITES!\n");
	 printStatus("Site   x          y          z          \n");
	 printStatus("------ ---------- ---------- -----------\n");
	 for (i = 0; i < numSites; i++)
		printStatus("%6s %10.6f %10.6f %10.6f\n", idealPositionLabel[i], lp_idx(idealPosition,i,0,LP_REAL),
			lp_idx(idealPosition,i,1,LP_REAL), lp_idx(idealPosition,i,2,LP_REAL));
	} //end single
}

void resetXYZInput() {
	int i,j;
	#pragma omp critical
	{
	 numSites = 0;
	 if (idealPosition != NULL) lp_free(idealPosition);
	 idealPosition = NULL;
	 for (i = 0; i < MAX_IDEAL_SITES; i++)
	   for (j = 0; j < MAX_IDEAL_SITES_LABEL_LEN; j++)
	    idealPositionLabel[i][j] = 0;
	} // end omp critical
}
