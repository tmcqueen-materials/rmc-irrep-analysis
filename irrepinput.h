/*
 * Read in irreducible representations, the "irrepinput", into global variables
 */
void readIrrepFile(int f, char *file) {
	int i,j,k,l,m;
	char *curLine;
	FTYPE irrepNorm;
	LAPACK *tmp;

        #pragma omp single
        {
         // allocate memory
         if (irrep != NULL) lp_free(irrep);
         irrep = lp_malloc(MAX_IRREPS, 3*MAX_IDEAL_SITES, false);
 	 /*
	  * Read in irreducible representations. They are stored, one per line,
	  * with (x,y,z) triples for each line in xyzFile in succession.
	  * Position 0 is the label for that irrep. Position 1 is an arbitrary
	  * number that is the same for all irreps that together form a degenerate set.
	  * We automatically normalize each irrep so that the norm (within a
	  * single unit cell) is unity, so that the coefficients can have
	  * real meaning.
	  */
	 openFile(f,file);
	 // erase any previous irrep information
	 for (i = 0; i < MAX_IRREPS; i++) {
	  for (j = 0; j < MAX_IRREPS+1; j++) {
		irrepByDegen[i][j] = -1;
	  }
	 }
         // we do not know how many irreps there are (although we expect numSites*3), read lines until end of file
	 numIrreps = 0;
         for (curLine = getNextLine(f); curLine != NULL && numIrreps < MAX_IRREPS; numIrreps++, curLine = getNextLine(f)) {
                i = sscanf_str(curLine, irrepLabel[numIrreps], MAX_IRREPS_LABEL_LEN);
                if (i < 1) printWarning("Warning: Could not read label of irrep from line %i: %s\n", numIrreps, curLine);
		else { k = 0; sscanf(curLine, " %*s%n", &k); curLine += k; }
		k = 0; m = -1;
		i = sscanf(curLine, " %i%n", &m, &k); curLine += k;
		if (i < 1) printWarning("Warning: Could not read degeneracy number of irrep from line %i: %s\n", numIrreps, curLine);
		else { // add to irreps-by-degen lookup table
		    l = 0;
		    while (irrepByDegen[l][0] != m && irrepByDegen[l][0] != -1) l++;
		    if (l >= MAX_IRREPS) { l = MAX_IRREPS-1;  printWarning("Warning: Overflow filling out irrepByDegen table!\n"); }
		    irrepByDegen[l][0] = m;
		    k = 1;
		    while (irrepByDegen[l][k] != -1) k++;
		    if (k > MAX_IRREPS) { k = MAX_IRREPS; printWarning("Warning: Overflow on k filling out irrepByDegen table!\n"); }
		    irrepByDegen[l][k] = numIrreps;
		}
		l = 0;
		do {
			k = 0;
			j = sscanf(curLine, " %lf %lf %lf%n", &(lp_idx(irrep,numIrreps,3*l+0,LP_REAL)), &(lp_idx(irrep,numIrreps,3*l+1,LP_REAL)),
				&(lp_idx(irrep,numIrreps,3*l+2,LP_REAL)), &k);
			curLine += k;
			if (j == 3) l++;
		} while (j == 3 && (l < numSites || numSites == 0));
                if (l != numSites && numSites > 0) 
			printWarning("Warning: Only read %i items from line %i: %s (expected %i)\n", 3*l+1, numSites, curLine, 3*numSites+1);
         }
	 closeFile(f);
	 numIrrepByDegen = l;
	 if (numIrreps != 3*numSites && numSites > 0) printWarning("Warning: numIrreps=%i is not equal to 3*numSites=%i\n", numIrreps, 3*numSites);
         printStatus("Total of %i irreps read in:\n", numIrreps);
         if (numIrreps >= MAX_IRREPS) printWarning("Warning: More irreps in file than allowed by MAX_IRREPS!\n");
         // cleanup irrep-by-degen 0th entries to number of irreps in each degeneracy set
         l = 0;
         while (irrepByDegen[l][0] != -1) {
            k = 0;
            while (irrepByDegen[l][k+1] != -1) k++;
            irrepByDegen[l][0] = k;
            l++;
         }
	 numIrrepByDegen = l;
	 printStatus("Irrep Degeneracy Table\n");
	 printStatus("Degen Irreps          \n");
	 printStatus("----- ----------------\n");
	 for (l = 0; l < numIrrepByDegen; l++) {
	    printStatus("%5i ", l);
	    for (k = 0; k < irrepByDegen[l][0]; k++)
		printStatus("%i ", irrepByDegen[l][k+1]);
	    printStatus("\n");
	 }
         printStatus("Irrep  Norm      ");
	 for (i = 0; i < numSites; i++) printStatus(" x          y          z          ");
	 printStatus("\n");
         printStatus("------ ----------");
	 for (i = 0; i < numSites; i++) printStatus(" ---------- ---------- -----------");
	 printStatus("\n");
	 tmp = lp_malloc(3,1,false);
	 if (LtoCnorm == NULL) printErrorAndExit(9, "LtoCnorm not defined! Has unit cell been specified yet?\n");
	 for (j = 0; j < numIrreps; j++) {
		// Normalize each irrep
		irrepNorm = 0.0;
		for (i = 0; i < numSites; i++) {
		    cblas_dgemv(CblasRowMajor, CblasNoTrans, lp_rows(LtoCnorm), lp_cols(LtoCnorm), 1.0, lp_ptr(LtoCnorm), lp_cols(LtoCnorm),
			&(lp_idx(irrep,j,3*i+0,LP_REAL)), 1, 0.0, lp_ptr(tmp), 1);
		    irrepNorm += cblas_ddot(3, lp_ptr(tmp), 1, lp_ptr(tmp), 1);
		}
		irrepNorm = sqrt(irrepNorm);
		for (i = 0; i < numSites; i++) {
			lp_idx(irrep,j,3*i+0,LP_REAL) /= irrepNorm;
                        lp_idx(irrep,j,3*i+1,LP_REAL) /= irrepNorm;
                        lp_idx(irrep,j,3*i+2,LP_REAL) /= irrepNorm;
		}
		printStatus("%6s %10.6f", irrepLabel[j], irrepNorm);
		for (i = 0; i < numSites; i++)
			printStatus(" %10.6f %10.6f %10.6f", lp_idx(irrep,j,3*i+0,LP_REAL), lp_idx(irrep,j,3*i+1,LP_REAL), lp_idx(irrep,j,3*i+2,LP_REAL));
		printStatus("\n");
	 }
	 lp_free(tmp);
	} // end omp single
}

/*
 * Read irrep order file into global variables
 */
void readIrrepOrderFile(int f, char *file) {
	char *curLine;
	int i,j,l,m,n;

       #pragma omp single
       {
        /*
         * Read in irreducible representation orders. Each order is stored in a column.
	 * The first three lines contain the arrangement of kx, ky, and kz wavevectors
	 * (only for degenerate line, zeros otherwise).
	 * Successive lines then contain the orders on a per-irrep basis. These
	 * must be in exactly the same order as the irrep.cfg file used to generate
	 * the input irrep coefficients.
	 *
         * Position 0 is the label for that irrep. Position 1 is an arbitrary number 
	 * the same for all irreps (or wavevectors) that together form a degenerate set.
	 * Position 2 is a number indicating whether it is an x (1), y (2), or z (3) component.
	 *
	 * After that is one or more columns of integers that specify how to order
	 * the irreps for different 'equivalent' choices. The first column is just
	 * a sequential set of integers (1,2,...,Nirrep). Additional columns give
	 * the ordering for each other 'equivalent' coordinate choice (if any). 
	 *
         * This code assumes that the a, b, and c basis vectors are orthonormal
         * (i.e. does not work for trigonal, hexagonal, monoclinic, triclinic !)
         * It can be easily modified to handle those cases, if one has the
         * basis vectors a, b, and c written in an orthonormal coordinate system
         * (like cartesian x,y,z).
         */
        openFile(f, file);
	printStatus(" Irrep   Degen  XYZ  Order0   Order1   Order2   Order3 \n");
	printStatus("------- ------- --- -------- -------- -------- --------\n");
	// first read in the kx,ky,kz degeneracy lines
	for (i = 0, curLine = getNextLine(f); curLine != NULL && i < 3; i++, curLine = getNextLine(f)) {
		j = 0; sscanf(curLine, " %*s%n", &j); curLine += j;
		// read degeneracy, XYZmap, and first ordering
		j = sscanf(curLine, " %i %i %*i%n", &(kvecDegen[i]), &m, &l);
		curLine += l;
		if (j < 2 || m-1 != i) printErrorAndExit(3, "Could not read first k-vector data from line %i: %s\n", i, curLine);
                // second+ order
                m = 1;
                while (sscanf(curLine, " %i%n", &j, &l) > 0 && m < MAX_IRREPS_ORDERS) {
                        m++;
                        curLine += l;
                }
                if (i == 0) numOrderSets = m;
                if (numOrderSets != m) printWarning("Incorrect number of orders from line %i (expected %i): %s\n", m, numOrderSets, curLine);

                printStatus("%6s%1i %7i %3s", "K", i, kvecDegen[i], "N/A");
                for (j = 0; j < numOrderSets && j < 4; j++) printStatus(" %8s", "N/A");
                printStatus("\n");
	}
        // we may not know how many irreps there are, read lines until end of file
        for (n = 0 /* read by previous for loop*/; curLine != NULL && n < MAX_IRREPS; n++, curLine = getNextLine(f)) {
                i = sscanf_str(curLine, irrepLabel[n], MAX_IRREPS_LABEL_LEN);
                if (i < 1) printWarning("Warning: Could not read label of irrep from line %i: %s\n", n+3, curLine);
                else { j = 0; sscanf(curLine, " %*s%n", &j); curLine += j; }
		// degeneracy, XYZMap, and first order
		j = sscanf(curLine, " %i %i %i%n", &(irrepDegen[n]), &(irrepXYZ[n]), &(irrepOrder[0][n]), &l);
		curLine += l;
                if (j < 3) printErrorAndExit(3, "Could not read irrep degeneracy or first order from line %i: %s\n", n+3, curLine);
		if (irrepXYZ[n] < 1 || irrepXYZ[n] > 3)
			printErrorAndExit(3, "XYZ Map numbers must be x=1, y=2, or z=3, not %i!\n", irrepXYZ[n]);
		if ((irrepOrder[0][n]-1) != n)
			printErrorAndExit(3, "Order 0 numbers do not sequentially increase by 1, starting at 1!\n");
		// second+ order
		i = 1;
		while (sscanf(curLine, " %i%n", &(irrepOrder[i][n]), &l) > 0 && i < MAX_IRREPS_ORDERS) {
			i++;
			curLine += l;
		}
		if (numOrderSets != i) printWarning("Incorrect number of orders from line %i: %i (expected %i): %s\n", n+3, i, numOrderSets, curLine);
		printStatus("%7s %7i %3i", irrepLabel[n], irrepDegen[n], irrepXYZ[n]);
		for (i = 0; i < numOrderSets && i < 4; i++) printStatus(" %8i", irrepOrder[i][n]);
		printStatus("\n");
        }
        printStatus("Total of %i orders. Up to the first four shown.\n", numOrderSets);
	if (n != numIrreps && numIrreps > 0) printErrorAndExit(5, "Number of irreps in order file does not match number of irreps from irrep file.\n");
	numIrreps = n;
	// Sanity check each order by making sure they sum to the expected value
	for (i = 0; i < numOrderSets; i++) {
		l = 0;
		for (j = 0; j < numIrreps; j++)
			l += abs(irrepOrder[i][j]);
		if (l != ((numIrreps+1)*numIrreps)/2)
			printWarning("Order %i has a mistake!\n", i);
	}
	numOrders = numOrderSets*NUM_SIGN_MAP;
        closeFile(f);
       } // end omp single
}

