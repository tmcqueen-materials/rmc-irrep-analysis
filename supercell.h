/*    For a given set of supercell indices (x,y,z or kx,ky,kz), compure the appropriate index in a
 *    full, but 1-D, representation.
 */
#define cellIndex(x,y,z) (mod((x),cellsPerSupercell[0])*cellsPerSupercell[1]*cellsPerSupercell[2] + mod((y),cellsPerSupercell[1])*cellsPerSupercell[2] + mod((z),cellsPerSupercell[2]))

/*
 * Give vector lengths and relative angles given a transform matrix
 */
#define magCol(ltc,c) ( cblas_dnrm2(lp_rows(ltc), &(lp_idx(ltc,0,c,LP_REAL)),lp_cols(ltc)) )
#define angCol(ltc,c1,c2) ( acos(cblas_ddot(lp_rows(ltc),&(lp_idx(ltc,0,c1,LP_REAL)),lp_cols(ltc),&(lp_idx(ltc,0,c2,LP_REAL)),lp_cols(ltc))/magCol(ltc,c1)/magCol(ltc,c2)) )

/*
 * Given two coordinates in relative lattice units, output the minimum distance between those two points,
 * obeying periodic boundary conditions. C = H*L, C an orthogonal coordinate system in desired units (here, angstroms)
 */
FTYPE rluMinDistance(FTYPE x1, FTYPE y1, FTYPE z1, FTYPE x2, FTYPE y2, FTYPE z2, LAPACK *H) {
	FTYPE dx = modp5(x1-x2);
	FTYPE dy = modp5(y1-y2);
	FTYPE dz = modp5(z1-z2);
	FTYPE t = 0.0;
	FTYPE rv = 0.0;
	int i;

	for (i = 0; i < 3; i++) {
	    t = lp_idx(H,i,0,LP_REAL)*dx + lp_idx(H,i,1,LP_REAL)*dy + lp_idx(H,i,2,LP_REAL)*dz;
	    rv += t*t;
	}

	return sqrt(rv);
}

/*
 * Print (sub/super)cell information
 */
void printCellsAndMatrices() {
	int i,j;
	printStatus("Supercell\n");
	printStatus("---------\n");
	printStatus("a = %12.5f  alpha = %9.3f\n", supercellEdge[0], cellAngles[0]*180.0/PI);
        printStatus("b = %12.5f   beta = %9.3f\n", supercellEdge[1], cellAngles[1]*180.0/PI);
        printStatus("c = %12.5f  gamma = %9.3f\n", supercellEdge[2], cellAngles[2]*180.0/PI);
	printStatus("Vol = %.5f\n", lp_det(ssLtoC));
	printStatus("ssLtoC = [");
	for (i = 0; i < 3; i++) {
	    printStatus("\n[ ");
	    for (j = 0; j < 3; j++)
		printStatus("%14.8e ", lp_idx(ssLtoC,i,j,LP_REAL));
	    printStatus(" ]");
	}
	printStatus("]\n");
        printStatus("ssCtoL = [");
        for (i = 0; i < 3; i++) {
            printStatus("\n[ ");
            for (j = 0; j < 3; j++)
                printStatus("%14.8e ", lp_idx(ssCtoL,i,j,LP_REAL));
            printStatus(" ]");
        }
        printStatus("]\n");
	printStatus("\nRecomputed Cell: [ %12.5f %12.5f %12.5f ] [ %9.3f %9.3f %9.3f ]\n\n", magCol(ssLtoC,0), magCol(ssLtoC,1), magCol(ssLtoC,2),
		angCol(ssLtoC,1,2)*180.0/PI, angCol(ssLtoC,0,2)*180.0/PI, angCol(ssLtoC,0,1)*180.0/PI);

        printStatus("\nSubcell\n");
        printStatus("---------\n");
        printStatus("a = %12.5f   alpha = %9.3f\n", cellEdge[0], cellAngles[0]*180.0/PI);
        printStatus("b = %12.5f    beta = %9.3f\n", cellEdge[1], cellAngles[1]*180.0/PI);
        printStatus("c = %12.5f   gamma = %9.3f\n", cellEdge[2], cellAngles[2]*180.0/PI);
        printStatus("Vol = %.5f\n", lp_det(LtoC));
        printStatus("LtoC = [");
        for (i = 0; i < 3; i++) {
            printStatus("\n[ ");
            for (j = 0; j < 3; j++)
                printStatus("%14.8e ", lp_idx(LtoC,i,j,LP_REAL));
            printStatus(" ]");
        }
        printStatus("]\n");
        printStatus("CtoL = [");
        for (i = 0; i < 3; i++) {
            printStatus("\n[ ");
            for (j = 0; j < 3; j++)
                printStatus("%14.8e ", lp_idx(CtoL,i,j,LP_REAL));
            printStatus(" ]");
        }
        printStatus("]\n");
        printStatus("\nRecomputed Cell: [ %12.5f %12.5f %12.5f ] [ %9.3f %9.3f %9.3f ]\n\n", magCol(LtoC,0), magCol(LtoC,1), magCol(LtoC,2),
                angCol(LtoC,1,2)*180.0/PI, angCol(LtoC,0,2)*180.0/PI, angCol(LtoC,0,1)*180.0/PI);

	return;
}

/* 
 * Construct LtoC and CtoL matrices given cell edges and angles
 */
void BuildLtoCandCtoL() {
	FTYPE V = 0.0;
	FTYPE fi;
	int i,j;

	/* Check inputs */
        if (cellsPerSupercell[0] <= 0)
                printErrorAndExit(6, "Cells per supercell x = %i is invalid\n", cellsPerSupercell[0]);
        if (cellsPerSupercell[1] <= 0)
                printErrorAndExit(6, "Cells per supercell y = %i is invalid\n", cellsPerSupercell[1]);
        if (cellsPerSupercell[2] <= 0)
                printErrorAndExit(6, "Cells per supercell z = %i is invalid\n", cellsPerSupercell[2]);
        if (supercellEdge[0] < ZERO_TOL)
                printErrorAndExit(6, "Supercell edge x = %.8f is invalid\n", supercellEdge[0]);
        if (supercellEdge[1] < ZERO_TOL)
                printErrorAndExit(6, "Supercell edge y = %.8f is invalid\n", supercellEdge[1]);
        if (supercellEdge[2] < ZERO_TOL)
                printErrorAndExit(6, "Supercell edge z = %.8f is invalid\n", supercellEdge[2]);
        if (cellAngles[0] < ZERO_TOL)
                printErrorAndExit(6, "Supercell angle alpha = %.8f (radians) is invalid\n", cellAngles[0]);
        if (cellAngles[1] < ZERO_TOL)
                printErrorAndExit(6, "Supercell angle beta = %.8f (radians) is invalid\n", cellAngles[1]);
        if (cellAngles[2] < ZERO_TOL)
                printErrorAndExit(6, "Supercell angle gamma = %.8f (radians) is invalid\n", cellAngles[2]);

	/* allocate matrices if needed */
	if (LtoC == NULL) { LtoC = lp_malloc(3,3,false); }
	if (ssLtoC == NULL) { ssLtoC = lp_malloc(3,3,false); }
	if (CtoL == NULL) { CtoL = lp_malloc(3,3,false); }
	if (ssCtoL == NULL) { ssCtoL = lp_malloc(3,3,false); }
	if (LtoCnorm == NULL) { LtoCnorm = lp_malloc(3,3,false); }

	/* Build LtoC matrix first (from Int. Tabl. Cryst. Vol. C) */
	lp_idx(LtoC,0,0,LP_REAL) = cellEdge[0];
	lp_idx(ssLtoC,0,0,LP_REAL) = supercellEdge[0];
        lp_idx(LtoC,1,0,LP_REAL) = 0.0;
        lp_idx(ssLtoC,1,0,LP_REAL) = 0.0;
        lp_idx(LtoC,2,0,LP_REAL) = 0.0;
        lp_idx(ssLtoC,2,0,LP_REAL) = 0.0;
        lp_idx(LtoC,0,1,LP_REAL) = cellEdge[1]*cos(cellAngles[2]);
        lp_idx(ssLtoC,0,1,LP_REAL) = supercellEdge[1]*cos(cellAngles[2]);
        lp_idx(LtoC,1,1,LP_REAL) = cellEdge[1]*sin(cellAngles[2]);
        lp_idx(ssLtoC,1,1,LP_REAL) = supercellEdge[1]*sin(cellAngles[2]);
        lp_idx(LtoC,2,1,LP_REAL) = 0.0;
        lp_idx(ssLtoC,2,1,LP_REAL) = 0.0;
        lp_idx(LtoC,0,2,LP_REAL) = cellEdge[2]*cos(cellAngles[1]);
        lp_idx(ssLtoC,0,2,LP_REAL) = supercellEdge[2]*cos(cellAngles[1]);
        lp_idx(LtoC,1,2,LP_REAL) = cellEdge[2]*(cos(cellAngles[0])-cos(cellAngles[1])*cos(cellAngles[2]))/sin(cellAngles[2]);
        lp_idx(ssLtoC,1,2,LP_REAL) = supercellEdge[2]*(cos(cellAngles[0])-cos(cellAngles[1])*cos(cellAngles[2]))/sin(cellAngles[2]);
	V = cellEdge[0]*cellEdge[1]*cellEdge[2]*sqrt(1.0-cos(cellAngles[0])*cos(cellAngles[0])-cos(cellAngles[1])*cos(cellAngles[1])
		-cos(cellAngles[2])*cos(cellAngles[2])+2.0*cos(cellAngles[0])*cos(cellAngles[1])*cos(cellAngles[2]));
        lp_idx(LtoC,2,2,LP_REAL) = V/(cellEdge[0]*cellEdge[1]*sin(cellAngles[2]));
	V /= cellEdge[0]*cellEdge[1]*cellEdge[2];
	V *= supercellEdge[0]*supercellEdge[1]*supercellEdge[2];
        lp_idx(ssLtoC,2,2,LP_REAL) = V/(supercellEdge[0]*supercellEdge[1]*sin(cellAngles[2]));
	
	/* Prepare CtoL by inverting LtoC */
	lp_invert(LtoC, &(CtoL));
	lp_invert(ssLtoC, &(ssCtoL));

	/* Prepare LtoCnorm by normalizing LtoC */
	fi = 0.0;
	cblas_dcopy(9,&fi,0,lp_ptr(LtoCnorm),1);
	cblas_daxpy(9,1.0/lp_det(LtoC),lp_ptr(LtoC),1,lp_ptr(LtoCnorm),1);

	/* Print Results, then Sanity Check */
	if (!LtoCandCtoLprinted) printCellsAndMatrices();
	LtoCandCtoLprinted = true;
	if (fabs(lp_det(ssLtoC)-V) > ZERO_TOL)
		printErrorAndExit(8, "Error: ssLtoC V=%.8f, expected V=%.8f\n", lp_det(ssLtoC),V);
	V /= supercellEdge[0]*supercellEdge[1]*supercellEdge[2];
	V *= cellEdge[0]*cellEdge[1]*cellEdge[2];
	if (fabs(lp_det(LtoC)-V) > ZERO_TOL)
		printErrorAndExit(8, "Error: LtoC V=%.8f, expected V=%.8f\n", lp_det(LtoC),V);
	for (i = 0; i < 3; i++) {
	    if (fabs(magCol(ssLtoC,0)-supercellEdge[0]) > ZERO_TOL)
		printErrorAndExit(8, "Error: Mismatch of vector %i in ssLtoC!\n", i);
            if (fabs(magCol(LtoC,0)-cellEdge[0]) > ZERO_TOL)
                printErrorAndExit(8, "Error: Mismatch of vector %i in LtoC!\n", i);
	    for (j = i+1; j < 3; j++) {
              if (fabs(angCol(ssLtoC,i,j)-cellAngles[3-i-j]) > ZERO_TOL)
                printErrorAndExit(8, "Error: Mismatch of angle (%i,%i) in ssLtoC!\n", i, j);
              if (fabs(angCol(LtoC,i,j)-cellAngles[3-i-j]) > ZERO_TOL)
                printErrorAndExit(8, "Error: Mismatch of angle (%i,%i) in LtoC!\n", i, j);
	    }
	}

	return;	
}

/* For speed, in return values, a (idx, sgn) pair are stored as sgn*idx (since sgn = +/-1)
 * These macros make the conversions transparent and automatic.
 */
#define CRidx int
#define CRidx_rv(idx, sgn) ((sgn)*(idx))
#define CRidx_idx(rv) (abs((rv)))
#define CRidx_sgn(rv) (((rv) >= 0) - ((rv) < 0))

/*
 * The range of unique wavevectors (frequencies) for the spectra are:
 * k_x = 0, ..., floor(Nx/2)   [*2*Pi/Nx]
 * k_y = 0, ..., floor(Ny/2)   [*2*Pi/Ny] if k_x = 0
 *     = 0, ..., Ny-1          [*2*Pi/Ny] otherwise
 * k_z = 0, ..., floor(Nz/2)   [*2*Pi/Nz] if k_x = 0 and k_y = 0
 *     = 0, ..., Nz-1          [*2*Pi/Nz] otherwise
 * Which corresponds to a total of Nunique = int(Nz/2)+1 + int(Ny/2)*Nz + int(Nx/2)*Ny*Nz
 * values when Nx, Ny, and Nz are odd. When Nx, Ny, and/or Nz are even, there are slightly fewer unique
 * components, but it requires complicated accounting to keep the variable lengths straight.
 * So we restrict ourselves to two representations:
 *    (i) "Full" representations, running from kx = 0...Nx-1, ky = 0...Ny-1, kz = 0...Nz-1
 *   (ii) "Compacted" representations, containing sequentially the frequency ranges above
 *        in a single 1-D array. This will exactly include all unique values when Nx,Ny,Nz
 *        are odd, and include only a few duplicates when one or more are even.
 * The bispectrum we use contains both auto (within irrep) and cross (between irrep) components,
 *  i.e. F*(irrep1,k1+k2)*F(irrep2,k2)*F(irrep1,k1).
 * To include all unique components, k1 spans a "compacted", unique set (ii), while k2 then
 *  spans a full set (i). We define this as the "compacted" bispectrum.
 * Here we have a series of index 'helper' functions to convert between these two representations.
 * The necessary transformations are:
 */

/*
 * 0. For a given kx,ky,kz, compure the appropriate index in a full, but 1-D, representation.
 *    The two components are directly related, so there is no sign to return.
 */
#define FtofullF(x,y,z) (cellIndex(x,y,z))

/*
 * 1. For a given full, but 1-D, index, compute the corresponding kx,ky,kz
 *    index. They are always directly related, so the sign relation is always +1.
 */
FRk fullFtoF(int idx) {
	FRk rv = {.kx = 0, .ky = 0, .kz = 0};
	rv.kx = idx/(cellsPerSupercell[1]*cellsPerSupercell[2]);
	rv.ky = (idx-rv.kx*cellsPerSupercell[1]*cellsPerSupercell[2])/cellsPerSupercell[2];
	rv.kz = (idx-rv.kx*cellsPerSupercell[1]*cellsPerSupercell[2]-rv.ky*cellsPerSupercell[2]);
	if (rv.kx < 0 || rv.kx > cellsPerSupercell[0]-1 || rv.ky < 0 || rv.ky > cellsPerSupercell[1]-1 || rv.kz < 0 || rv.kz > cellsPerSupercell[2]-1)
		printErrorAndExit(5, "Error in fullFtoF! idx=%i gives invalid value(s): k=(%i,%i,%i)\n", idx, rv.kx, rv.ky, rv.kz);
	return rv;
}

/*
 * 2. For a given kx,ky,kz, compute the appropriate index in the compacted representation.
 *    The two components may be conjugates, so we return the appropriate sign (+1 or -1)
 *    to complete the mapping.
 */
CRidx FtoargF(int x, int y, int z) {
	int kx = mod(x,cellsPerSupercell[0]);
	int ky = mod(y,cellsPerSupercell[1]);
	int kz = mod(z,cellsPerSupercell[2]);
	int sgn = 1;
	int idx = 0;
	if (kx == 0 && ky == 0 && ( kz/((cellsPerSupercell[2]/2)+1) ) > 0) {
		kz = mod(abs(kz-cellsPerSupercell[2]),cellsPerSupercell[2]);
		sgn = -1;
	} else if (kx == 0 && ( ky/((cellsPerSupercell[1]/2)+1) ) > 0) {
		ky = mod(abs(ky-cellsPerSupercell[1]),cellsPerSupercell[1]);
		kz = mod(abs(kz-cellsPerSupercell[2]),cellsPerSupercell[2]);
		sgn = -1;
	} else if ( ( kx/((cellsPerSupercell[0]/2)+1) ) > 0) {
		kx = mod(abs(kx-cellsPerSupercell[0]),cellsPerSupercell[0]);
		ky = mod(abs(ky-cellsPerSupercell[1]),cellsPerSupercell[1]);
		kz = mod(abs(kz-cellsPerSupercell[2]),cellsPerSupercell[2]);
		sgn = -1;
	}
	idx = MAX(0,kx-1)*cellsPerSupercell[1]*cellsPerSupercell[2] + MIN(1,kx)*((cellsPerSupercell[2]/2)+1+(cellsPerSupercell[1]/2)*cellsPerSupercell[2]);
	if (kx == 0) { idx += MAX(0,ky-1)*cellsPerSupercell[2] + MIN(1,ky)*((cellsPerSupercell[2]/2)+1); }
	else { idx += ky*cellsPerSupercell[2]; }
	idx += kz;
	return CRidx_rv(idx, sgn);
}

/*
 * 3. For a given index in the compacted representation, compute the corresponding kx,ky,kz
 *    index. We always choose the one that is directly related, not the conjugate, so the
 *    sign relation is always +1.
 */
FRk argFtoF(int idx) {
	FRk rv = {.kx = 0, .ky = 0, .kz = 0};
	int idx2, idx3;
	if (idx < ((cellsPerSupercell[2]/2)+1)) { rv.kz = idx; return rv; }
	idx2 = idx - ((cellsPerSupercell[2]/2)+1);
	if (idx2 < ((cellsPerSupercell[1]/2)*cellsPerSupercell[2])) { rv.kz = mod(idx2,cellsPerSupercell[2]); rv.ky = 1+(idx2/cellsPerSupercell[2]); return rv; }
	idx3 = idx2 - ((cellsPerSupercell[1]/2)*cellsPerSupercell[2]);
	rv.kx = 1+(idx3/(cellsPerSupercell[1]*cellsPerSupercell[2]));
	rv.ky = (mod(idx3,(cellsPerSupercell[1]*cellsPerSupercell[2]))/cellsPerSupercell[2]);
	rv.kz = mod(mod(idx3,(cellsPerSupercell[1]*cellsPerSupercell[2])),cellsPerSupercell[2]);
	return rv;
}

FRk argFtoFMapped(int idx, COO map) {
	FRk rv = argFtoF(idx);
	rv.kx *= map.ksign[0];
	rv.ky *= map.ksign[1];
	rv.kz *= map.ksign[2];
	return rv;
}

/* 
 * 4. For a given index to the compacted representation of the bispectrum, compute the
 *    corresponding kx1,ky1,kz1 , kx2,ky2,kz2 indices. As with 2, we always choose the one
 *    that is directly related, not a conjugate.
 */
FRk2 argBtoB(int idx) {
	int idx1 = (idx/numCells);
	int idx2 = mod(idx,numCells);
	FRk2 rv;
	rv.k1 = argFtoF(idx1);
	rv.k2 = fullFtoF(idx2);
	return rv;
}

/*
 * 5. For a given pair of kx1,ky1,kz1 , kx2,ky2,kz2 indices, compute the index in the
 *    compacted representation of the bispectrum. The two components may be conjugates, 
 *    so we return the appropriate sign (+1 or -1) to complete the mapping.
 */
CRidx BtoargB(x1,y1,z1,x2,y2,z2) {
	CRidx k1 = FtoargF(x1,y1,z1);
	int k2 = FtofullF(x2,y2,z2);
	int idx = CRidx_idx(k1)*numCells+k2;

	if (CRidx_sgn(k1) > 0) { return CRidx_rv(idx, 1); }
	// otherwise, use B(k1,k2) = B*(-k1,-k2)
	k1 = FtoargF(-x1,-y1,-z1);
	k2 = FtofullF(-x2,-y2,-z2);
	idx = CRidx_idx(k1)*numCells+k2;
        if (CRidx_sgn(k1) > 0) { return CRidx_rv(idx, -1); }

	/* Should never get here! */
	printErrorAndExit(5, "Error in BtoargB! No solution found: %i,(%i,%i,%i),%i,(%i,%i,%i)!\n", k1,x1,y1,z1,k2,x2,y2,z2);
	/* satisfy return requirement (we never get here due to previous line) */
	return CRidx_rv(0,0);
}

/*
 * 6. When cellsPerSupercell[0], cellsPerSupercell[1], and/or cellsPerSupercell[2] are even, the compacted representation still contains a few duplicates
 * Here is a function that checks whether the given entry is a duplicate and, if so, returns
 * the index of the non-duplicate entry
 */

CRidx argFunique(aFs) {
	if (mod(cellsPerSupercell[0],2) != 0 && mod(cellsPerSupercell[1],2) != 0 && mod(cellsPerSupercell[2],2) != 0) { return CRidx_rv(aFs,1); }
	FRk aftf = argFtoF(CRidx_idx(aFs));
	CRidx nv;
        /* cellsPerSupercell[0], cellsPerSupercell[0]+cellsPerSupercell[1], or cellsPerSupercell[0]+cellsPerSupercell[2] even */
        if (mod((2*aftf.kx),cellsPerSupercell[0]) == 0) {
                if (mod((2*aftf.ky),cellsPerSupercell[1]) == 0) {
                        if (aftf.kz <= (cellsPerSupercell[2]/2)) { return aFs; }
                        else { nv = FtoargF(aftf.kx, aftf.ky, cellsPerSupercell[2]-aftf.kz); return CRidx_rv(CRidx_idx(nv), -1*CRidx_sgn(nv)); }
		}
                if (mod((2*aftf.kz),cellsPerSupercell[2]) == 0) {
                        if (aftf.ky <= (cellsPerSupercell[1]/2)) { return aFs; }
                        else { nv = FtoargF(aftf.kx, cellsPerSupercell[1]-aftf.ky, aftf.kz); return CRidx_rv(CRidx_idx(nv), -1*CRidx_sgn(nv)); }
                }
		if (aftf.ky <= (cellsPerSupercell[1]/2)) { return aFs; }
		else { nv = FtoargF(aftf.kx, mod((cellsPerSupercell[1]-aftf.ky),cellsPerSupercell[1]), mod((cellsPerSupercell[2]-aftf.kz),cellsPerSupercell[2])); return CRidx_rv(CRidx_idx(nv), -1*CRidx_sgn(nv)); }
	}
        /* cellsPerSupercell[1], cellsPerSupercell[1]+cellsPerSupercell[2], cellsPerSupercell[2] even not an issue since kx only goes up to int(cellsPerSupercell[0]/2) always (so for (x,cellsPerSupercell[1]/2,z) , (-x,cellsPerSupercell[1]/2,-z) is not stored) */
        return CRidx_rv(aFs,1);
}

#define IsargFunique(aFs) ( ( (argFunique(aFs)) == aFs) ? (true) : (false) )

// Select appropriate bispectrum block based on compacted (k1) + full (k2) pointers
#define bsidx(k1,k2) (CRidx_idx((k1))*numCells+(k2))

