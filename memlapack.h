// BEGIN SAFE MALLOC

// To allow for accurate memory keeping, the size of an allocated memory block is stored in the sizeof(UINT64)
// bytes *before* the pointer that we return, so we can get at it on free.

void *safe_malloc(UINT64 siz) {
	void *rv = NULL;
	#pragma omp critical
	{
	 if (memAllocated+siz+(UINT64)sizeof(UINT64) > MAX_ALLOC_BYTES)
		printErrorAndExit(5, "Maximum allowed memory of %li exceeded trying to allocate a %li byte array\n", 
			MAX_ALLOC_BYTES, siz+(UINT64)sizeof(UINT64));
	 memAllocated += siz+(UINT64)sizeof(UINT64);
	 if (memAllocated > maxMemAllocated) maxMemAllocated = memAllocated;
	} // end omp critical
	rv = calloc(siz+sizeof(UINT64),1);
	if (!rv) printErrorAndExit(5, "Out of memory trying to allocate a %li byte array\n", siz+(UINT64)sizeof(UINT64));
        ((UINT64 *)rv)[0] = siz+(UINT64)sizeof(UINT64);
	return (void *)&(((UINT64 *)rv)[1]);
}

void safe_free(void *mem) {
	UINT64 siz;
	if (mem != NULL) {
		#pragma omp critical
		{
		 siz = ((UINT64 *)mem)[-1];
		 if (siz > memAllocated)
			printWarning("Freeing a %li byte block of memory, but only %li bytes were thought to be allocated!\n",
				siz, memAllocated);
		 memAllocated -= siz;
		} // end omp critical
		free(&(((UINT64 *)mem)[-1]));
	}
	mem = NULL;
}

// END SAFE MALLOC

// BEGIN Lapack
// This code assumes the BLAS and LAPACK library use 64-bit integers when compiled on
// 64-bit systems (so-called "ILP64", NOT "LP64"). Be sure to compile and link against
// ILP64 BLAS and LAPACK, NOT LP64.
// This requirement stems from the fact that even though lp_malloc, etc, have a maximum
// single dimension of 2^31-1, the two dimensional array can be as big as 2^62, which
// requires 64-bit numbers to index when being worked on by BLAS Level 1 functions.
#include "mkl.h"
#ifdef HAVE_LAPACKE_CONFIG_H
#include "lapacke_alt/lapacke/include/lapacke.h"
#endif

LAPACK *lp_malloc(int m, int n, bool cmplx) {
	LAPACK *rv = NULL;
	rv = (LAPACK *)safe_malloc((UINT64)sizeof(LAPACK));
	rv->stride = n;
	rv->rows = m;
	if (cmplx == true) rv->num_per_entry = 2;
	else rv->num_per_entry = 1;

	#pragma omp critical
	{
         if (memAllocated+(UINT64)m*(UINT64)n*(UINT64)sizeof(FTYPE)*(UINT64)rv->num_per_entry > MAX_ALLOC_BYTES)
                printErrorAndExit(5, "Maximum allowed memory of %li exceeded trying to allocate a %ix%i array (%li bytes)\n",
                        MAX_ALLOC_BYTES, m, n, (UINT64)m*(UINT64)n*(UINT64)sizeof(FTYPE)*(UINT64)rv->num_per_entry);
	 memAllocated += (UINT64)m*(UINT64)n*(UINT64)sizeof(FTYPE)*(UINT64)rv->num_per_entry;
	 if (memAllocated > maxMemAllocated) maxMemAllocated = memAllocated;
	} // end omp critical

	rv->arr = (FTYPE *)mkl_malloc((UINT64)m*(UINT64)n*(UINT64)sizeof(FTYPE)*(UINT64)(rv->num_per_entry), 64);
	if (!(rv->arr)) printErrorAndExit(5, "Out of memory trying to allocate a %ix%i array (%li bytes)\n", m, n, (UINT64)m*(UINT64)n*(UINT64)sizeof(FTYPE)*(UINT64)rv->num_per_entry);

	// Warn about very large arrays
	if ((UINT64)m*(UINT64)n*(UINT64)(rv->num_per_entry) > (UINT64)2147483647)
	    printWarning("Warning: a very large array with %li elements, more than 2^31-1, has been allocated. The entire code base has not yet been validated for use with such large arrays, so treat output with care.\n", (UINT64)m*(UINT64)n*(UINT64)(rv->num_per_entry));

	// Zero array
	cblas_dscal((UINT64)m*(UINT64)n*(UINT64)(rv->num_per_entry), 0.0, rv->arr, 1);

	return rv;
}

void lp_free(LAPACK *array) {
	if (array != NULL) {
	 if (array->arr) {
		mkl_free(array->arr);
		#pragma omp critical
		{
                 if ((UINT64)array->stride*(UINT64)array->rows*(UINT64)sizeof(FTYPE)*(UINT64)array->num_per_entry > memAllocated)
                        printWarning("Freeing a %ix%i array (%li bytes), but only %li bytes were thought to be allocated!\n",
                                array->rows, array->stride, (UINT64)array->stride*(UINT64)array->rows*(UINT64)sizeof(FTYPE)*(UINT64)array->num_per_entry, memAllocated);
		 memAllocated -= (UINT64)array->stride*(UINT64)array->rows*(UINT64)sizeof(FTYPE)*(UINT64)array->num_per_entry;
		} // end omp critical
	 }
	 safe_free((void *)array);
	}
}

#define LP_REAL 0
#define LP_IMAG 1
#define lp_idx(array,m,n,cmplx) ((array)->arr[(array)->num_per_entry*((array)->stride*(m)+(n))+cmplx])
#define lp_mag2(array,m,n) (lp_idx(array,m,n,LP_REAL)*lp_idx(array,m,n,LP_REAL)+lp_idx(array,m,n,LP_IMAG)*lp_idx(array,m,n,LP_IMAG))
#define lp_phase(array,m,n) (atan2(lp_idx(array,m,n,LP_IMAG),lp_idx(array,m,n,LP_REAL)))
#define lp_rows(array) ((array)->rows)
#define lp_cols(array) ((array)->stride)
#define lp_ptr(array) (&(lp_idx((array),0,0,LP_REAL)))
#define lp_ptr_row(array,row) (&(lp_idx((array),(row),0,LP_REAL)))

#ifndef lapack_int
#define lapack_int int
#endif

// END Lapack

void freeOrientationGlobals() {
	int i;
        if (orientMatrix != NULL) {
                for (i = 0; i < numOrientationMatrix; i++) lp_free(orientMatrix[i]);
                safe_free(orientMatrix);
        }
        orientMatrix = NULL;
        if (orientShiftMatrix != NULL) {
                for (i = 0; i < numOrientationMatrix; i++) lp_free(orientShiftMatrix[i]);
                safe_free(orientShiftMatrix);
        }
        orientShiftMatrix = NULL;
	if (orientToIdentity != NULL) safe_free(orientToIdentity);
	orientToIdentity = NULL;
	if (orientations != NULL) {
		for (i = 0; i < numOrientations; i++) safe_free(orientations[i]);
		safe_free(orientations);
	}
	orientations = NULL;
}

void freeGlobals() {
	if (atomMatrix != NULL) lp_free(atomMatrix);
	atomMatrix = NULL;
	if (atomMatrixLabels != NULL) safe_free(atomMatrixLabels);
	atomMatrixLabels = NULL;
        if (irrepCoeffs != NULL) lp_free(irrepCoeffs);
        irrepCoeffs = NULL;
	if (irrepOffsets != NULL) lp_free(irrepOffsets);
	irrepOffsets = NULL;
        if (irrep != NULL) lp_free(irrep);
        irrep = NULL;
        if (idealPosition != NULL) lp_free(idealPosition);
        idealPosition = NULL;
	if (LtoC != NULL) lp_free(LtoC);
	LtoC = NULL;
	if (CtoL != NULL) lp_free(CtoL);
	CtoL = NULL;
	if (LtoCnorm != NULL) lp_free(LtoCnorm);
	LtoCnorm = NULL;
        if (ssLtoC != NULL) lp_free(ssLtoC);
        ssLtoC = NULL;
        if (ssCtoL != NULL) lp_free(ssCtoL);
        ssCtoL = NULL;
	LtoCandCtoLprinted = false;

	freeOrientationGlobals();
}


/* Invert matrix */
void lp_invert(LAPACK *m, LAPACK **minv) {
	int i;
	lapack_int *ipiv = NULL;
	LAPACK *mi = *(minv);
	if (m != mi) 
        	cblas_dcopy(lp_rows(m)*lp_cols(m),lp_ptr(m),1,lp_ptr(mi),1);

	ipiv = safe_malloc(sizeof(lapack_int)*MAX(1,MIN(lp_rows(mi),lp_cols(mi))));
        i = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, lp_rows(mi), lp_cols(mi), lp_ptr(mi), lp_cols(mi), ipiv);
	if (i != 0) { safe_free(ipiv); printErrorAndExit(8, "Error in lp_invert on LU factorization: %i\n", i); }
	i = LAPACKE_dgetri(LAPACK_ROW_MAJOR, lp_rows(mi), lp_ptr(mi), lp_cols(mi), ipiv);
	if (i != 0) { safe_free(ipiv); printErrorAndExit(8, "Error in lp_invert on inversion: %i\n", i); }
	safe_free(ipiv);

	return;
}

/* Matrix determinant */
FTYPE lp_det(LAPACK *m) {
	int i;
	FTYPE det = 1.0;
        lapack_int *ipiv = NULL;

        ipiv = safe_malloc(sizeof(lapack_int)*MAX(1,MIN(lp_rows(m),lp_cols(m))));
        i = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, lp_rows(m), lp_cols(m), lp_ptr(m), lp_cols(m), ipiv);
        if (i != 0) { safe_free(ipiv); printErrorAndExit(8, "Error in lp_det on LU factorization: %i\n", i); }

	for (i = 0; i < lp_rows(m); i++) {
	    det *= lp_idx(m,i,i,LP_REAL);
	    if (ipiv[i] != (lapack_int)(i+1)) det *= -1.0;
	}
        safe_free(ipiv);

	return det;
}
