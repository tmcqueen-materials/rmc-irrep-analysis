// HARD LIMITS (fixed by encoding, CPU, or compiler limits; code will likely break if these changed without modifications elsewhere)
// These are appropriate values for 64-bit systems
// MAX_ALLOC_BYTES can be increased safely, if you have the RAM
#define MAX_SUPERCELLS_EDGE 2147483647  // practically limited to less by next constraint in most cases
#define MAX_WAVEVECTORS     2147483647
#define MAX_ALLOC_BYTES     5368709120  // 5 GB, but can be increased (total RAM should be at least this number plus 500 MB)
// For 32-bit compiles, must use these values instead!
// #define MAX_SUPERCELLS_EDGE 32767
// #define MAX_WAVEVECTORS     32767
// #define MAX_ALLOC_BYTES     2147483647 // Hard limit on 32-bit systems

// COMPILED-IN DEFAULTS (could change to read from cfg file if desired)
#define MAX_IDEAL_SITES 256
#define MAX_IDEAL_SITES_LABEL_LEN 7
#define MAX_IDEAL_ACTUAL_DISTANCE 1.0 // in angstrom
#define MAX_SYSTEMATIC_SHIFT 0.1 // in angstrom
#define MAX_IRREPS 768 // should be MAX_IDEAL_SITES*3
#define MAX_IRREPS_LABEL_LEN 7
#define MAX_IRREPS_ORDERS 192
#define ZERO_TOL 1e-11
#define MAX_OPEN_FILES 2
#define MAX_LINE_LEN 65536 //65535+1 for terminating character
#define MAX_ORIENTATIONS_TO_IDENTITY 6
#define MAX_ORIENTATIONS 192 // This is more than will be present in any (non-magnetic) space group
#define NUM_SIGN_MAP 1

// Constants
#define PI 3.14159265358979323846
#define FTYPE double

// Select principle value of an angle (in radians)
#define pri2Pi(x) (((x) < 0.0) ? ( (2.0*PI)*((int)(((x)-PI)/(2.0*PI+ZERO_TOL))) ) : ( (2.0*PI)*((int)(((x)+PI)/(2.0*PI-ZERO_TOL))) ) )
#define mod2Pi(x) ( (x) - pri2Pi((x)) )
// only to pi variant
#define modPi(x) (((x) < 0.0) ? ((x) - ((PI)*((int)(((x)-PI/2.0)/(PI+ZERO_TOL))))  ) : ( (x) - ((PI)*((int)(((x)+PI/2.0)/(PI-ZERO_TOL)))) ) )

// Convert to range [0.0,1.0)
#define mod1(x) ( (x) - (FTYPE)pri1((x)) )
#define pri1(x) ( ((x) < 0.0) ? ( ((int)((x-1.0)/1.0)) ) : ( ((int)((x)/1.0)) ) )
// Convert to range [-0.5, 0.5)
#define modp5(x) ( mod1((x)+0.5) - 0.5 )
#define prip5(x) ( pri1((x)+0.5) )

// Standard C defines the modulus as having the same sign as the dividend
// but we need the modulus of a negative number to always the same sign as divisor.
// This macro implements this version of modulus (n must be positive),
// assuming the compiler implements it as same sign as the dividend OR
// same sign as the divisor.
#define mod(a,n) ( ((a)%(n) < 0) ? ((a)%(n) + n) : ((a)%(n)) )

// Define a function returning the sign of an integer (use fsgn for FTYPEs)
#define sgn(rv) ( ((rv) >= 0) - ((rv) < 0) )

// Define boolean if needed
#ifndef bool
#define bool int
#define true 1
#define false 0
#endif

// Define UINT64 if needed
#define UINT64 uint64_t

// Define min and max if needed
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif

// Define a complex value equal to 1
FTYPE CMPLX1p0i[2] = { 1.0, 0.0 };
// Define a complex value equal to 0
FTYPE CMPLX0p0i[2] = { 0.0, 0.0 };

// Define a matrix type
typedef struct {
        int stride;
        int rows;
        int num_per_entry; // 1 = real, 2 = complex
        FTYPE *arr;
} LAPACK;

/* A kx,ky,kz tuple is stored as a structure, as are pairs of kx,ky,kz tuples.
 */
typedef struct {
        int kx;
        int ky;
        int kz;
} FRk;

/* Two kx,ky,kz tuples are needed for higher order spectral analysis
 */
typedef struct {
        FRk k1;
        FRk k2;
} FRk2;

/* An index and phase offset is needed for centering/order orientation
 */
typedef struct {
	int ir;       // irrep
	FTYPE ph;     // constant (k independent) phase offset
	int tph[3];   // translational phase offsets (in 2*PI/cellsPerSupercell units)
	int ksign[3]; // sign (if any) to apply to kx,ky,kz components
} COO;
