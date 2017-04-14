// memlapack.h
UINT64 memAllocated = 0; // in bytes
UINT64 maxMemAllocated = 0; // in bytes

// io.h
char curFileLine[MAX_OPEN_FILES][MAX_LINE_LEN] = { { 0 } };
FILE *curFile[MAX_OPEN_FILES] = { NULL };

// xyzinput.h
int numSites = 0;
LAPACK *idealPosition = NULL; // MAX_IDEAL_SITES x 3 (=x,y,z) [only first numSites rows used, real]
char idealPositionLabel[MAX_IDEAL_SITES][MAX_IDEAL_SITES_LABEL_LEN] = { { 0 } };

// irrepinput.h
int numIrreps = 0;
LAPACK *irrep = NULL; // MAX_IRREPS x MAX_IDEAL_SITES*3 (=x,y,z) [only first numSites*3 columns of first numIrreps rows used, real]
char irrepLabel[MAX_IRREPS][MAX_IRREPS_LABEL_LEN] = { { 0 } };
int kvecDegen[3] = { -1, -1, -1 };
int numIrrepByDegen = -1;
int irrepByDegen[MAX_IRREPS][MAX_IRREPS+1] = { { -1 } }; // [x][0] is number in that degeneracy set
int irrepDegen[MAX_IRREPS] = { -1 };
int irrepXYZ[MAX_IRREPS] = { -1 };
int irrepOrder[MAX_IRREPS_ORDERS][MAX_IRREPS] = { { -1 } };
int numOrderSets = 0, numOrders = 0;
int signMap[NUM_SIGN_MAP][3] = {
        {  1,  1,  1 }};// (+,+,+) [ONLY ENTRY MUST BE THIS]

// orientation.h
int numOrientationMatrix = 0;
LAPACK **orientMatrix = NULL;      // Pointers to numOrientationMatrix arrays, each of size 3x3 (real)
LAPACK **orientShiftMatrix = NULL; // Pointers to numOrientationMatrix arrays, each of size 1x3 (real)
int *orientToIdentity = NULL;      // Pointer to numOrientationMatrix integers, each specifying how many times the corresponding
				   // orientMatrix must be applied to return to the identity operator.
int numOrientations = 0;
int **orientations;                // Pointers to numOrientations arrays, each of size numOrientationMatrix, that specify how many
				   // times to apply each orientMatrix to generate that orientation.

// eyefile.h
LAPACK *atomMatrix = NULL;	   // in subcell relative lattice coordinates
int *atomMatrixLabels = NULL;

// supercell (ircontrib.h, eyefile.h, and supercell.h)
int cellsPerSupercell[3] = { 0 }; // Nx, Ny, Nz
int numCells = 0; // Nfull = Nx*Ny*Nz
int numAtoms = 0; // numSites*Nx*Ny*Nz
int numUnique = 0;
FTYPE supercellEdge[3] = { 0.0 };
FTYPE cellEdge[3] = { 0.0 };
FTYPE cellAngles[3] = { 0.0 };
LAPACK *LtoC = NULL, *CtoL = NULL; // 3x3, real, convert between cartesian, angstrom (C) and relative, lattice (L) coordinates (subcell!)
LAPACK *LtoCnorm = NULL; // same as LtoC, but normalized so det(LtoCnorm) = 1
LAPACK *ssLtoC = NULL, *ssCtoL = NULL; // 3x3, real, convert between cartesian, angstrom (C) and relative, lattice (L) coordinates (supercell!)
bool LtoCandCtoLprinted = false;

// ircontrib.h
LAPACK *irrepCoeffs = NULL; // numCells x numIrreps, complex
bool irrepOffset = false; // global offsets removed?
bool irrepOffsetPerSite = false; // where the global offsets on a per-site basis?
LAPACK *irrepOffsets = NULL; // 1 x 3 (global offset), numSites x 3 (site offset), NULL (no offsets), real
