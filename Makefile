CC = icc
CFLAGS = -Wall -O3 -xHost -ip -unroll -no-prec-div -opt-prefetch -scalar_rep -parallel -par-threshold80 -openmp -I/opt/intel/mkl/include -DMKL_ILP64
LIB_LPK = -lmkl_lapack95_ilp64 -Wl,--start-group -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread
# If your version of MKL does not include LAPACKE calls, you must run "make lapacke" in the lapacke_alt/lapacke directory and use these lines instead
#CFLAGS = -Wall -O3 -xHost -ip -unroll -no-prec-div -opt-prefetch -scalar_rep -parallel -par-threshold80 -openmp -I/opt/intel/mkl/include -DMKL_ILP64 -DHAVE_LAPACKE_CONFIG_H -DUSE_MKL_LAPACK
#LIB_LPK = ./lapacke_alt/liblapacke.a -lmkl_lapack95_ilp64 -Wl,--start-group -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread


all: rmc-irrep-decompose rmc-irrep-rebuild rmc-breakup rmc-symmetry
headers: breakup.h eyefile.h io.h irrepdecompose.h irreprebuild.h orientation.h xyzinput.h defs.h globals.h ircontrib.h irrepinput.h memlapack.h supercell.h

rmc-irrep-decompose : rmc-irrep-decompose.c headers
	$(CC) -I. -I$(srcdir) $(CFLAGS) rmc-irrep-decompose.c -o rmc-irrep-decompose $(LIB_LPK)
rmc-irrep-rebuild : rmc-irrep-rebuild.c headers
	$(CC) -I. -I$(srcdir) $(CFLAGS) rmc-irrep-rebuild.c -o rmc-irrep-rebuild $(LIB_LPK)
rmc-breakup: rmc-breakup.c headers
	$(CC) -I. -I$(srcdir) $(CFLAGS) rmc-breakup.c -o rmc-breakup $(LIB_LPK)
rmc-symmetry: rmc-symmetry.c headers
	$(CC) -I. -I$(srcdir) $(CFLAGS) rmc-symmetry.c -o rmc-symmetry $(LIB_LPK)

clean:
	rm -f rmc-irrep-decompose rmc-irrep-rebuild rmc-breakup rmc-symmetry rmc-irrep-decompose.o rmc-irrep-rebuild.o rmc-breakup.o rmc-symmetry.o
