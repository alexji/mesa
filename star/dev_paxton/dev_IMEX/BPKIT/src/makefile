CCC = g++
F77 = gfortran # for GCC 4.x
#F77 = g77      # for GCC 3.x

CFLAGS = -O3 -ansi -pedantic
FFLAGS = -O3

OBJ = $(SGESVD) \
      BpFortran.o \
      BpC.o \
      FGMRES.o \
      HBTMat.o \
      BlockMat.o \
      BRelax.o \
      BILUK.o \
      BTIF.o \
      BlockVec.o \
      CSRMat.o \
      CSRmath.o \
      CSRapinvs.o \
      CSRgmres.o \
      CSRilut.o \
      SparseUtil.o \
      DenseMat.o \
      BpResource.o \
      ilut.o \
      readhb.o \
      spmm.o \
      spsm.o

lib libbpkit.a: $(OBJ)
	ar rcv libbpkit.a $(OBJ)

.SUFFIXES: .cpp .f
.cpp.o:
	$(CCC) $(CFLAGS) -c $*.cpp
.f.o:
	$(F77) $(FFLAGS) -c $*.f

clean:
	rm -f *.o libbpkit.a
