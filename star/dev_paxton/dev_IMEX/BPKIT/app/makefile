CCC = g++
CC  = gcc
FC  = gfortran # for GCC 4.x
#FC  = g77      # for GCC 3.x

INCDIRS = -I../src
CFLAGS  = $(INCDIRS) -O3
FFLAGS  = -O3
LDFLAGS = 

LALIB  = -llapack -lblas
FC_LIB =

LIBS = ../src/libbpkit.a $(LALIB) $(FC_LIB) -lm -lgfortran # for GCC 4.x
#LIBS = ../src/libbpkit.a $(LALIB) $(FC_LIB) -lm            # for GCC 3.x

# Fortran examples require you to explicitly link Fortran libs on your system
# EXES = bp simple f1example f2example f3example c_example b22_test

EXES = bp simple c_example b22_test

all: $(EXES)

bp: bp.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o bp bp.o $(LIBS)

simple: simple.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o simple simple.o $(LIBS)

f1example: f1example.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o f1example f1example.o $(LIBS)

f2example: f2example.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o f2example f2example.o $(LIBS)

f3example: f3example.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o f3example f3example.o $(LIBS)

c_example: c_example.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o c_example c_example.o $(LIBS)

b22_test: b22_test.o B22.o ../src/libbpkit.a
	$(CCC) $(LDFLAGS) -o b22_test b22_test.o B22.o $(LIBS)

clean:
	rm -f *.o core $(EXES)

.SUFFIXES: .cpp .c .f .o
.cpp.o:
	$(CCC) $(CFLAGS) -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<
