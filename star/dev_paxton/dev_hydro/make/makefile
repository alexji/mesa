
MESA_DIR = ../../../..

STAR = star

# work_standard_makefile

# change this as necessary.  see utils/makefile_header for definitions.

include $(MESA_DIR)/utils/makefile_header

LOAD_MESA = $(LOAD_EXTRAS1) $(LOAD_MESA_STAR) $(LOAD_EXTRAS)

#################################################################

ifndef STAR
STAR = star
endif

#################################################################

INITIAL  = hse3
# full set: bill chevalier dust ejecta hse hse2 hse3 isentropic messy readtable shocktube uniform 

COBJS = dispatch.o profiler.o report.o readpar.o domain.o gridsetup.o geometry.o misc.o \
timestep.o onestep.o riemann.o boundary.o gravity.o plm.o $(INITIAL).o ascii.o euler.o

# STAR_WORK_OBJS defined in makefile_header
OBJS = $(COBJS) $(STAR_WORK_OBJS) run.o

Clib = -lm 

#################################################################

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_JOB_DIR = $(MESA_DIR)/star/job

$(STAR) : $(OBJS)
	$(LOADER) $(FCopenmp) $(Clib) -o $(WORK_DIR)/$(STAR) $(OBJS) $(LOAD_MESA)

#################################################################

# WORK_COMPILE defined in makefile_header

ifeq ($(SKIP_NAN_TRAPS),YES)
run.o: $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $<
else
run.o: $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $(FCtrapNANs) $<
endif

ifeq ($(MESA_CONTRIB_DIR),)
run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f90
	$(WORK_COMPILE) $<
else
run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f90
	$(WORK_COMPILE) -I $(MESA_CONTRIB_DIR)/hooks $<
endif

%.o: $(WORK_SRC_DIR)/%.f90
	$(WORK_COMPILE) $<

%.o: $(STAR_JOB_DIR)/%.f90
	$(WORK_COMPILE) $<

ifeq ($(USE_OPENMP),YES)
CCopenmp = -fopenmp
else
CCopenmp = 
endif
CC = gcc $(CCopenmp) -O3 -g -fno-common -fexceptions
CCchecks = -fbounds-check -Wuninitialized -Warray-bounds
CCwarn = 
WORK_CC_COMPILE = $(CC) $(CCchecks) $(CCwarn) -I$(WORK_SRC_DIR) -c $

%.o: $(WORK_SRC_DIR)/%.c
	$(WORK_CC_COMPILE) $<

clean:
	-@rm -f *.o *.mod *.smod $(WORK_DIR)/$(STAR)








