.SUFFIXES: .F90 .o
TARGET = ../sis

#------ gfortran
FC = gfortran
CPP = -cpp
#OPT = -O3 -Wall -ffree-line-length-none -fbounds-check -fno-range-check -fallow-invalid-boz -frecursive
OPT = -fopenmp -O3 -Wall -ffree-line-length-none -fbounds-check -fno-range-check -fallow-invalid-boz -frecursive

# ------------------------------------------
#  Git:
GIT_VER = $(shell git log -1 --format="%h")
GIT_BRA = $(shell git status | head -1| cut -c11-)
DATESEQ = $(shell date -u +%Y%m%d%H%M)  # Date
VERSION = -DVERGIT -DVERBUILD="\"$(GIT_VER)\"" -DVERBRANCH="\"$(GIT_BRA)\"" -DVERDATE="\"$(DATESEQ)\""

# ------------------------------------------
FLAG = $(CPP) $(OPT)

# ------------------------------------------

# Constant modules 
FMOD1 = const.o const_idx.o const_phys.o

# Variable modules
FMOD2 = var_top.o var_state.o var_potential.o var_io.o

# File modules 
FMOD3 = dcd.o

# subroutines
FOBJ1 = \
      list_local.o list_bp.o \
	energy.o\
	energy_bond.o energy_angl.o energy_bp.o\
	read_sisinfo.o read_force_field.o

OBJS = $(FMOD1) $(FMOD2) $(FMOD3) $(FOBJ1)

all : $(TARGET)

$(TARGET): $(OBJS) main.o
	$(FC) -o $@ $(FLAG) $(OBJS) main.o

clean:
	/bin/rm -f *.o *.mod $(TARGET)

main.o:
	$(FC) $(VERSION) $(FLAG) -c main.F90

# suffix rule
.F90.o:
	$(FC) $(FLAG) -c $< 
