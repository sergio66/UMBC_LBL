#
# Makefile make_sarta_pclsam
#
# ------------
# Object files
# ------------
OBJS1 = cpf3.o cpf.o loop_pCqSDHC.o loop_pCqSDHC_plus.o pCqSDHC.o q_loop_SDV.o\
       qsdhc.o qSDV.o SDvhh2RI.o driver_HTP.o

OBJS2 = cpf3.o cpf.o loop_pCqSDHC.o loop_pCqSDHC_plus.o pCqSDHC.o q_loop_SDV.o\
       qsdhc.o qSDV.o SDvhh2RI.o driver_HTP_all.o

# ------------------------------------------------------------------------
# Compilers
# ------------------------------------------------------------------------

# INTEL ifort F77
# https://software.intel.com/en-us/node/522723 for profiling
# ------------------------------------------------------------------------
# -u  : turn off implicit typing of variables
# -g  : generate debugging information (turns off optimization)
# -C  : do run time subscript range checking
# -w0 : inform about unused variables
# -O3 : heavy optimization
# -64 : 64-bit objects (libraries must match)
# -us : puts trailing underscore
F77 = ifort
FLIB = -mkl

OPT = -O2 -extend-source 132 -W1

# Portland Group F77:
# ------------------------------------------------------------------------
#    -byteswapio  read big-endian binary files (needed for coefficient files)
#    -C           array bounds checking
#    -g           include symbol table info
#    -O           moderate optimization
#FLIB = 
#OPT = -O -Mextend
#F77 = /cm/shared/apps/pgi/current/linux86-64/16.5/bin/pgf95

# ------------------------------------------------------------------------
# Build
# ------------------------------------------------------------------------
htp:	$(OBJS1)
	$(F77) $(OBJS1) $(FLIB)

htp_all:$(OBJS2)
	$(F77) $(OBJS2) $(FLIB)

# Compile source to object files
.F.o:
	$(F77) -c $(OPT) $*.F

clean:
	rm -f *.o *.x *.out
