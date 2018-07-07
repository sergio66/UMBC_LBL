#
## ------------------------ Linux Absoft f77 --------------------------
# Linux with Absoft Fortran ... does not support fseek the way we want!
# and so it SHOULD NOT be used to compile the F77 readers. However it 
# supports structures, and so it can be used to compile makeRTPfile.x
# Absoft compiler options
# -f    fold all names to lower case
# -N109 fold all names to upper case
# -W    wide source file
# -w    suppress warning messages (absoft is more fussy than SGI or g77)
# -A    alignment warnings
# -C    check array bounds
# -O    some optimizations
# -N3   add record info to unformatted files
# -s    static allocation
# -N2   force intrinsic double functions
# -N113 force double precision 
# -N114 force untyped variables as warnings
# -N15  so that LAPACK lib can be included; adds trailing underscores 
# -f    fold all names to lower case 
F77   = f77
F77=/asl/opt/absoft/absoft10.0/bin/af77
FLAGS = -w -W -s -f -C -O2 -N114 -N15 -N3 -B108
LIBS = -lU77 -llapack -lblas
LIBS = -lU77 -L/asl/opt/absoft/absoft10.0/extras/lapack/precompiled/64bit/ -llapack -lblas

PARAMS = parameters.inc
OBJS   = LM_calc_CO2_new-coorected.o
OBJS2   = aux_LM_calc_CO2_new-coorected.o loop_code_umbclbl_fast.o
OBJS3   = aux_LM_calc_CO2_new-coorected.o Total_program_win.o

EXECS  = LM_calc_CO2_new-coorected.x loop_code_umbclbl_fast.x Total_program_win.x
all : $(EXECS)

basic: LM_calc_CO2_new-coorected.x loop_code_umbclbl_fast.x Total_program_win.x

############# these are the main execs   ###################
LM_calc_CO2_new-coorected.x: $(OBJS) $(PARAMS)
	$(F77) -o LM_calc_CO2_new-coorected.x $(FLAGS) $(OBJS) $(LIBS)

loop_code_umbclbl_fast.x: $(OBJS2) $(PARAMS)
	$(F77) -o loop_code_umbclbl_fast.x $(FLAGS) $(OBJS2) $(LIBS)

Total_program_win.x: $(OBJS3) $(PARAMS)
	$(F77) -o Total_program_win.x $(FLAGS) $(OBJS3) $(LIBS)


############# these are the main execs   ###################

LM_calc_CO2_new-coorected.o : LM_calc_CO2_new-coorected.f $(PARAMS)
	$(F77) -c $(FLAGS) LM_calc_CO2_new-coorected.f

aux_LM_calc_CO2_new-coorected.o : aux_LM_calc_CO2_new-coorected.f $(PARAMS)
	$(F77) -c $(FLAGS) aux_LM_calc_CO2_new-coorected.f

loop_code_umbclbl_fast.o : loop_code_umbclbl_fast.f $(PARAMS)
	$(F77) -c $(FLAGS) loop_code_umbclbl_fast.f

Total_program_win.o : Total_program_win.f $(PARAMS)
	$(F77) -c $(FLAGS) Total_program_win.f

clean:  
	rm *.o *~ 2> /dev/null || true 

