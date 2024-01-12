# 
# makefile to build
#   - hselect
#   - read_hitran.mex

HDIR = SAMPLE_HITREAD
#MBIN = /asl/opt/matlab/R2009b/bin
MBIN = /usr/ebuild/software/MATLAB/2021b/bin/

all: hselect mex

hselect: hselect.o hutils.o
	cc -o hselect hselect.o hutils.o

hselect.o: hselect.c hdefs.h
	cc -c hselect.c

hutils.o: hutils.c hdefs.h
	cc -c hutils.c

mex: read_hitran.c hutils.c hdefs.h
	$(MBIN)/mex -compatibleArrayDims read_hitran.c hutils.c 

clean:
	rm *.o hselect 2> /dev/null || true

dist:
	 cd .. ;				\
	 tar -zcf $(HDIR)/$(HDIR).tgz		\
	 $(HDIR)/README				\
	 $(HDIR)/read_hitr.txt			\
	 $(HDIR)/Makefile			\
	 $(HDIR)/hdefs.h			\
	 $(HDIR)/hutils.c			\
	 $(HDIR)/hselect.c			\
	 $(HDIR)/read_hitran.c			\
	 $(HDIR)/read_hitran2.m                 \
	 $(HDIR)/hitread.m                      \
	 $(HDIR)/translate2oldHITparams.m       \
	 $(HDIR)/driver_findlines_plot.m        \
	 $(HDIR)/g5.dat              

