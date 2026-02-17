# 
# makefile to build
#   - hselect
#   - read_hitran.mex

MBIN = /usr/ebuild/installs/software/MATLAB/2023b/bin/
# MBIN = /usr/ebuild/installs/software/MATLAB/2025b/bin/

all: hselect mex

hselect: hselect.o hutils.o
	cc -o hselect hselect.o hutils.o

hselect.o: hselect.c hdefs.h
	cc -c hselect.c

hutils.o: hutils.c hdefs.h
	cc -c hutils.c

mex: read_hitran.c hutils.c hdefs.h
	$(MBIN)/mex read_hitran.c hutils.c 

clean:
	rm *.o hselect 2> /dev/null || true

