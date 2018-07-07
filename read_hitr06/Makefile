# 
# makefile to build
#   - hselect
#   - read_hitran.mex

HDIR = read_hitr06
# MBIN = /usr/local/matlab/bin
# MBIN = /usr/local/matlabR2006b/bin
MBIN = /asl/opt/matlab/R2009b/bin

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
	 $(HDIR)/read_hitran2.m

