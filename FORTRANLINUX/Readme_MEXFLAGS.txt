in Matlab prompt
>> mex -setup fortran
in Linux prompt
>> xnw ~/.matlab/R2021b/mex_FORTRAN_glnxa64.xml
    add in eg FFLAGS = -fallow-invalid-boz  into the list
    add in /usr/ebuild/software/GCCcore/8.2.0/bin/gfortran
by replacing 10.2.0 with 8.2.0 in 
       LDDEBUGFLAGS=&quot;$LDDEBUGFLAGS&quot;" GFORTRAN_INSTALLDIR="/usr/ebuild/software/GCCcore/8.2.0/bin" GFORTRAN_LIBDIR="/umbc/ebuild-soft/skylake/software/GCCcore/8.2\
.0/lib64" GFORTRAN_VERSION="8.2.0"/>


OR AT UNIX COMMAND LINE do eg
  mex -v loop_chi2.F
to see what mex is doing

can try to reset compiler using eg
 mex -v FC="/usr/ebuild/software/GCCcore/8.2.0/bin/gfortran" loop_chi2.F
but hmmm, dos this really work?????

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

https://jponttuset.cat/matlab2014-mex-flags/

Adding flags to the mex Matlab 2014 compiler

Up until Matlab 2013b, one could modify the C/C++ mex compilation
flags by editing the mexopts.sh file, usually found in
~/.matlab/mexopts.sh. A common example nowadays is adding the
-std=c++11 option to compile C++ code from the C++11 standard (which
has some awesome features). To do that in Matlab 2013 and below, you
would simply add the flag to the CXXFLAGS variable.


In Matlab 2014, however, things got a little bit more complicated than
that. Googling a little I found that now the options were stored in
XML files called something like mex_C_glnxa64.xml or
mex_C++_glnxa64.xml, but these were nowhere to be found in my system,
until I finally realized you have to call mex setup before to create
them.


The final steps are as follow:

Run mex -setup C++ or mex -setup C, depending on which language you
are gonna compile. This will create the corresponding XML file.  Open
the file ~/.matlab/mex_C++_glnxa64.xml or
~/.matlab/mex_C_glnxa64.xml. In some versions of Matlab it’s in
~/.matlab/R2014a/mex_C++_glnxa64.xml.  Edit the CXXFLAGS line to
something like: CXXFLAGS="-ansi -fexceptions -fPIC
-fno-omit-frame-pointer -pthread -std=c++11" And you’re good to go!


You might also need to change the compiler used to a non-default
one. You can do so by calling mex with mex GCC='/path/to/bin/g++'
you_file.cpp
