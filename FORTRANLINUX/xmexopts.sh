#
# mexopts.sh	Shell script for configuring MEX-file creation script,
#               mex.  These options were tested with the specified compiler.
#
# usage:        Do not call this file directly; it is sourced by the
#               mex shell script.  Modify only if you don't like the
#               defaults after running mex.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_MEX_OPT: Template Options file for building MEX-files via the system ANSI compiler
#
# Copyright 1984-2000 The MathWorks, Inc.
# $Revision: 1.77 $  $Date: 2000/09/21 20:59:54 $
#----------------------------------------------------------------------------
#
    TMW_ROOT="$MATLAB"
    MFLAGS=''
    MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex"
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        alpha)
#----------------------------------------------------------------------------
#           cc -V | grep UNIX
#           DEC C V5.9-008 on Digital UNIX V4.0 (Rev. 1229)
#           Digital UNIX Compiler Driver 3.11
            CC='cc'
            CFLAGS='-shared -ieee -pthread -std1'
            CLIBS="$MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
#           f77 -what
#           Compaq Fortran 77 Driver V5.3-11
#           Compaq Fortran 77 V5.3-189-449BB
            FC='f77'
            FFLAGS='-shared -fpe3 -pthread'
            FLIBS="$MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-pthread -shared -Wl,-expect_unresolved,'*',-hidden,-exported_symbol,$ENTRYPOINT,-exported_symbol,mexVersion"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        hpux)
#----------------------------------------------------------------------------
#           what `which cc`
#           HP92453-01 A.11.01.00 HP C Compiler
#            $ PATCH/11.00:PHCO_95167  Oct  1 1998 13:46:32 $
            CC='cc'
            CFLAGS='+Z +DA2.0 -D_POSIX_C_SOURCE=199506L -Wp,-H65535 -Ae'
            CLIBS="$MLIBS -lm -lc"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
#           what `which f90`
#           HP-UX f90 20000117 (163844)  B3907DB/B3909DB B.11.01.11
#           HP F90 v2.4
#            $ PATCH/11.00:PHCO_95167  Oct  1 1998 13:46:32 $
            F90LIBDIR=`which f90 | sed -n -e '1s|bin/f90|lib/pa2.0|p'`
            FC='f90'
            FFLAGS='+Z +DA2.0'
            FLIBS="$MLIBS -lm -L$F90LIBDIR -lF90 -lcl -lc -lisamstub"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS="-b +e $ENTRYPOINT +e mexVersion"
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        hp700)
#----------------------------------------------------------------------------
#           what `which cc`
#           HP92453-01 A.10.32.30 HP C Compiler
            CC='cc'
#           Remove +DAportable from CFLAGS if you wish to optimize
#           for target machine
            CFLAGS='+Z -Ae +DAportable -Wp,-H65535'
            CLIBS="$MLIBS -lc"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
#           what `which f90`
#           HP-UX f90 20000107 (183817)  B3907DB/B3909DB B.10.20.19
#           HP F90 v2.4
            F90LIBDIR=`which f90 | sed -n -e '1s|bin/f90|lib|p'`
            FC='f90'
            FFLAGS='+Z +DAportable'
            FLIBS="$MLIBS -L$F90LIBDIR -lF90 -lcl -lc -lisamstub"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS="-b +e $ENTRYPOINT +e mexVersion +e errno"
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        ibm_rs)
#----------------------------------------------------------------------------
#           lslpp -l ibmcxx.cmp
#           3.6.6.0  COMMITTED  IBM C and C++ Compilers
            CC='cc'
            CFLAGS='-D_THREAD_SAFE -D_ALL_SOURCE -qchars=signed -qlanglvl=ansi'
            CLIBS="$MLIBS -lmatlb -lmat -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
#           lslpp -cl xlfcmp
#           7.1.0.0  COMMITTED  I XL Fortran Compiler
            FC='f77'
            FFLAGS=''
            FLIBS="$MLIBS -lmat -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-bE:$TMW_ROOT/extern/lib/$Arch/$MAPFILE -bM:SRE -bnoentry"
            LDOPTIMFLAGS='-O -Wl,-s'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        glnx86)
#----------------------------------------------------------------------------
            RPATH="-Wl,--rpath-link,$TMW_ROOT/extern/lib/$Arch,--rpath-link,$TMW_ROOT/bin/$Arch"
#           gcc -v
#           gcc version 2.95.2 19991024 (release)
            CC='gcc'
            CFLAGS='-fPIC -ansi -D_GNU_SOURCE -pthread'
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
#           g77 -v -xf77-version 
#           g77 version 2.95.2 19991024 (release) 
#           (from FSF-g77 version 0.5.25 19991024 (release))
#           NOTE: g77 is not thread safe
            FC='g77'
            FFLAGS='-fPIC'
            FLIBS="$RPATH $MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-pthread -shared -Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sgi)
#----------------------------------------------------------------------------
#           cc -version
#           MIPSpro Compilers: Version 7.3.1.2m
            CC='cc'
            CFLAGS='-n32 -signed -D_POSIX_C_SOURCE=199506L -D__EXTENSIONS__ -D_XOPEN_SOURCE -D_XOPEN_SOURCE_EXTENDED'
            CLIBS="-dont_warn_unused $MLIBS -lm -lc"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
#           f77 -version
#           MIPSpro Compilers: Version 7.3.1.2m
            FC='f77'
            FFLAGS='-n32'
            FLIBS="-dont_warn_unused $MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-n32 -shared -exported_symbol $ENTRYPOINT -exported_symbol mexVersion"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sol2)
#----------------------------------------------------------------------------
#           cc -V
#           WorkShop Compilers 5.0 98/12/15 C 5.0
            CC='cc'
            CFLAGS='-KPIC -dalign -xlibmieee -Xc -D__EXTENSIONS__ -D_POSIX_C_SOURCE=199506L -mt'
            SUNLIBDIR=`which cc | sed -n -e '1s|bin/cc|lib|p'`
            CLIBS="$MLIBS -L$SUNLIBDIR -lm -lc"

            COPTIMFLAGS='-xO3 -xlibmil -DNDEBUG'
            CDEBUGFLAGS='-g'
#
#           f77 -V
#           WorkShop Compilers 5.0 99/09/16 FORTRAN 77 5.0 patch 107596-03
            FC='f77'
            FFLAGS='-KPIC -dalign -mt'
            FLIBS="$MLIBS -lF77 -lM77 -lsunmath -lm -lcx -lc"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-G -mt -M$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
