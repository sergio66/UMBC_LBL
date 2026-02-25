#/usr/ebuild/installs/software/MATLAB/2023b/bin/mexopts.sh

#            FC='/sw/bin/g77'
#            FFLAGS='-fno-common'
#            FLIBS="-L/usr/lib/ -L/sw/lib -lg2c $MLIBS -lm -compatibleArrayDims"
#            FOPTIMFLAGS='-O'
#            FDEBUGFLAGS='-g'

echo 'mexing CKD... '
echo $FFLAGS $LDFLAGS $FLIBS

## set the mex compiler
## --------------------
mexF77=/usr/local/matlab/bin/mex
mexF77=/usr/cluster/matlab/r2013a/bin/mex
mexF77=/usr/cluster/matlab/r2016b/bin/mex  ## this worked way back
mexF77=/asl/opt/matlab/R2009b/bin/mex
mexF77=/usr/ebuild/software/MATLAB/2023b/bin/mex
mexF77=/usr/ebuild/software/MATLAB/2021b/bin/mex
mexF77=/usr/ebuild/installs/software/MATLAB/2023b/bin/mex
echo "setting mexF77 = " $mexF77
echo " "

echo "########################################################################"
echo "use calconwater_locg_ckd2p5.sc for MT CKD 2.5"
echo "use calconwater_locg_ckd3p2.sc for MT CKD 3.2"
echo "########################################################################"

#        NCL = /usr/lib
#        NCI = /usr/include
#        NCI = /usr/ebuild/installs/software/netCDF-Fortran/4.6.1-iimpi-2023a/include
#        INCLUDES:= -I. -I$(NCI)
#	NETCDF = yes

NCI=/usr/ebuild/installs/software/netCDF-Fortran/4.6.1-iimpi-2023a/include
FINC="-I$NCI"
FC='gfortran90'
FFLAGS='-fexceptions -fbackslash'
FFLAGS="$FFLAGS -fno-underscoring -fPIC -fno-omit-frame-pointer "
FFLAGS="$FINC $FFLAGS -fno-underscoring -fPIC -fno-omit-frame-pointer "
FLIBS="$RPATH $MLIBS -lm"
FOPTIMFLAGS='-O'
FDEBUGFLAGS='-g'

echo $FFLAGS

## make CKD4.3
## ------------
echo "making CKD4.3"
echo " "
rm calcon_loc_mtckd_43_wrap.o ../calconwater_loc_ckd4p3.mexa64

### https://www.mathworks.com/matlabcentral/answers/334280-how-should-i-mex-a-fortran-code-main-function-that-calls-other-fortran-codes-sub-function-that-in
### $mexF77  calconwater_loc_ckd3p2.F calconwater_loc_ckd3p2g.F    FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
#ln -s MT_CKD3.2/cntnm/build/lblparams.mod     lblparams.mod
#ln -s MT_CKD3.2/cntnm/build/phys_consts.mod   phys_consts.mod
#ln -s MT_CKD3.2/cntnm/build/planet_consts.mod planet_consts.mod

ln -s MT_CKD_H2O-4.3/build/mt_ckd_h2o.mod      mt_ckd_h2o.mod
ln -s MT_CKD_H2O-4.3/build/read_file.mod       read_file.mod
# ln -s MT_CKD_H2O-4.3/build/phys_consts.mod     phys_consts.mod
# diff MT_CKD3.2/cntnm/build/phys_consts.mod  MT_CKD_H2O-4.3/build/phys_consts.mod

#### calconwater_loc_ckd4p3.F90 is the gatewy function to calconwater_loc_ckd4p3.F90
#### this is testing we can compile it, but it does not make .o
##   $mexF77 calcon_loc_mtckd_43_wrap.F90                               FFLAGS='$FFLAGS' 
##   even easier : gfortran calcon_loc_mtckd_43_wrap.F90

#### use these two for run8watercontinnum
#$mexF77 -v -c MT_CKD_H2O-4.3/build/ calcon_loc_mtckd_43_wrap.F90                         FFLAGS='$FFLAGS' 
#$mexF77 -I'MT_CKD_H2O-4.3/build/' calconwater_loc_ckd4p3.F90 calcon_loc_mtckd_43_wrap.o      FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
$mexF77 -v -c  calcon_loc_mtckd_43_wrap.F90                         FFLAGS='$FFLAGS' 
$mexF77  calconwater_loc_ckd4p3.F90 calcon_loc_mtckd_43_wrap.o      FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'

echo "########################################################################"

#mv calconwater_loc_ckd2p5.mexa64 ../.
#mv calconwater_loc_ckd3p2.mexa64 ../.
mv calconwater_loc_ckd4p3.mexa64 ../.

ls -lt calcon_loc_mtckd_43_wrap.o ../calconwater_loc_ckd4p3.mexa64

