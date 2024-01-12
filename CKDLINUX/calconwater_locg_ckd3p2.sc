/usr/cluster/matlab/2016b/bin/mexopts.sh

#            FC='/sw/bin/g77'
#            FFLAGS='-fno-common'
#            FLIBS="-L/usr/lib/ -L/sw/lib -lg2c $MLIBS -lm -compatibleArrayDims"
#            FOPTIMFLAGS='-O'
#            FDEBUGFLAGS='-g'

echo 'mexing CKD... '
echo $FFLAGS $LDFLAGS $FLIBS

mexF77=/usr/local/matlab/bin/mex
mexF77=/usr/cluster/matlab/r2013a/bin/mex
mexF77=/usr/cluster/matlab/r2016b/bin/mex  ## this worked way back
mexF77=/asl/opt/matlab/R2009b/bin/mex
mexF77=/usr/ebuild/software/MATLAB/2023b/bin/mex
mexF77=/usr/ebuild/software/MATLAB/2021b/bin/mex

#$mexF77  convec.F convecg.F                   FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
$mexF77  calcon.F calcong.F                   FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
$mexF77  calconwater.F calconwaterg.F         FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
$mexF77  calconwater_loc.F calconwater_locg.F FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
$mexF77  calconwater_loc_ckd2p5.F calconwater_loc_ckd2p5g.F FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'

FC='gfortran90'
FFLAGS='-fexceptions -fbackslash'
FFLAGS="$FFLAGS -fno-underscoring -fPIC -fno-omit-frame-pointer "
FLIBS="$RPATH $MLIBS -lm"
FOPTIMFLAGS='-O'
FDEBUGFLAGS='-g'

### https://www.mathworks.com/matlabcentral/answers/334280-how-should-i-mex-a-fortran-code-main-function-that-calls-other-fortran-codes-sub-function-that-in
### $mexF77  calconwater_loc_ckd3p2.F calconwater_loc_ckd3p2g.F FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
ln -s MT_CKD3.2/cntnm/build/lblparams.mod     lblparams.mod
ln -s MT_CKD3.2/cntnm/build/phys_consts.mod   phys_consts.mod
ln -s MT_CKD3.2/cntnm/build/planet_consts.mod planet_consts.mod
$mexF77 -v -c calcon_loc_mtckd_32_wrap.F90                          FFLAGS='$FFLAGS' 
$mexF77 calconwater_loc_ckd3p2.F90 calcon_loc_mtckd_32_wrap.o FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'

mv calcon.mexa64                 ../.
mv calconwater.mexa64            ../.
mv calconwater_loc.mexa64        ../.
mv calconwater_loc_ckd2p5.mexa64 ../.
mv calconwater_loc_ckd3p2.mexa64 ../.
