/asl/opt/matlab/R2009b/bin/mexopts.sh

### see https://www.mathworks.com/help/matlab/matlab_external/build-fortran-mex-file.html
###     https://www.mathworks.com/help/matlab/matlab_external/upgrading-mex-files-to-use-64-bit-api.html
###     https://www.mathworks.com/help/matlab/matlab_external/additional-steps-to-update-fortran-source-code.html

FC='/sw/bin/g77'
FFLAGS='-fno-common -fallow-argument-mismatch -w'
FLIBS="-L/usr/lib/ -L/sw/lib -lg2c $MLIBS -lm -compatibleArrayDims"
FOPTIMFLAGS='-O'
FDEBUGFLAGS='-g'

echo 'mexing calcon and CKD 2.5... '
echo $FFLAGS $LDFLAGS $FLIBS

mexF77=/usr/local/matlab/bin/mex
mexF77=/usr/cluster/matlab/r2013a/bin/mex
mexF77=/asl/opt/matlab/R2009b/bin/mex               ## this worked way back, also had  /asl/opt/matlab/R2009b/bin/mexopts.sh
mexF77=/usr/ebuild/software/MATLAB/2021b/bin/mex
mexF77=/usr/ebuild/software/MATLAB/2023b/bin/mex

####$mexF77  convec.F convecg.F                   FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
$mexF77  calcon.F calcong.F                   FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
$mexF77  calconwater.F calconwaterg.F         FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
$mexF77  calconwater_loc.F calconwater_locg.F FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
$mexF77  calconwater_loc_ckd2p5.F calconwater_loc_ckd2p5g.F FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
# 
mv calcon.mexa64                 ../.
mv calconwater.mexa64            ../.
mv calconwater_loc.mexa64        ../.
mv calconwater_loc_ckd2p5.mexa64 ../.
