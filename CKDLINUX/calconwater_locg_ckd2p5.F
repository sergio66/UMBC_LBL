/asl/opt/matlab/R2009b/bin/mexopts.sh

            FC='/sw/bin/g77'
            FFLAGS='-fno-common'
            FLIBS="-L/usr/lib/ -L/sw/lib -lg2c $MLIBS -lm -compatibleArrayDims"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'

echo 'mexing CKD... '
echo $FFLAGS $LDFLAGS $FLIBS

#/asl/opt/matlab/R2009b/bin/mex  convec.F convecg.F                   FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
/asl/opt/matlab/R2009b/bin/mex  calcon.F calcong.F                   FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
/asl/opt/matlab/R2009b/bin/mex  calconwater.F calconwaterg.F         FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
/asl/opt/matlab/R2009b/bin/mex  calconwater_loc.F calconwater_locg.F FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
/asl/opt/matlab/R2009b/bin/mex  calconwater_loc_ckd2p5.F calconwater_loc_ckd2p5g.F FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'

mv calcon.mexa64                 ../.
mv calconwater.mexa64            ../.
mv calconwater_loc.mexa64        ../.
mv calconwater_loc_ckd2p5.mexa64 ../.
