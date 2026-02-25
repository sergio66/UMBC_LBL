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
echo "########################################################################"

# echo "########################################################################"
# 
# ### the : <<'END'  and  END are commenting this out for now since you can separately call calconwater_locg_ckd2p5.sc
# ### the : <<'END'  and  END are commenting this out for now since you can separately call calconwater_locg_ckd2p5.sc
# ### the : <<'END'  and  END are commenting this out for now since you can separately call calconwater_locg_ckd2p5.sc
# 
# : <<'END'
# 
# ## make CKD2.5
# ## ------------
# echo "making CKD2.5"
# echo " "
# #$mexF77  convec.F convecg.F                                 FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
# $mexF77  calcon.F calcong.F                                  FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
# $mexF77  calconwater.F calconwaterg.F                        FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
# $mexF77  calconwater_loc.F calconwater_locg.F                FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
# $mexF77  calconwater_loc_ckd2p5.F calconwater_loc_ckd2p5g.F  FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
# 
# END
# 
# echo "########################################################################"

FC='gfortran90'
FFLAGS='-fexceptions -fbackslash'
FFLAGS="$FFLAGS -fno-underscoring -fPIC -fno-omit-frame-pointer "
FLIBS="$RPATH $MLIBS -lm"
FOPTIMFLAGS='-O'
FDEBUGFLAGS='-g'

## make CKD3.2
## ------------
echo "making CKD3.2"
echo " "
rm calcon_loc_mtckd_32_wrap.o ../calconwater_loc_ckd3p2.mexa64

### https://www.mathworks.com/matlabcentral/answers/334280-how-should-i-mex-a-fortran-code-main-function-that-calls-other-fortran-codes-sub-function-that-in
### $mexF77  calconwater_loc_ckd3p2.F calconwater_loc_ckd3p2g.F    FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'
ln -s MT_CKD3.2/cntnm/build/lblparams.mod     lblparams.mod
ln -s MT_CKD3.2/cntnm/build/phys_consts.mod   phys_consts.mod
ln -s MT_CKD3.2/cntnm/build/planet_consts.mod planet_consts.mod
#### calconwater_loc_ckd3p2.F90 is the gatewy function to calconwater_loc_ckd3p2.F90
#$mexF77 calcon_loc_mtckd_32_wrap.F90                               FFLAGS='$FFLAGS' 
$mexF77 -v -c calcon_loc_mtckd_32_wrap.F90                         FFLAGS='$FFLAGS' 
$mexF77 calconwater_loc_ckd3p2.F90 calcon_loc_mtckd_32_wrap.o      FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'

echo "########################################################################"

# mv calcon.mexa64                 ../.
# mv calconwater.mexa64            ../.
# mv calconwater_loc.mexa64        ../.
# mv calconwater_loc_ckd2p5.mexa64 ../.
mv calconwater_loc_ckd3p2.mexa64 ../.

ls -lt calcon_loc_mtckd_32_wrap.o ../calconwater_loc_ckd3p2.mexa64
