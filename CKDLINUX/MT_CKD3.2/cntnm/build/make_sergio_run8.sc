## clean as much as possible
make -f makefile.common clean

## remove the .mod files in the */obj directories
/bin/rm -R *.obj/*

## remove the directories themselves
/bin/rmdir *.obj

## make!!
gmake -f make_cntnm_sergio_run8 linuxGNUdbl
