#!/bin/bash

## these come from running eg loop_filelist_gN.m
## -rw-rw-r-- 1 sergio pi_strow   10 Dec 13 17:56 g1_ir_list.txt
## -rw-rw-r-- 1 sergio pi_strow   60 Dec 13 17:56 gN_ir_list.txt
## -rw-r--r-- 1 sergio pi_strow 7340 Feb 13  2014 gN_ir_xseclist.txt

## N1 = number of lines in g1_ir_list.txt
## NN = number of lines in gN_ir_list.txt
## NX = number of lines in gN_ir_xseclist.txt

# run this with sbatch --array=1-N1   sergio_matlab_makegas.sbatch 1   WV
# run this with sbatch --array=1-NN   sergio_matlab_makegas.sbatch 10  GID 2-40
# run this with sbatch --array=1-NX   sergio_matlab_makegas.sbatch 51  GID 51-80

#  Name of the job:
#SBATCH --job-name=RUN8_rtp_profile

#  N specifies that 1 job step is to be allocated per instance of
#matlab
#SBATCH -N1

#  This specifies the number of cores per matlab session will be
#available for parallel jobs
#SBATCH --cpus-per-task 1

#  Specify the desired partition develop/batch/prod
##SBATCH --partition=batch
#SBATCH --partition=cpu2021

#### clust_runXtopts_savegasN_file : actually does the whole /home/sergio/SPECTRA/run8 shebang, slow and tedious
##SBATCH --qos=medium+
#SBATCH --time=11:59:00 
##SBATCH --time=15:59:00 
##SBATCH --mem-per-cpu=12000

## for most gases
##SBATCH --qos=short+
##SBATCH --time=0:59:00 

## for eg HN03 and SO2
#SBATCH --qos=medium+
#SBATCH --time=11:59:00 

#SBATCH --mem-per-cpu=12000

if [ $# -gt 0 ]; then
  echo "Your command line contains $# arguments"
elif [ $# -eq 0 ]; then
  echo "Your command line contains no arguments"
fi

####### calling srun twice for these paralel jobs might be asking for too much overhead
if [[ "$1" -eq "" ]]; then
  # this is default, make the spectral individually into gX.dat
  echo "cmd line arg = DNE, clust_runXtopts_savegas1_file"
  matlab -nodisplay -r "clust_runXtopts_savegas1_file; exit"
elif [[ "$1" -eq "1" ]]; then
  # this is default, make the spectral individually into gX.dat
  echo "cmd line arg = 1, clust_runXtopts_savegas1_file"
  matlab -nodisplay -r "clust_runXtopts_savegas1_file; exit"
elif [[ "$1" -eq "10" ]]; then
  # this is default, make the spectral individually into gX.dat
  echo "cmd line arg = 10, clust_runXtopts_savegasN_file"
  matlab -nodisplay -r "clust_runXtopts_savegasN_file; exit"
elif [[ "$1" -eq "51" ]]; then
  # this is default, put together the spectra into abs.dat
  echo "cmd line arg = 51, xscclust_runXtopts_mkgNvfiles"
  matlab -nodisplay -r "xscclust_runXtopts_mkgNvfiles; exit"  
else
  echo "error : use command line argument none,1  or 10 or 51"
fi
