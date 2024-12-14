#!/bin/bash

# run this with sbatch --array=1-10000 --output='testslurm' sergio_matlab_.sbatch 1   WV
# run this with sbatch --array=1-10000                      sergio_matlab_.sbatch 10  GID 2-40
# run this with sbatch --array=1-10000                      sergio_matlab_.sbatch 50  GID 51-80
#

#  Name of the job:
#SBATCH --job-name=RUN8_rtp_profile

#  N specifies that 1 job step is to be allocated per instance of
#matlab
#SBATCH -N1

#  This specifies the number of cores per matlab session will be
#available for parallel jobs
#SBATCH --cpus-per-task 1

#  Specify the desired partition develop/batch/prod
#SBATCH --partition=batch

#### clust_runXtopts_savegasN_file : actually does the whole /home/sergio/SPECTRA/run8 shebang, slow and tedious
#SBATCH --qos=medium_prod+
#SBATCH --time=11:59:00 
##SBATCH --time=15:59:00 
#SBATCH --mem-per-cpu=12000

#### clust_runXtopts_mkgNvfiles : puts together the Toffset files so you can start compressing, pretty fast
##SBATCH --qos=short+
##SBATCH --time=0:59:00 
##SBATCH --mem-per-cpu=12000

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
fi
