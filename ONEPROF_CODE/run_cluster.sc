## this is set up for glist = 1 2 3 4 5 6 9 12] 
## and 1205-1405 cm-1 == 8 chunks   .. so runwater 1-8, other 7 gases * 8 chunks

/bin/rm slurm*.out; /bin/rm slurm*.err
sbatch --array=1-8  sergio_matlab_makegas.sbatch 1
sbatch --array=1-56 sergio_matlab_makegas.sbatch 10
