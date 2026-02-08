The file TIPS_2024_1p2.pdf contains an overview of the TIPS for each molecule in 
HITRAN2024. The data contained is based on version 1.2 of the compilation and also
available from Zenodo (https://zenodo.org/records/17191976).

The file TIPS_2024_FORTRAN_v1p2.zip is a compressed file of the FORTRAN codes needed to 
compile TIPS_2024_v1p2.FOR.  TIPS_2024_v1p2.FOR is a stand-alone program to compute the 
Total Internal Partition Sums (TIPS) for all the molecules present on HITRAN2024, plus 
TIPS for C3H8.  TIPS_2024_v1p2.FOR contains the subroutine BD_TIPS_2024 that users can 
insert into their codes to obtain Q(T) for their applications. Using Intel FORTRAN the 
compilation is ifort -w -132 TIPS_2024_v1p2.FOR -o TIPS_2024_v1p2.EXE

The file TIPS_2024_v1p2_py.zip contains the file TIPS_2024_v1p2.py, which is the 
python language version of the TIPS_2024_v1p2.FOR code. and a folder "data", which contains 
the python dictionaries needed to run the code.  Unzip the file to create a folder 
TIPS_2024_v1p2_python.  In a command window go to that folder and type 
python3 TIPS_2024_v1p2.py to execute it.

The reference for the codes is Partition sums for molecules and their isotopologues for 
HITRAN2024:
-- Robert R. Gamache, Bastien Vispoel, Jonathan Tennyson, Sergei N. Yurchenko, 
   Oleg L. Polyansky, Iouli E. Gordon, Robert J. Hargreaves, Xinchuan Huang, Journal of 
   Quantitative Spectroscopy & Radiative Transfer 345 (2025) 109568.  
   https://doi.org/10.1016/j.jqsrt.2025.109568
-- The Corrigendum is provided at Journal of Quantitative Spectroscopy & Radiative 
   Transfer 346, (2025) 109618.
   https://doi.org/10.1016/j.jqsrt.2025.109618
 