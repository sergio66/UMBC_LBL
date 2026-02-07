#!/bin/bash

STR="compiling TIPS_2024_allisotopes.f"
echo $STR
# f77 -o TIPS_2024_allisotopes.x -N109 -W -O3 TIPS_2024_allisotopes.f or
# /asl/opt/absoft/absoft10.0/bin/af77 -o TIPS_2024_allisotopes.x -N109 -W -O3 TIPS_2024_allisotopes.f
# # -N109 fold all names to upper case
# # -W    wide source file
# # -O    some optimizations
ifort -o TIPS_2024_allisotopes.x -names uppercase -extend-source 132 -O3 TIPS_2024_allisotopes.f 

########################################################################
