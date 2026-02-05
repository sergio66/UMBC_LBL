#!/bin/sh
#
# build the local HITRAN datbase from HITRAN08 plus updates
# run in /asl/data/hitran

tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/01_hit08.par     > g1.dat
tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/02_hit08_f53.par > g2_origH08.dat
tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/03_hit08_f53.par > g3.dat
tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/04_hit08.par     > g4.dat
tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/05_hit08.par     > g5.dat
tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/06_hit08_f53.par > g6.dat
tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/07_hit08.par     > g7.dat
tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/08_hit08.par     > g8.dat
tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/09_hit08.par     > g9.dat
tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/10_hit08.par     > g10.dat
tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/11_hit08.par     > g11.dat
tr -d '\015\032' < /asl/data/hitran/H2008/HITRAN2008/By-Molecule/Uncompressed-files/12_hit08.par     > g12.dat
