#!/bin/sh
#
# build the local HITRAN datbase from HITRAN04 plus updates
# run in /asl/data/hitran

tr -d '\015\032' < HITRAN04/Lines/01_hit04.par   > h04.by.gas/g1.dat
tr -d '\015\032' < HITRAN04/Lines/02_hit04.par   > h04.by.gas/g2.dat
tr -d '\015\032' < HITRAN04/Lines/03_hit04.par   > h04.by.gas/g3.dat
tr -d '\015\032' < HITRAN04/Lines/04_hit04.par   > h04.by.gas/g4.dat
tr -d '\015\032' < HITRAN04/Lines/05_hit04.par   > h04.by.gas/g5.dat
tr -d '\015\032' < HITRAN04/Lines/06_hit04.par   > h04.by.gas/g6.dat
tr -d '\015\032' < HITRAN04/Lines/07_hit04.par   > h04.by.gas/g7.dat
tr -d '\015\032' < HITRAN04/Lines/08_hit04.par   > h04.by.gas/g8.dat
tr -d '\015\032' < HITRAN04/Lines/09_hit04.par   > h04.by.gas/g9.dat
tr -d '\015\032' < HITRAN04/Lines/10_hit04.par   > h04.by.gas/g10.dat
tr -d '\015\032' < HITRAN04/Lines/11_hit04.par   > h04.by.gas/g11.dat
tr -d '\015\032' < HITRAN04/Lines/12_hit04.par   > h04.by.gas/g12.dat
tr -d '\015\032' < HITRAN04/Lines/13_hit04.par   > h04.by.gas/g13.dat
tr -d '\015\032' < HITRAN04/Lines/14_hit04.par   > h04.by.gas/g14.dat
tr -d '\015\032' < HITRAN04/Lines/15_hit04.par   > h04.by.gas/g15.dat
tr -d '\015\032' < HITRAN04/Lines/16_hit04.par   > h04.by.gas/g16.dat
tr -d '\015\032' < HITRAN04/Lines/17_hit04.par   > h04.by.gas/g17.dat
tr -d '\015\032' < HITRAN04/Lines/18_hit04.par   > h04.by.gas/g18.dat
tr -d '\015\032' < HITRAN04/Lines/19_hit04.par   > h04.by.gas/g19.dat
tr -d '\015\032' < HITRAN04/Lines/20_hit04.par   > h04.by.gas/g20.dat
tr -d '\015\032' < HITRAN04/Lines/21_hit04.par   > h04.by.gas/g21.dat
tr -d '\015\032' < HITRAN04/Lines/22_hit04.par   > h04.by.gas/g22.dat
tr -d '\015\032' < HITRAN04/Lines/23_hit04.par   > h04.by.gas/g23.dat
tr -d '\015\032' < HITRAN04/Lines/24_hit04.par   > h04.by.gas/g24.dat
tr -d '\015\032' < HITRAN04/Lines/25_hit04.par   > h04.by.gas/g25.dat
tr -d '\015\032' < HITRAN04/Lines/26_hit04.par   > h04.by.gas/g26.dat
tr -d '\015\032' < HITRAN04/Lines/27_hit04.par   > h04.by.gas/g27.dat
tr -d '\015\032' < HITRAN04/Lines/28_hit04.par   > h04.by.gas/g28.dat
tr -d '\015\032' < HITRAN04/Lines/29_hit04.par   > h04.by.gas/g29.dat
tr -d '\015\032' < HITRAN04/Lines/30_hit04.par   > h04.by.gas/g30.dat
tr -d '\015\032' < HITRAN04/Lines/31_hit04.par   > h04.by.gas/g31.dat
tr -d '\015\032' < HITRAN04/Lines/32_hit04.par   > h04.by.gas/g32.dat
tr -d '\015\032' < HITRAN04/Lines/33_hit04.par   > h04.by.gas/g33.dat
tr -d '\015\032' < HITRAN04/Lines/34_hit04.par   > h04.by.gas/g34.dat
tr -d '\015\032' < HITRAN04/Lines/35_hit04.par   > h04.by.gas/g35.dat
tr -d '\015\032' < HITRAN04/Lines/36_hit04.par   > h04.by.gas/g36.dat
tr -d '\015\032' < HITRAN04/Lines/37_hit04.par   > h04.by.gas/g37.dat
tr -d '\015\032' < HITRAN04/Lines/38_hit04.par   > h04.by.gas/g38.dat
tr -d '\015\032' < HITRAN04/Lines/39_hit04.par   > h04.by.gas/g39.dat

tr -d '\015\032' < HITRAN04/Updates/01_hit06.par > h04.by.gas/g1.dat
tr -d '\015\032' < HITRAN04/Updates/04_hit06.par > h04.by.gas/g4.dat
tr -d '\015\032' < HITRAN04/Updates/12_hit06.par > h04.by.gas/g12.dat
tr -d '\015\032' < HITRAN04/Updates/19_hit06.par > h04.by.gas/g19.dat
tr -d '\015\032' < HITRAN04/Updates/27_hit05.par > h04.by.gas/g27.dat
tr -d '\015\032' < HITRAN04/Updates/28_hit06.par > h04.by.gas/g28.dat
tr -d '\015\032' < HITRAN04/Updates/36_hit06.par > h04.by.gas/g36.dat

