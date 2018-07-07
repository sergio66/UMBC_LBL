c  ****************************************
      Block Data ISO_2002
c  ****************************************
      implicit DOUBLE PRECISION (a-h,o-z)
c
      INCLUDE 'Species_2002.cmn'
      INCLUDE 'Isotops_2002.cmn'
c
c    The number of isotopes for a particular molecule:
      DATA (ISONM(I),I=1,NMOL)/
c     H2O, CO2, O3, N2O, CO, CH4, O2,
     +  6,   9,  18,   5,  6,   3,  3,
c      NO, SO2, NO2, NH3, HNO3, OH, HF, HCl, HBr, HI,
     +  3,   2,   1,   2,    1,  3,  1,   2,   2,  1,
c     ClO, OCS, H2CO, HOCl, N2, HCN, CH3Cl, H2O2, C2H2, C2H6, PH3
     +  2,   5,    3,    2,  1,   3,     2,    1,    2,    1,   1,       4/24/97
c     COF2, SF6, H2S, HCOOH, HO2, O, ClONO2,  NO+, HOBr, C2H4
     +   1,   1,   3,     1,   1, 1,      2,    1,    2,    2/
c

c       H2O
      DATA (ISO82(1,J),J=1,6)/
     +  161,181,171,162,182,172/

      DATA (ISO82(2,J),J=1,9)/
c       CO2
     +  626,636,628,627,638,637,828,728,727/
c       O3
      DATA (ISO82(3,J),J=1,18)/
     +  666,668,686,667,676,886,868,678,768,
     +  786,776,767,888,887,878,778,787,777/

      DATA (ISO82(4,J),J=1,5)/
c       N2O
     +  446,456,546,448,447/

      DATA (ISO82(5,J),J=1,6)/
c       CO,                 
     +  26,36,28,27,38,37/  

      DATA (ISO82(6,J),J=1,3)/
c       CH4
     +  211,311,212/

      DATA (ISO82(7,J),J=1,3)/
c       O2        
     +  66,68,67/  

      DATA (ISO82(8,J),J=1,3)/
c       NO
     +  46,56,48/

      DATA (ISO82(9,J),J=1,2)/
c       SO2
     +  626,646/

      DATA (ISO82(10,J),J=1,1)/
c      NO2 
     + 646/  

      DATA (ISO82(11,J),J=1,2)/
c       NH3
     +  4111,5111/

      DATA (ISO82(12,J),J=1,1)/
c       HNO3
     +  146/

      DATA (ISO82(13,J),J=1,3)/
c       OH
     +  61,81,62/

      DATA (ISO82(14,J),J=1,1)/
c       HF
     +  19/

      DATA (ISO82(15,J),J=1,2)/
c       HCl
     +  15,17/

      DATA (ISO82(16,J),J=1,2)/
c       HBr
     +  19,11/

      DATA (ISO82(17,J),J=1,1)/
c       HI
     +  17/

      DATA (ISO82(18,J),J=1,2)/
c       ClO,  
     +  56,76/

      DATA (ISO82(19,J),J=1,5)/
c       OCS              
     +  622,624,632,623,822/

      DATA (ISO82(20,J),J=1,3)/
c       H2CO
     +  126,136,128/

      DATA (ISO82(21,J),J=1,2)/
c       HOCl,    
     +  165,167/

      DATA (ISO82(22,J),J=1,1)/
c       N2
     +  44/

      DATA (ISO82(23,J),J=1,3)/
c       HCN
     +  124,134,125/

      DATA (ISO82(24,J),J=1,2)/
c       CH3Cl
     +  215,217/

      DATA (ISO82(25,J),J=1,1)/
c       H2O2
     +  1661/

      DATA (ISO82(26,J),J=1,2)/
c       C2H2       
     +  1221,1231/

      DATA (ISO82(27,J),J=1,1)/
c       C2H6
     +  1221/ 

      DATA (ISO82(28,J),J=1,1)/
c       PH3
     =  1111/

      DATA (ISO82(29,J),J=1,1)/
c     COF2
     + 269/ 

      DATA (ISO82(30,J),J=1,1)/
c       SF6
     +  29/

      DATA (ISO82(31,J),J=1,3)/
c       H2S         
     +  121,141,131/

      DATA (ISO82(32,J),J=1,1)/
c       HCOOH 
     +  126/

      DATA (ISO82(33,J),J=1,1)/
c       HO2
     +  166/

      DATA (ISO82(34,J),J=1,1)/
c       O  
     +  6/

      DATA (ISO82(35,J),J=1,2)/
c       ClONO2     
     +  5646,7646/

      DATA (ISO82(36,J),J=1,1)/
c       NO+
     +  46/

      DATA (ISO82(37,J),J=1,2)/
c      HOBr
     + 169,161/

      DATA (ISO82(38,J),J=1,2)/
c       C2H4
     +  221,231/

      end
