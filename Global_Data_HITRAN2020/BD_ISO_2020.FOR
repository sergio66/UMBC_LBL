c  ****************************************
      Block Data ISO_2020
c  ****************************************
      implicit DOUBLE PRECISION (a-h,o-z)
c
c...date last changed 1 December, 2020
      INCLUDE 'SPECIES_2020.CMN'
      INCLUDE 'ISOTOPS.CMN'
c
c    The number of isotopes for a particular molecule:
      DATA (ISONM(I),I=1,NMOL)/
c  #    1    2   3    4   5    6   7 
c     H2O, CO2, O3, N2O, CO, CH4, O2,
     +  9,  13, 18,   5,  9,  4,  6,
     
c  #    8    9   10   11    12  13  14   15   16  17
c      NO, SO2, NO2, NH3, HNO3, OH, HF, HCl, HBr, HI,
     +  3,   4,   2,   2,    2,  3,  2,   4,   4,  2,
     
c  #   18   19    20    21  22   23     24    25    26    27   28
c     ClO, OCS, H2CO, HOCl, N2, HCN, CH3Cl, H2O2, C2H2, C2H6, PH3
     +  2,   6,    3,    2,  3,   3,     2,    1,    3,    3,   1,  
     
c  #    29,  30   31     32   33 34      35    36    37    38    
c     COF2, SF6, H2S, HCOOH, HO2, O, ClONO2,  NO+, HOBr, C2H4
     +   2,   1,   3,     1,   1, 1,      2,    1,    2,    3,
     
c  #     39     40     41   42    43    44    45    46   47   
c     CH3OH, CH3Br, CH3CN, CF4, C4H2, HC3N,   H2,   CS,  SO3,
     +   1,      2,     4,   1,    1,    6,    2,    4,   1, 

c  #     48     49     50     51    52     53    54     55   56    57
c      C2N2,  COCL2,   SO,  CH3F   GeH4,  CS2,  CH3I,  NF3, C3H4, CH3
     +   2,      2,     3,     1,    5,     4,     2,    1,    1,   1/
c
c

c 1     H2O
      DATA (ISO82(1,J),J=1,9)/
     +  161,181,171,162,182,172,262,282,272/

      DATA (ISO82(2,J),J=1,13)/
c 2     CO2
     +  626,636,628,627,638,637,828,827,727,838,837,737,646/
c 3     O3
      DATA (ISO82(3,J),J=1,18)/
     +  666,668,686,667,676,688,868,678,768,
     +  687,677,767,888,788,878,778,787,777/

      DATA (ISO82(4,J),J=1,5)/
c 4     N2O
     +  446,456,546,448,447/

      DATA (ISO82(5,J),J=1,9)/
c 5     CO,                 
     +  26,36,28,27,38,37,46,48,47/  

      DATA (ISO82(6,J),J=1,4)/
c 6    CH4
     +  211,311,212,312/

      DATA (ISO82(7,J),J=1,6)/
c 7     O2        
     +  66,68,67,88,87,77/  

      DATA (ISO82(8,J),J=1,3)/
c 8     NO
     +  46,56,48/

      DATA (ISO82(9,J),J=1,4)/
c 9     SO2
     +  626,646,636,628/

      DATA (ISO82(10,J),J=1,2)/
c10    NO2 
     + 646,656/  

      DATA (ISO82(11,J),J=1,2)/
c11     NH3
     +  4111,5111/

      DATA (ISO82(12,J),J=1,2)/
c12     HNO3
     +  146,156/

      DATA (ISO82(13,J),J=1,3)/
c13     OH
     +  61,81,62/

      DATA (ISO82(14,J),J=1,2)/
c14     HF
     +  19,29/

      DATA (ISO82(15,J),J=1,4)/
c15     HCl
     +  15,17,25,27/

      DATA (ISO82(16,J),J=1,4)/
c16     HBr
     +  19,11,21,29/

      DATA (ISO82(17,J),J=1,2)/
c17     HI
     +  17,27/

      DATA (ISO82(18,J),J=1,2)/
c18     ClO,  
     +  56,76/

      DATA (ISO82(19,J),J=1,6)/
c19     OCS              
     +  622,624,632,623,822,634/

      DATA (ISO82(20,J),J=1,3)/
c20     H2CO
     +  126,136,128/

      DATA (ISO82(21,J),J=1,2)/
c21     HOCl,    
     +  165,167/

      DATA (ISO82(22,J),J=1,3)/
c22     N2
     +  44,45,55/

      DATA (ISO82(23,J),J=1,3)/
c23     HCN
     +  124,134,125/

      DATA (ISO82(24,J),J=1,2)/
c24     CH3Cl
     +  215,217/

      DATA (ISO82(25,J),J=1,1)/
c25     H2O2
     +  1661/

      DATA (ISO82(26,J),J=1,3)/
c26     C2H2       
     +  1221,1231,1222/

      DATA (ISO82(27,J),J=1,3)/
c27    C2H6
     +  1221,1231,1222/ 

      DATA (ISO82(28,J),J=1,1)/
c28     PH3
     =  1111/

      DATA (ISO82(29,J),J=1,2)/
c29   COF2
     + 269,369/ 

      DATA (ISO82(30,J),J=1,1)/
c30     SF6
     +  29/

      DATA (ISO82(31,J),J=1,3)/
c31     H2S         
     +  121,141,131/

      DATA (ISO82(32,J),J=1,1)/
c32     HCOOH 
     +  126/

      DATA (ISO82(33,J),J=1,1)/
c33     HO2
     +  166/

      DATA (ISO82(34,J),J=1,1)/
c34     O  
     +  6/

      DATA (ISO82(35,J),J=1,2)/
c35     ClONO2     
     +  5646,7646/

      DATA (ISO82(36,J),J=1,1)/
c36     NO+
     +  46/

      DATA (ISO82(37,J),J=1,2)/
c37    HOBr
     + 169,161/

      DATA (ISO82(38,J),J=1,3)/
c38     C2H4
     +  112211,112311,112212/

      DATA (ISO82(39,J),J=1,1)/
c39     CH3OH
     +  2161/

      DATA (ISO82(40,J),J=1,2)/
c40     CH3Br
     +  219,211/

      DATA (ISO82(41,J),J=1,4)/
c41     CH3CN
     +  2124,2134,3124,3134/

      DATA (ISO82(42,J),J=1,1)/
c42     CF4
     +  29/

      DATA (ISO82(43,J),J=1,1)/
c43     C4H2
     +  2211/

      DATA (ISO82(44,J),J=1,6)/
c44     HC3N
     +  12224,12225,12234,12324,13224,22224/

      DATA (ISO82(45,J),J=1,2)/
c45     H2
     +  11,12/

      DATA (ISO82(46,J),J=1,4)/
c46     CS
     +  22,24,32,23/

      DATA (ISO82(47,J),J=1,1)/
c47       SO3
     +  26/

      DATA (ISO82(48,J),J=1,2)/
c48     C2N2
     +  4224,5225/

c...new 2020
      DATA (ISO82(49,J),J=1,2)/
c49       COCl2
     +  2655,2657/

      DATA (ISO82(50,J),J=1,3)/
c50     SO
     +  26,46,28/

      DATA (ISO82(51,J),J=1,1)/
c51     CH3F
     +  219/

      DATA (ISO82(52,J),J=1,5)/
c52     GeH4
     +  7411,7211,7011,7311,7611/

      DATA (ISO82(53,J),J=1,4)/
c53       CS2
     +  222,224,223,232/

      DATA (ISO82(54,J),J=1,2)/
c54       CH3I
     +  217,317/

      DATA (ISO82(55,J),J=1,1)/
c55     NF3
     +  49/

       DATA (ISO82(56,J),J=1,1)/
c56     C3H4
     +  1221/

      DATA (ISO82(57,J),J=1,1)/
c57     CH3
     +  2111/




      end
