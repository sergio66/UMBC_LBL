c WARNING : IF YOU EDIT THIS, THEN EDIT IncludeMaxer.m as well!!!!!

c this is max length of arrays that can be used in the Mex Files  
c this number came out of 
c   200000 = max number of elements in mesh 
c              eg (755-655)/0.0005 = 200000
c        4 = number tacked on to arrays so boxint(y,5) can be done  
      integer MaxLen
c      parameter(MaxLen=200010)
      parameter(MaxLen=2500010)

c assume max number of any of P,Q,R lines = 300
      integer MaxPQR
      parameter(MaxPQR=300)

c assume max number of any of layers = 100
      integer kMaxLayer
      parameter(kMaxLayer=200)

c assume max number of lines in a band = 10000
      integer kMaxBandLines
      parameter(kMaxBandLines=2500000)
c      parameter(kMaxBandLines=50000)

c assume max number of isotopes per molcule = 20
      integer kMaxIsotopes
      parameter(kMaxIsotopes=20)

c this is for MTCKD 2.5 : 20000 pts can be output at 10 cm-1 res
      integer kCKDcoarse
      parameter (kCKDcoarse = 25000)
