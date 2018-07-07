c WARNING : IF YOU EDIT THIS, THEN EDIT IncludeMaxer.m as well!!!!!

c this is max length of arrays that can be used in the Mex Files  
c this number came out of 
c   200000 = max number of elements in mesh 
c              eg (755-655)/0.0005 = 200000
c              eg (755-705)/0.0005 = 100000
c              eg (755-730)/0.0005 = 50000
c        4 = number tacked on to arrays so boxint(y,5) can be done  
      integer MaxLen,MaxLenJMH
      parameter(MaxLen=200010)
      parameter(MaxLenJMH=4450010)

c assume max number of any of P,Q,R lines = 300
      integer MaxPQR
      parameter(MaxPQR=300)

c assume max number of any of layers = 100
      integer kMaxLayer
      parameter(kMaxLayer=200)

c assume max number of lines in a band = 10000
      integer kMaxBandLines
c      parameter(kMaxBandLines=50000)   !! before Nov 2007
      parameter(kMaxBandLines=200000)

c assume max number of isotopes per molcule = 20
      integer kMaxIsotopes
      parameter(kMaxIsotopes=20)
