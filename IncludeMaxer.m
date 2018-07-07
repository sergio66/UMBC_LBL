% copied from max.inc
% this is max length of arrays that can be used in the Mex Files  
% this number came out of 
%   200000 = max number of elements in mesh 
%              eg (755-655)/0.0005 = 200000
%        4 = number tacked on to arrays so boxint(y,5) can be done  
%      integer MaxLen
%      parameter(MaxLen=200010)
MaxLen=200010;

% assume max number of any of P,Q,R lines = 300
%      integer MaxPQR
%      parameter(MaxPQR=300)
MaxPQR=300;

% assume max number of any of layers = 100
%      integer kMaxLayer
%      parameter(kMaxLayer=100)
kMaxLayer=200;

% assume max number of lines in a band = 10000
%      integer kMaxBandLines
%      parameter(kMaxBandLines=50000)
kMaxBandLines=50000;

% assume max number of isotopes per molcule = 20
%      integer kMaxIsotopes
%      parameter(kMaxIsotopes=20)
kMaxIsotopes=20;
