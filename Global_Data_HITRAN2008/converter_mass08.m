%% this reads in molparam1.txt, which is a "simpler" version of molparam.txt
%% it then outputs a text file, in the same fashion as massXY.dat

clear all

dd = load('molparam1.txt');

%% there have to be 42 molecules, see molparam.txt

ix5 = dd(:,5); boo = find(ix5 == -1); 

if length(boo) ~= 42
  error('oops : did not find 42 gases!!!')
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the number of isotopes per molecule
for ii = 1 : length(boo)-1
  gasid(ii)       = ii;
  numisotopes(ii) = boo(ii+1)-boo(ii)-1;
  neg1(ii)        = -1';
  zer(ii)         = 0;
  end
ii = 42;
gasid(ii)       = ii;
numisotopes(ii) = max(size(dd)) - boo(ii);
neg1(ii)        = -1;
zer(ii)         = 0;

matrixA = [gasid; numisotopes; neg1; zer; zer];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the [masses abundances gj Q296] of the isotopes
ij = 0;
for ii = 1 : length(boo)-1
  woof = boo(ii)+1 : boo(ii+1)-1;
  woof = dd(woof,[5 2 4 3]);
  for jj = 1 : matrixA(2,ii)
    ij = ij + 1;
    gasid(ij)  = ii;
    molarmass(ij) = woof(jj,1);
    abundance(ij) = woof(jj,2);
    gj(ij)        = woof(jj,3);
    q296(ij)      = woof(jj,4);
    end
  end
ii = 42;
  woof = boo(ii)+1 : max(size(dd));
  woof = dd(woof,[5 2 4 3]);
  for jj = 1 : matrixA(2,ii)
    ij = ij + 1;
    gasid(ij)  = ii;
    molarmass(ij) = woof(jj,1);
    abundance(ij) = woof(jj,2);
    gj(ij)        = woof(jj,3);
    q296(ij)      = woof(jj,4);
    end

matrixB = [gasid; molarmass; abundance; gj; q296];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('mass08.dat','w');
%fprintf(fid,' \% see /home/sergio/SPECTRA/Global_Data_HITRAN2008/converter_mass08.m \n');
fprintf(fid,'     %4i    %4i   %10.6e %10.6e %10.6e\n',matrixA);
fprintf(fid,'  %4i    %10.6f   %10.6e  %4i   %10.6e   \n',matrixB);
fclose(fid);

