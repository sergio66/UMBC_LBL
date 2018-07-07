function [line,num_band,PQR,bandtype,lengthless]=...
     removeCO2linesNEW(band,lineORIG,low,high,num_band,PQR,bandtype,...
                    vers,strengthM,homepath,hitlin_fname,xfar,...
                    PQRallowed,exchangelinecenters); 

% this function removes the relevant band lines if the band lies between 
% wavenumbers low,high
% NOTE THAT low  = low_user  - 25
%           high = high_user + 25
%num_bad tells us how many bands to do
%PQR is a code that tells us which bands to add in  01 = Q delt pi
%                                                   02 = Q sig  pi
%                                                  -14 = P sig  pi
%                                                  -15 = P delt pi
%                                                  +14 = R sig  pi
%                                                  +15 = R delt pi

%                                                  -11 = P sig  sig
%                                                  -12 = P delt delt
%                                                  -13 = P pi   pi

%                                                  +11 = R sig  sig
%                                                  +12 = R delt delt
%                                                  +13 = R pi   pi

% so PQRallowed tells us which ones of the above to allow

%bandtype tells you which band to do eg 740 etc

%  Q_deltpi        =  [668 740 2093];
%  Q_deltpi_lower  =  [2   4   2   ];
%  Q_deltpi_upper  =  [4   8   14  ];
%  Q_isotope       =  [1   1   1   ];

%  Q_sigpi        =  [618 648 662 667 720 791 1932 2080 2129];
%  Q_sigpi_lower  =  [2   1   1   1   2   3   1    1    2 ];
%  Q_sigpi_upper  =  [3   2   2   2   5   8   6    8    15];
%  Q_isotope      =  [1   2   3   1   1   1   1    1    1 ];

%  P_sigsig        =  [2350 2351 2352 2353 2354];
%  P_sigsig_lower  =  [1    1    1     3    5];
%  P_sigsig_upper  =  [9    9    9     23   25];
%  P_isotope      =   [1    2    3     1    1];

%  P_deltdelt        =  [2310 2311];
%  P_deltdelt_lower  =  [4    4];
%  P_deltdelt_upper  =  [24   24];
%  P_isotope         =  [1    2];

%  P_pipi        =  [2320 2321 2322];
%  P_pipi_lower  =  [2    2    2];
%  P_pipi_upper  =  [16   16   16];
%  P_isotope      = [1    2    3];

lengthless = 0;

global quiet p2311_21_jmax p2350_jmax r2350_jmax pr2351_jmax

iii    = 0;
iFound = -1;          %assume nothing found

ugh = -1; %do not do the PR 15 um band line mixing
ugh =  0; %only do Q720,741,791 bands
ugh = +1; %do everything 

ugh = +1;      %%%%%%% <<<<<<<< ---------------------

if (ugh > 0)
  %ZZ
  %this would only do the R sig sig branch!!!!!!
  %IndicesToSearch=[01 02     11 -12 12 -13 13 -14 14 -15 15];
  %this would only do the P pi pi branch!!!!!!
  % IndicesToSearch=[01 02 -11 11 -12 12 -13    -14 14 -15 15];
  %do not do PR sigsig, R pipi, R deltdelt
  %IndicesToSearch=[01 02        -12    -13    -14 14 -15 15];
  %this is *****all****
  IndicesToSearch = [01 02 -11 11 -12 12 -13 13 -14 14 -15 15];
elseif (ugh < 0)
  IndicesToSearch = [01 02 -11 11 -12 12 -13 13];
elseif(ugh ==  0)
  IndicesToSearch = [01 02];
end

IndicesToSearch = PQRallowed;

ident = [01];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  Q_deltpi        =  [668 740 2093];
  Q_deltpi_lower  =  [2   4   2   ];
  Q_deltpi_upper  =  [4   8   14  ];
  Q_isotope       =  [1   1   1   ];
  str.name(iii,1:10) = 'Q_deltpi  '; 
  str.ID(iii,1)      = 01;   
  str.PQRtype(iii,1) = 'Q';
  str.number(iii,1)  = 3;
  str.iso(iii,1:str.number(iii,1))      = Q_isotope;
  str.bbb(iii,1:str.number(iii,1))      = Q_deltpi;       
  str.thelower(iii,1:str.number(iii,1)) = Q_deltpi_lower; 
  str.theupper(iii,1:str.number(iii,1)) = Q_deltpi_upper;
end

ident = [02];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  Q_sigpi        =  [618 648 662 667 720 791 1932 2080 2129];
  Q_sigpi_lower  =  [2   1   1   1   2   3   1    1    2 ];
  Q_sigpi_upper  =  [3   2   2   2   5   8   6    8    15];
  Q_isotope      =  [1   2   3   1   1   1   1    1    1 ];
  str.name(iii,1:10) = 'Q_sigpi   '; 
  str.ID(iii,1)      = 02;   
  str.PQRtype(iii,1) = 'Q';
  str.number(iii,1)  = 9;
  str.iso(iii,1:str.number(iii,1))      = Q_isotope;
  str.bbb(iii,1:str.number(iii,1))      = Q_sigpi;       
  str.thelower(iii,1:str.number(iii,1)) = Q_sigpi_lower; 
  str.theupper(iii,1:str.number(iii,1)) = Q_sigpi_upper;
end

ident = [-14];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  P_sigpi        =  [618 648 662 667 720 791 ];
  P_sigpi_lower  =  [2   1   1   1   2   3   ];
  P_sigpi_upper  =  [3   2   2   2   5   8   ];
  P_isotope      =  [1   2   3   1   1   1   ];
  str.name(iii,1:10) = 'P_sigpi   '; 
  str.ID(iii,1)      = -14;   
  str.PQRtype(iii,1) = 'P';
  str.number(iii,1)  = 6;
  str.iso(iii,1:str.number(iii,1))      = P_isotope;
  str.bbb(iii,1:str.number(iii,1))      = P_sigpi;
  str.thelower(iii,1:str.number(iii,1)) = P_sigpi_lower; 
  str.theupper(iii,1:str.number(iii,1)) = P_sigpi_upper;
end

ident = [+14];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  R_sigpi        =  [618 648 662 667 720 791 ];
  R_sigpi_lower  =  [2   1   1   1   2   3   ];
  R_sigpi_upper  =  [3   2   2   2   5   8   ];
  R_isotope      =  [1   2   3   1   1   1   ];
  str.name(iii,1:10) = 'R_sigpi   '; 
  str.ID(iii,1)      = +14;   
  str.PQRtype(iii,1) = 'R';
  str.number(iii,1)  = 6;
  str.iso(iii,1:str.number(iii,1))      = R_isotope;
  str.bbb(iii,1:str.number(iii,1))      = R_sigpi;
  str.thelower(iii,1:str.number(iii,1)) = R_sigpi_lower; 
  str.theupper(iii,1:str.number(iii,1)) = R_sigpi_upper;
end

ident = [-15];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  P_deltpi        =  [668 740];
  P_deltpi_lower  =  [2   4  ];
  P_deltpi_upper  =  [4   8  ];
  P_isotope      =   [1   1  ];
  str.name(iii,1:10) = 'P_deltpi  '; 
  str.ID(iii,1)      = -15;   
  str.PQRtype(iii,1) = 'P';
  str.number(iii,1)  = 2;
  str.iso(iii,1:str.number(iii,1))      = P_isotope;
  str.bbb(iii,1:str.number(iii,1))      = P_deltpi;       
  str.thelower(iii,1:str.number(iii,1)) = P_deltpi_lower; 
  str.theupper(iii,1:str.number(iii,1)) = P_deltpi_upper;
end

ident = [+15];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  R_deltpi        =  [668 740];
  R_deltpi_lower  =  [2   4  ];
  R_deltpi_upper  =  [4   8  ];
  R_isotope      =   [1   1  ];
  str.name(iii,1:10) = 'R_deltpi  '; 
  str.ID(iii,1)      = +15;   
  str.PQRtype(iii,1) = 'R';
  str.number(iii,1)  = 2;
  str.iso(iii,1:str.number(iii,1))      = R_isotope;
  str.bbb(iii,1:str.number(iii,1))      = R_deltpi;       
  str.thelower(iii,1:str.number(iii,1)) = R_deltpi_lower; 
  str.theupper(iii,1:str.number(iii,1)) = R_deltpi_upper;
end

ident = [-11];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  P_sigsig        =  [2350 2351 2352 2353 2354];
  P_sigsig_lower  =  [1    1    1     3    5  ];
  P_sigsig_upper  =  [9    9    9     23   25 ];
  P_isotope      =   [1    2    3     1    1  ];
  str.name(iii,1:10) = 'P_sigsig  '; 
  str.ID(iii,1)      = -11;   
  str.PQRtype(iii,1) = 'P';
  str.number(iii,1)  = 5;
  str.iso(iii,1:str.number(iii,1))      = P_isotope;
  str.bbb(iii,1:str.number(iii,1))      = P_sigsig;       
  str.thelower(iii,1:str.number(iii,1)) = P_sigsig_lower; 
  str.theupper(iii,1:str.number(iii,1)) = P_sigsig_upper;
end

ident = [+11];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  R_sigsig        =  [2350 2351 2352 2353 2354];
  R_sigsig_lower  =  [1    1    1    3     5  ];
  R_sigsig_upper  =  [9    9    9    23    25 ];
  R_isotope      =   [1    2    3    1     1  ];
  str.name(iii,1:10) = 'R_sigsig  '; 
  str.ID(iii,1)      = +11;   
  str.PQRtype(iii,1) = 'R';
  str.number(iii,1)  = 5;
  str.iso(iii,1:str.number(iii,1))      = R_isotope;
  str.bbb(iii,1:str.number(iii,1))      = R_sigsig;       
  str.thelower(iii,1:str.number(iii,1)) = R_sigsig_lower; 
  str.theupper(iii,1:str.number(iii,1)) = R_sigsig_upper;
end

%ident = [-11];
%if (sum(ismember(ident,IndicesToSearch)) == 1)
%  iii = iii+1;
%  P_sigsig        =  [2350 2351 2352 2353 2354 2355];
%  P_sigsig_lower  =  [1    1    1     3    5    1];
%  P_sigsig_upper  =  [9    9    9     23   25   9];
%  P_isotope      =   [1    2    3     1    1    4];
%  str.name(iii,1:10) = 'P_sigsig  '; 
%  str.ID(iii,1)      = -11;   
%  str.PQRtype(iii,1) = 'P';
%  str.number(iii,1)  = 6;
%  str.iso(iii,1:str.number(iii,1))      = P_isotope;
%  str.bbb(iii,1:str.number(iii,1))      = P_sigsig;       
%  str.thelower(iii,1:str.number(iii,1)) = P_sigsig_lower; 
%  str.theupper(iii,1:str.number(iii,1)) = P_sigsig_upper;
%end

%ident = [+11];
%if (sum(ismember(ident,IndicesToSearch)) == 1)
%  iii = iii+1;
%  R_sigsig        =  [2350 2351 2352 2353 2354 2355];
%  R_sigsig_lower  =  [1    1    1    3     5    1];
%  R_sigsig_upper  =  [9    9    9    23    25   9];
%  R_isotope      =   [1    2    3    1     1    4];
%  str.name(iii,1:10) = 'R_sigsig  '; 
%  str.ID(iii,1)      = +11;   
%  str.PQRtype(iii,1) = 'R';
%  str.number(iii,1)  = 6;
%  str.iso(iii,1:str.number(iii,1))      = R_isotope;
%  str.bbb(iii,1:str.number(iii,1))      = R_sigsig;       
%  str.thelower(iii,1:str.number(iii,1)) = R_sigsig_lower; 
%  str.theupper(iii,1:str.number(iii,1)) = R_sigsig_upper;
%end

ident = [-12];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  P_deltdelt        =  [2310 2311];
  P_deltdelt_lower  =  [4    4];
  P_deltdelt_upper  =  [24   24];
  P_isotope         =  [1    2];
  str.name(iii,1:10) = 'P_deltdelt'; 
  str.ID(iii,1)      = -12;   
  str.PQRtype(iii,1) = 'P';
  str.number(iii,1)  = 2;
  str.iso(iii,1:str.number(iii,1))      = P_isotope;
  str.bbb(iii,1:str.number(iii,1))      = P_deltdelt;       
  str.thelower(iii,1:str.number(iii,1)) = P_deltdelt_lower; 
  str.theupper(iii,1:str.number(iii,1)) = P_deltdelt_upper;
end

ident = [+12];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  R_deltdelt        =  [2310 2311];
  R_deltdelt_lower  =  [4    4];
  R_deltdelt_upper  =  [24   24];
  R_isotope        =   [1    2];
  str.name(iii,1:10) = 'R_deltdelt'; 
  str.ID(iii,1)      = +12;   
  str.PQRtype(iii,1) = 'R';
  str.number(iii,1)  = 2;
  str.iso(iii,1:str.number(iii,1))      = R_isotope;
  str.bbb(iii,1:str.number(iii,1))      = R_deltdelt;       
  str.thelower(iii,1:str.number(iii,1)) = R_deltdelt_lower; 
  str.theupper(iii,1:str.number(iii,1)) = R_deltdelt_upper;
end

ident = [-13];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  P_pipi        =  [2320 2321 2322];
  P_pipi_lower  =  [2    2    2];
  P_pipi_upper  =  [16   16   16];
  P_isotope      = [1    2    3];
  str.name(iii,1:10) = 'P_pipi    '; 
  str.ID(iii,1)      = -13;   
  str.PQRtype(iii,1) = 'P';
  str.number(iii,1)  = 3;
  str.iso(iii,1:str.number(iii,1))      = P_isotope;
  str.bbb(iii,1:str.number(iii,1))      = P_pipi;       
  str.thelower(iii,1:str.number(iii,1)) = P_pipi_lower; 
  str.theupper(iii,1:str.number(iii,1)) = P_pipi_upper;
end

ident = [+13];
if (sum(ismember(ident,IndicesToSearch)) == 1)
  iii = iii+1;
  R_pipi        =  [2320 2321 2322];
  R_pipi_lower  =  [2    2    2];
  R_pipi_upper  =  [16   16   16];
  R_isotope     =  [1    2    3];
  str.name(iii,1:10) = 'R_pipi    '; 
  str.ID(iii,1)      = +13;   
  str.PQRtype(iii,1) = 'R';
  str.number(iii,1)  = 3;
  str.iso(iii,1:str.number(iii,1))      = R_isotope;
  str.bbb(iii,1:str.number(iii,1))      = R_pipi;       
  str.thelower(iii,1:str.number(iii,1)) = R_pipi_lower; 
  str.theupper(iii,1:str.number(iii,1)) = R_pipi_upper;
end

if (iii ~= length(IndicesToSearch))
  error('you have made a counting mistake in removeCO2lines')
end

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
for iii = 1:length(IndicesToSearch)
  zstr = str.name(iii,:);                     %type of band
  ID = str.ID(iii,1);                         %code tells which type of band
  bbb = str.bbb(iii,1:str.number(iii,1));     %which bands there actually are
  PQRtype = str.PQRtype(iii,1);               %P Q or R band
  isotope  = str.iso(iii,1:str.number(iii,1));      %isiotope
  thelower = str.thelower(iii,1:str.number(iii,1)); %lower state vib quant nums
  theupper = str.theupper(iii,1:str.number(iii,1)); %upper state vib quant nums

%%%%% if you only want ONE band, ONE branch, then
%%%%% 1) in run4co2quick.m, wher you have band=union(blah1,blah2), 
%%%%%     eventually only have            band=[720] or whatever
%%%%% 2) adjust the if condition with PQRtype as necessary : 
%%%%% if ((length(intersect(band,bbb)) == 1) & (low <= band) & (band  <= high)
%%%%%      ( (PQRtype == 'R'))

  if ((length(intersect(band,bbb)) == 1) & (low <= band) & (band  <= high))
    iFound = 1;          %a band has been found

    fname1 = [hitlin_fname];
    tempPQRlines = makeDAVEhitlin(band,vers,strengthM,homepath,...
                                  fname1,xfar,exchangelinecenters);
    
    num_band           = num_band+1;
    PQR(num_band)      = ID;
    bandtype(num_band) = band;

    whichone = find(band == bbb);
    v_l      = thelower(whichone);
    v_u      = theupper(whichone);
    isotype  = isotope(whichone);

    iso                    = lineORIG.iso';
    v_lower                = lineORIG.ilsgq';
    v_upper                = lineORIG.iusgq';
    j_lower_state_in_words = lineORIG.bslq;
    j_PQR                  = j_lower_state_in_words(:,5);
    j_level                = str2num(j_lower_state_in_words(:,6:8));

    index1 = find(j_PQR == PQRtype); 
    index2 = find(v_lower == v_l);   
    index3 = find(v_upper == v_u);   
    index4 = find(iso == isotype);   
    index  = intersect(intersect(intersect(index1,index2),index3),index4);
    if (length(index) > 0)
      fprintf(1,'looking to remove %4i lines from %s BAND %4i ... \n',...
              length(index),PQRtype,band)
    end

    tt = ['J numbers for ' num2str(band) ' PQR = ' num2str(PQR(num_band))];
    %plot(lineORIG.wnum(index),j_level(index)); title(tt); pause(0.1);

%%%%%%% run8co2_linemixUMBC.mhas the global settings of p2350_jmax r2350_jmax etc
%%%%%%% run8co2_linemixUMBC.mhas the global settings of p2350_jmax r2350_jmax etc
%%%%%%% run8co2_linemixUMBC.mhas the global settings of p2350_jmax r2350_jmax etc

%%%%%%%% if find_strongest > 0, find the 50 strongest lines from here!!!!!!!!!
    find_strongest = -1;
    if ((find_strongest > 0) & (length(index) >=  find_strongest))
      before = length(index);
      strengths_of_lines = lineORIG.stren(index);
      [yyy,iii] = sort(strengths_of_lines);
      iii = fliplr(iii);
      index = index(iii(1:find_strongest));
      after = length(index);
      fprintf(1,'a before = %3i after = %3i band = %3i \n',before,after,band);
    end

%%%%%%%% find the 60 lowest j lines for the P,R_sigsig 2350 !!!!!!!!!
    pband_iso = [2350];
    inside = length(intersect(band,pband_iso));
    if (PQR(num_band) < 0)
      find_lowest = p2350_jmax;
    elseif (PQR(num_band) > 0)
      find_lowest = r2350_jmax;
    end
    %%recall PQR(xxx) = -1 ==> P branch
    %%recall PQR(xxx) = 0  ==> Q branch
    %%recall PQR(xxx) = +1 ==> R branch
    if ((find_lowest > 0) & (inside > 0))
      jval = j_level(index);
      [yyy,iii] = sort(jval);   %from lowest to highest
      index = index(iii);       %this is now resorted from lo to hi
      jval = jval(iii);

      before = length(index);
      iii = find(jval <= find_lowest);  %%%find the lowest j values
      index = index(iii);
      after = length(index);
      fprintf(1,'b before = %3i after = %3i band = %3i \n',before,after,band);
    end

%%%%%%%% find the 60 lowest j lines for the PR_sigsig 2351 !!!!!!!!!
%%%%%%%% these are the strongest isotope bands and give problems
    pband_iso = [2351];
    inside = length(intersect(band,pband_iso));
    find_lowest = pr2351_jmax;
    %%recall PQR(xxx) = -1 ==> P branch
    %%recall PQR(xxx) = 0  ==> Q branch
    %%recall PQR(xxx) = +1 ==> R branch
    if ((find_lowest > 0) & (PQR(num_band) < 0) & (inside > 0))
      jval = j_level(index);
      [yyy,iii] = sort(jval);   %from lowest to highest
      index = index(iii);       %this is now resorted from lo to hi
      jval = jval(iii);

      before = length(index);
      iii = find(jval <=  find_lowest);  %%%find the lowest j values
      index = index(iii);
      after = length(index);
      fprintf(1,'c before = %3i after = %3i band = %3i \n',before,after,band);
    end

%%%%%%%% find the 60 lowest j lines for P_pipi 2321 or P_deltdelt 2311
%%%%%%%% these are the strongest isotope bands and give problems
    pband_iso = [2311 2321];
    inside = length(intersect(band,pband_iso));
    find_lowest = p2311_21_jmax;
    %%recall PQR(xxx) = -1 ==> P branch
    %%recall PQR(xxx) = 0  ==> Q branch
    %%recall PQR(xxx) = +1 ==> R branch
    if ((find_lowest > 0) & (PQR(num_band) < 0) & (inside > 0))
      jval = j_level(index);
      [yyy,iii] = sort(jval);   %from lowest to highest
      index = index(iii);       %this is now resorted from lo to hi
      jval = jval(iii);

      before = length(index);
      iii = find(jval <= find_lowest);  %%%find the lowest j values
      index = index(iii);
      after = length(index);
      fprintf(1,'d before = %3i after = %3i band = %3i \n',before,after,band);
    end

%%%%%%%%%%%%%%%%%%%%%%%% we are done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    doplot = quiet;
    if ((length(index) > 0)&(doplot>0))
      semilogy(lineORIG.wnum(index),lineORIG.stren(index),'+'); 
      title('these are the PQR lines');
      PQRtype
      band
      pause(0.01);
      fprintf(1,'number of lines = %3i\n',length(index));
    end

    %now that we have plucked out the lines, discard them!!!!!!!!!!
    %note this means we OVERWRITE index!!!!!!!!
    origlines  = 1:lineORIG.linct;
    linctindex = ~ismember(origlines,index);
    lengthless = lengthless + length(index);

    index       = find(linctindex == 1);

    line.linct   = sum(linctindex);
    line.igas   = lineORIG.igas; 
    line.iso    = lineORIG.iso(index);  
    line.wnum   = lineORIG.wnum(index); 
    line.stren  = lineORIG.stren(index); 
    line.tprob  = lineORIG.tprob(index); 
    line.abroad = lineORIG.abroad(index); 
    line.sbroad = lineORIG.sbroad(index); 
    line.els    = lineORIG.els(index); 
    line.abcoef = lineORIG.abcoef(index); 
    line.tsp    = lineORIG.tsp(index); 
    line.iusgq  = lineORIG.iusgq(index); 
    line.ilsgq  = lineORIG.ilsgq(index);
    line.gasid  = lineORIG.gasid(index);
    mmm         = 1:length(index);
    line.uslq   = lineORIG.uslq(index,:);
    line.bslq   = lineORIG.bslq(index,:); 
    line.ai     = lineORIG.ai(index,:); 
    line.ref    = lineORIG.ref(index,:);

    %reset lineORIG so it has "fewer" lines for the next check
    lineORIG = line;
  end
end
line = lineORIG;

%-------------------------------------------------------------------------
%%%else initialize line from lineORIG assuming nothing found
if (iFound == -1)
  fprintf(1,'nothing found in this bandset ......!! %4i \n',band);

%  index=1:lineORIG.linct;
%  line.linct   = lineORIG.linct; 
%  line.igas   = lineORIG.igas; 
%  line.iso    = lineORIG.iso(index);  
%  line.wnum   = lineORIG.wnum(index); 
%  line.stren  = lineORIG.stren(index); 
%  line.tprob  = lineORIG.tprob(index); 
%  line.abroad = lineORIG.abroad(index); 
%  line.sbroad = lineORIG.sbroad(index); 
%  line.els    = lineORIG.els(index); 
%  line.abcoef = lineORIG.abcoef(index); 
%  line.tsp    = lineORIG.tsp(index); 
%  line.iusgq  = lineORIG.iusgq(index); 
%  line.ilsgq  = lineORIG.ilsgq(index);
%  for mm=1:length(index)
%    line.ZUSLQ(mm,1:9)   = lineORIG.ZUSLQ(index(mm),1:9);
%    line.ZBSLQ(mm,1:9)   = lineORIG.ZBSLQ(index(mm),1:9); 
%    line.ZAI(mm,1:3)   = lineORIG.ZAI(index(mm),1:3); 
%    line.ZREF(mm,1:6)  = lineORIG.ZREF(index(mm),1:6);
%  end
%  line.gasid  = lineORIG.gasid(index);

  line=lineORIG;
end

%-------------------------------------------------------------------------
