function [line1,line1A] = co2lines(band,hitran_version);

%%%this file finds the lines you want for the band in question

if nargin == 1
  hitran_version = 'h12';
  hitran_version = 'h16';
  hitran_version = 'h20';
  hitran_version = 'h24';    
end

gasID=2;

fnamePRE     = '/asl/data/hitran/h92.by.gas/g';
fnamePRE     = '/salsify/scratch4/h96.by.gas/g';
fnamePRE     = '/asl/data/hitran/h98.by.gas/g';
fnamePRE     = '/asl/data/hitran/h08.by.gas/g';
fnamePRE     = ['/asl/data/hitran/' hitran_version '.by.gas/g'];
fnamePRE     = [hitranpath  hitran_version '.by.gas/g'];

fnamePOST    = '.dat';
fnameIN      = int2str(gasID);
hitlin_fname = [fnamePRE fnameIN fnamePOST];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assume isotope ======= 1 unless stated
isotope=1;

%assume we have to make a hit* file

do_makefile = 1;

%find all PQR lines for isotope 1 (or 2 or 3) : these are PQR_sigpi
if (band == 618)
  v_l = 2;v_u = 3;
elseif (band == 648)
  v_l = 1;v_u = 2;isotope = 2;
elseif (band == 662)
  v_l = 1;v_u = 2;isotope = 3;
elseif (band  ==  667)
  v_l = 1;  v_u = 2; 
elseif (band == 720)
  v_l = 2;  v_u = 5; 
elseif (band == 791)
  v_l = 3;v_u = 8;
elseif (band == 1932)
  v_l = 1; v_u = 6;  
elseif (band == 2080)
  v_l = 1;v_u = 8;
elseif (band == 2129)
  v_l = 2; v_u = 15;  

%find all PQR lines for isotope 1 : these are PQR_deltpi
elseif (band == 668)
  v_l = 2;v_u = 4;
elseif (band == 740)
  v_l = 4;v_u = 8;
elseif (band == 2093)
  v_l = 2;v_u = 14;

%find all PQR lines for isotope 1 : these are PQR_sigsig
elseif (band == 2350)
  v_l = 1;v_u = 9;
elseif (band == 2351)
  v_l = 1;v_u = 9; isotope = 2;
elseif (band == 2352)
  v_l = 1;v_u = 9; isotope = 3;
elseif (band == 2353)
  v_l = 3;v_u = 23; 
elseif (band == 2354)
  v_l = 5;v_u = 25; 
%elseif (band == 2355)  isotope = 4;
%  v_l = 1; v_u = 9; 

%find all PQR lines for isotope 1 : these are PQR_pipi
elseif (band == 2320) 
  v_l = 2;v_u = 16; 
elseif (band == 2321) 
  v_l = 2;v_u = 16; isotope = 2;
elseif (band == 2322) 
  v_l = 2;v_u = 16; isotope = 3;

%find all PQR lines for isotope 1 : these are PQR_deltdelt
elseif (band == 2310) 
  v_l = 4;v_u = 24;
elseif (band == 2311) 
  v_l = 4;v_u = 24; isotope = 2;

else
  fprintf('no need to generate hitBLAH.mat file!!');
  do_makefile =  -1;
end

if (do_makefile > 0)
  %%%%%%some of these CO2 bands go quite far!!!!!! eg P2350 is from 2160-2350
  %call hittomat to get the lines

  stretch = 500;
  start = band-stretch;  
  stop = band+stretch;

  %%line = hitread(start,stop,1.0e-28,gasID,hitlin_fname);
  line = hitread(start,stop,0,gasID,hitlin_fname,-1);

  %jj = find(line.wnum <= 2183.3 & line.wnum >= 2183.0); 
  %jj = find(line.wnum <= 2186.7 & line.wnum >= 2186.5); 
  %jj = find(line.wnum <= 2190.05 & line.wnum >= 2189.95); 
  %line.wnum(jj);
  %line.tsp(jj);

  line = should_I_translate2oldHITparams(line,2,hitran_version);

  %get the branch info P Q or R
  pqr = line.bslq(:,5);

  ind = find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
      (line.iso == isotope) & ...
      ((pqr' == 'P') | (pqr' == 'Q') | (pqr' == 'R'))  );
  ii = length(ind);
  fprintf(1,'number of PQR lines selected for band %5i = %5i \n',band,ii);

  line1.ai    = line.ai(ind,1:3);
  line1.trpob = line.tprob(ind)';    %transition probablility  
  line1.els   = line.els(ind)';
  line1.wnum  = line.wnum(ind)';
  line1.iso   = line.iso(ind)';
  line1.bslq  = line.bslq(ind,1:9);
  line1.uslq  = line.uslq(ind,1:9);
  line1.stren = line.stren(ind)';
  line1.ilsgq = line.ilsgq(ind)';
  line1.iusgq = line.iusgq(ind)';

  line1A.accuracy    = line.ai(ind,1:3);
  line1A.dipole      = line.tprob(ind)';    %transition probablility  
  line1A.elower      = line.els(ind)';
  line1A.freq        = line.wnum(ind)';
  line1A.gas_id      = gasID;
  line1A.iso         = line.iso(ind)';
  line1A.j_lower     = line.bslq(ind,1:9);
  line1A.j_upper     = line.uslq(ind,1:9);
  line1A.p_shift     = line.tsp(ind)';
  line1A.reference   = line.ref(ind,1:6);
  line1A.stren       = line.stren(ind)';
  line1A.v_lower     = line.ilsgq(ind)';
  line1A.v_upper     = line.iusgq(ind)';
  line1A.w           = line.abroad(ind)';  %air broadenend widths
  line1A.w_s         = line.sbroad(ind)';  %self broadened widths
  line1A.w_temp      = line.abcoef(ind)';  %tempr correct for air broadened widths
end  

fprintf(1,'lower quanta = %3i upper quanta = %3i \n',v_l,v_u);
fprintf(1,'lower bound = %10.5f upper bound = %10.5f \n',min(line1A.freq),max(line1A.freq));
