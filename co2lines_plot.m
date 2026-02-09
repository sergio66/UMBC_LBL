%function [line1]=co2lines_plot(band,hitran_version);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% see HITRAN bd_vibs.for     for the spectroscopic notation  %%%%%%
%%%%%%            /asl/data/hitran/HITRAN2k/BD_Vibs.for           %%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% this file finds the lines you want for the band in question
%%% copied from co2lines.m

%% iHIT = input('enter HITRAN version (1992,1996,1998,2000,2004,2008,2012) : ');
iHIT = hitran_version;
if iHIT == 2000
  hitran_version = 'h2k';
else
  hitran_version = num2str(iHIT); 
  hitran_version = ['h' hitran_version(3:4)];
end

gasID = 2;
%%fnamePRE = '/asl/data/hitran/h92.by.gas/g';
fnamePRE = '/salsify/scratch4/h96.by.gas/g';
fnamePRE = '/asl/data/hitran/h98.by.gas/g';
fnamePRE = ['/asl/data/hitran/' hitran_version '.by.gas/g'];
fnamePRE = [hitranpath  hitran_version '.by.gas/g'];

fnamePOST = '.dat';
fnameIN = int2str(gasID);
hitlin_fname = [fnamePRE fnameIN fnamePOST];

start = 600;
stop = 2800;

%%line = hitread(start,stop,1.0e-28,gasID,hitlin_fname);
line = hitread(start,stop,0,gasID,hitlin_fname);

% band = 2350;
% band = 667;

if (band == 618)
  v_l=2;v_u=3;
elseif (band == 648)
  v_l=1;v_u=2;isotope=2;
elseif (band == 662)
  v_l=1;v_u=2;isotope=3;
elseif (band == 667)
  v_l=1;  v_u=2; 
elseif (band == 720)
  v_l=2;  v_u=5; 
elseif (band == 791)
  v_l=3;v_u=8;
elseif (band == 2080)
  v_l=1;v_u=8;

%find all PQR lines for isotope 1 : these are PQR_deltpi
elseif (band == 668)
  v_l=2;v_u=4;
elseif (band == 740)
  v_l=4;v_u=8;
elseif (band == 2093)
  v_l=2;v_u=14;

%find all PQR lines for isotope 1 : these are PQR_sigsig
elseif (band == 2350)
  v_l=1;v_u=9;
elseif (band == 2351)
  v_l=1;v_u=9; isotope=2;
elseif (band == 2352)
  v_l=1;v_u=9; isotope=3;
elseif (band == 2353)
  v_l=3;v_u=23; 
elseif (band == 2354)
  v_l=5;v_u=25; 

%find all PQR lines for isotope 1 : these are PQR_pipi
elseif (band == 2320) 
  v_l=2;v_u=16; 
elseif (band == 2321) 
  v_l=2;v_u=16; isotope=2;
elseif (band == 2322) 
  v_l=2;v_u=16; isotope=3;

%find all PQR lines for isotope 1 : these are PQR_deltdelt
elseif (band == 2310) 
  v_l=4;v_u=24;
elseif (band == 2311) 
  v_l=4;v_u=24; isotope=2;
  end

%[v1 v2 L v3 r] <----- [v1 v2 L v3 r]
%class 5    CO2 
%c                                       .......1   .......2   .......3    
%      DATA (AVIB(I),I=171,229)/        '   00001','   01101','   10002', 
%c 
%c      .......4   .......5   .......6   .......7   .......8   .......9 
%     +'   02201','   10001','   11102','   03301','   11101','   00011', 
%c 
%c      ......10   ......11   ......12   ......13   ......14   ......15 
%%     +'   20003','   12202','   20002','   04401','   12201','   20001', 
%c 
%c      ......16   ......17   ......18   ......19   ......20   ......21 
%     +'   01111','   21103','   13302','   21102','   05501','   13301', 
%c 
%c      ......22   ......23   ......24   ......25   ......26   ......27 
%     +'   21101','   10012','   02211','   10011','   30004','   22203', 
%c 
%c      ......28   ......29   ......30   ......31   ......32   ......33 
%     +'   14402','   30003','   22202','   06601','   30002','   14401', 
%c 
%c      ......34   ......35   ......36   ......37   ......38   ......39 
%     +'   22201','   30001','   11112','   03311','   11111','   00021', 
%c 
%c      ......40   ......41   ......42   ......43   ......44   ......45 
%     +'   31104','   31103','   31102','   20013','   12212','   23301', 
%c 
%c      ......46   ......47   ......48   ......49   ......50   ......51 
%     +'   31101','   04411','   20012','   12211','   20011','   01121', 
%c 
%c      ......52   ......53   ......54   ......55   ......56   ......57 
%     +'   40004','   32203','   21113','   40002','   13312','   05511', 
%c 
%c      ......58   ......59 
%     +'   21112','   13311'/ 

%get the branch info P Q or R
pqr=line.bslq(:,5);

isotope = 1;
v_l = 1;

ind=find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
      (line.iso == isotope) & ...
      ((pqr' == 'P') | (pqr' == 'Q') | (pqr' == 'R'))  );
ii=length(ind);
fprintf(1,'number of PQR lines selected = %5i \n',ii);

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

semilogy(line1.wnum,line1.stren);
bar(line1.wnum,line1.stren,2);

%%for 2350 band
axis([2275 2400 0 4e-18]); xlabel('Wavenumber cm-1'); ylabel('LineStrength');
text(2320,3.7e-18,'P branch','Fontsize',15)
text(2360,3.7e-18,'R branch','Fontsize',15)
title('(\nu_{1} \nu_{2} \nu_{3}) = (000) \leftarrow (001)')
printfig(1,'/home/sergio/PAPERS/SPECIALTALKS/JOB06/FIGS/v3','jpg')

%%for 667 band
axis([625 725 0 3e-19]); xlabel('Wavenumber cm-1'); ylabel('LineStrength');
text(640,  2e-19,'P branch','Fontsize',15)
text(663.5,2.6e-19,'Q branch','Fontsize',15)
text(680,  2e-19,'R branch','Fontsize',15)
title('(\nu_{1} \nu_{2} \nu_{3}) = (000) \leftarrow (010)')
printfig(1,'/home/sergio/PAPERS/SPECIALTALKS/JOB06/FIGS/v2','jpg')

iDo = -1;
if iDo > 0
  accuracy    = line.ai(ind,1:3);
  dipole      = line.tprob(ind)';    %transition probablility  
  elower      = line.els(ind)';
  freq        = line.wnum(ind)';
  gas_id      = gasID;
  iso         = line.iso(ind)';
  j_lower     = line.bslq(ind,1:9);
  j_upper     = line.uslq(ind,1:9);
  p_shift     = line.tsp(ind)';
  reference   = line.ref(ind,1:6);
  stren       = line.stren(ind)';
  v_lower     = line.ilsgq(ind)';
  v_upper     = line.iusgq(ind)';
  w           = line.abroad(ind)';  %air broadenend widths
  w_s         = line.sbroad(ind)';  %self broadened widths
  w_temp      = line.abcoef(ind)';  %tempr correct for air broadened widths
  end  

fprintf(1,'lower quanta = %3i upper quanta = %3i \n',v_l,v_u);
fprintf(1,'lower bound = %10.5f upper bound = %10.5f \n',min(freq),max(freq));
