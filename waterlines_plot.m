function [line1] = waterlines_plot(wn1,wn2);

addpath /home/sergio/SPECTRA/read_hitr06

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% see HITRAN bd_vibs.for     for the spectroscopic notation  %%%%%%
%%%%%%            /asl/data/hitran/HITRAN2k/BD_Vibs.for           %%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% this file finds the lines you want for the band in question
%%% copied from co2lines.m

iHIT = input('enter HITRAN version (1992,1996,1998,2000,2004,2008,2012,2016) : ');
if iHIT == 2000
  strHIT = 'h2k';
else
  strHIT = num2str(iHIT); 
  strHIT = ['h' strHIT(3:4)];
end

gasID=1;
%%fnamePRE='/asl/data/hitran/h92.by.gas/g';
fnamePRE = '/salsify/scratch4/h96.by.gas/g';
fnamePRE = '/asl/data/hitran/h98.by.gas/g';
fnamePRE = ['/asl/data/hitran/' strHIT '.by.gas/g'];
fnamePOST ='.dat';
fnameIN = int2str(gasID);
hitlin_fname = [fnamePRE fnameIN fnamePOST];

if nargin <= 1
  wn1 = 0;
  wn2 = 4000;
end

start = wn1;
stop  = wn2;

%%line = hitread(start,stop,1.0e-28,gasID,hitlin_fname);
line = hitread(start,stop,0,gasID,hitlin_fname,-1);

figure(1); clf; semilogy(line.wnum,line.stren); title('ALL lines');

woo = find(line.iso == 4);
figure(2); clf; semilogy(line.wnum(woo),line.stren(woo),'.','markersize',5); title('HDO lines');

fprintf(1,'found %8i total lines and %8i HDO lines \n',length(line.wnum),length(woo));

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

%fprintf(1,'lower quanta = %3i upper quanta = %3i \n',v_l,v_u);
%fprintf(1,'lower bound = %10.5f upper bound = %10.5f \n',min(freq),max(freq));
