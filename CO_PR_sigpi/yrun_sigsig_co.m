function [f,voi,fc,theRATIO]=yrun(temperature,band,ptotal,pself,layeramt,...
     prb,LVF,IO,birn,boxcar,hiResFreq,outwave,theRATIO)

%[f,y]=yrun(temperature,band,ptotal,pself,layeramt,frequser,prb,LVF,IO,birn)
%this function computes CO2 R-line mixing due to the 2350 band
%disp('bands: 2350')

if ((band ~= 2150))
  error('need 2150 CO band!!!!!!!!!')
  end

if (prb == 'p')
  prb='P';
  end
if (prb == 'r')
  prb='R';
  end

MGC=8.314674269981136; %%%%%%%%%%% correct value 
%layer amt in kmoles/cm^2
%density = n/V = partpress/(RT)
%density * path length = layer amt
%path length = layer amt/ density
path_length=layeramt*1000*10000;      %change from kilomoles/cm2 to moles/m2
den=101325*pself/MGC/temperature;    %density in no of moles/m^3
path_length=path_length/den*100;      %path length changed from m to cm

%for testing data in /salsify/data/Strow/Tobin/Co2q_B_sigpie 
%dont forget to change pressure from torr to atm (p/760) 
%path_length=layeramt; 

[jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr,stuff]=...
    loader_co(temperature,band,path_length,ptotal,pself,prb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this is testing LOOP (f77) vs XLOOP (MATLAB, which calls f77 code anyway)

addpath /home/sergio/SPECTRA
addpath /home/sergio/SPECTRA/FORTRANLINUX
addpath /home/sergio/MATLABCODE

%junkvect = loop(very_near.iso,mass_iso,brd,strength,freq,fine_freq,...
%                   tempr,very_near.linct,length(fine_freq),lvgNum);
very_near_wow = load('/home/sergio/SPECTRA/very_near_co.mat');
lvgNum = 1;
junkvect = loop(very_near_wow.very_near.iso,very_near_wow.mass_iso,very_near_wow.brd,...
                very_near_wow.strength,very_near_wow.freq,very_near_wow.fine_freq,...
                very_near_wow.tempr,very_near_wow.very_near.linct,length(very_near_wow.fine_freq),lvgNum);
xjunkvect = xloop(very_near_wow.very_near.iso,very_near_wow.mass_iso,very_near_wow.brd,...
                very_near_wow.strength,very_near_wow.freq,very_near_wow.fine_freq,...
                very_near_wow.tempr,very_near_wow.very_near.linct,length(very_near_wow.fine_freq),lvgNum);
sum(junkvect-xjunkvect)

mass_amu = 12+16; %% CO
mass_amu = 27.99491;

%brd      = broad(p,ps,1.0,for_brd,self_brd,pwr,tempr,gasID);
brd      = broad(ptotal,pself,1.0,w_forr,w_selfr,0.5,temperature,5);
%{
bah = very_near_wow.very_near.iso;  whos bah;   bah = 1.0d0*ones(size(jr)); whos bah
bah = very_near_wow.mass_iso;       whos bah;   bah = 1.0d0*[mass_amu mass_amu+1]; whos bah
bah = very_near_wow.brd;            whos bah;   bah = brd; whos bah
bah = very_near_wow.strength;       whos bah;   bah = strenrt; whos bah
bah = very_near_wow.freq;           whos bah;   bah = freqr; whos bah
bah = very_near_wow.fine_freq;      whos bah;   bah = outwave; whos bah
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is lorentz
%% xoutvect  = xloop(1.0d0*ones(size(jr')),1.0d0*[mass_amu mass_amu+1],brd',strenrt',freqr',...
%%                  outwave,temperature,length(strenrt),length(outwave),-1);

outvect = loop(1.0d0*ones(size(jr')),1.0d0*[mass_amu mass_amu+1],brd',strenrt',freqr',...
                 outwave,temperature,length(strenrt),length(outwave),+1);
xoutvect = xloop(1.0d0*ones(size(jr')),1.0d0*[mass_amu mass_amu+1],brd',strenrt',freqr',...
                 outwave,temperature,length(strenrt),length(outwave),+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
 input linecenter cm-1, line strength (cm2/molecule)
 input mass(amu) T(K)
 input SpeedAvgWidth(Lorentz) SpeedDependenceWidth cm-1
 input SpeedAvgShift SpeedDependaceShift cm-1
 input VelChangeFreq cm-1 and CorrelationParam
 input Sig0,SigF,dSig
%}

frq_dependent  = [];
spd_dependentR = [];
spd_dependentI = [];
homexdir = pwd;
cd /home/sergio/SPECTRA/FORTRANLINUX/SpeedDependentVoigt/
fipspeed = mktemp('fipspd.txt');
fopspeed = mktemp('fopspd.txt');

tic
for iix = 1 : length(jr)
  fprintf(1,'line %4i of %4i \n',iix,length(jr));
  fid = fopen(fipspeed,'w');
  fprintf(fid,' %10.6f %10.6e \n',freqr(iix),strenrt(iix));
  fprintf(fid,' %10.6f %10.6e \n',mass_amu*1.0d0,temperature);
  fprintf(fid,' %10.6f %10.6e \n',brd(iix),0.0);
  fprintf(fid,' %10.6f %10.6e \n',0.0,0.0);
  fprintf(fid,' %10.6f %10.6e \n',0.0,0.0);  
  fprintf(fid,' %10.6f %10.6f %10.6e \n',min(outwave),max(ceil(outwave)),mean(abs(diff(outwave))));  
  fclose(fid);
  spder = ['!driver_HTP.x < ' fipspeed ' > ' fopspeed];
  eval(spder);
  sedder = ['!sed -i -e 1,11d ' fopspeed];
  eval(sedder)
  gah = load(fopspeed);
  if iix == 1
    frq_dependent = gah(:,2);
    spd_dependentR = gah(:,3);
    spd_dependentI = gah(:,4);    
  else
    spd_dependentR = spd_dependentR + gah(:,3);
    spd_dependentI = spd_dependentI + gah(:,4);    %% really .. you will need + Y(iix)*gah(:,3);
  end
end
toc

rmer = ['!/bin/rm ' fipspeed ' ' fopspeed]; eval(rmer)
cder = ['cd ' homexdir]; eval(cder);
spd_dependentR = spd_dependentR';
spd_dependentI = spd_dependentI';
frq_dependent  = frq_dependent';

whos outvect xoutvect spd_dependentR spd_dependentI outwave frq_dependent

plot(outwave,outvect,'b',outwave,xoutvect,'r',frq_dependent,spd_dependentR,'k.-');
plot(outwave,outvect./xoutvect,outwave(1:end-2),outvect(1:end-2)./spd_dependentR(3:length(outwave)));
plot(outwave,outvect,'b.-',outwave,xoutvect,'r',outwave(1:end-2),spd_dependentR(3:length(outwave)),'k.-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
 input linecenter cm-1, line strength (cm2/molecule)
 input mass(amu) T(K)
 input SpeedAvgWidth(Lorentz) SpeedDependenceWidth cm-1
 input SpeedAvgShift SpeedDependaceShift cm-1
 input VelChangeFreq cm-1 and CorrelationParam
 input Sig0,SigF,dSig
%}

zfrq_dependent  = [];
zspd_dependentR = [];
zspd_dependentI = [];
homexdir = pwd;
cd /home/sergio/SPECTRA/FORTRANLINUX/SpeedDependentVoigt/
fipspeed = mktemp('fipspd.txt');
fopspeed = mktemp('fopspd.txt');

tic
for iix = 1 : length(jr)
  fprintf(1,'line %4i of %4i \n',iix,length(jr));
  fid = fopen(fipspeed,'w');
  fprintf(fid,' %10.6f %10.6e \n',freqr(iix),strenrt(iix));
  fprintf(fid,' %10.6f %10.6e \n',mass_amu*1.0d0,temperature);
  fprintf(fid,' %10.6f %10.6e \n',brd(iix),brd(iix)/2.0);
  fprintf(fid,' %10.6f %10.6e \n',0.0,0.0);
  fprintf(fid,' %10.6f %10.6e \n',0.0,0.0);  
  fprintf(fid,' %10.6f %10.6f %10.6e \n',min(outwave),max(ceil(outwave)),mean(abs(diff(outwave))));  
  fclose(fid);
  spder = ['!driver_HTP.x < ' fipspeed ' > ' fopspeed];
  eval(spder);
  sedder = ['!sed -i -e 1,11d ' fopspeed];
  eval(sedder)
  gah = load(fopspeed);
  if iix == 1
    zfrq_dependent = gah(:,2);
    zspd_dependentR = gah(:,3);
    zspd_dependentI = gah(:,4);    
  else
    zspd_dependentR = zspd_dependentR + gah(:,3);
    zspd_dependentI = zspd_dependentI + gah(:,4);    %% really .. you will need + Y(iix)*gah(:,3);
  end
end
toc
						      
cder = ['cd ' homexdir]; eval(cder);
zspd_dependentR = zspd_dependentR';
zspd_dependentI = zspd_dependentI';
zfrq_dependent  = zfrq_dependent';

plot(outwave,outvect,'b',outwave,xoutvect,'r',frq_dependent,spd_dependentR,'k.-');
plot(outwave,outvect./xoutvect,outwave(1:end-2),outvect(1:end-2)./spd_dependentR(3:length(outwave)));
plot(outwave,outvect,'b.-',outwave,xoutvect,'r',...
     outwave(1:end-2),spd_dependentR(3:length(outwave)),'k.-',outwave(1:end-2),zspd_dependentR(3:length(outwave)),'g*-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

keyboard_nowindow

z2frq_dependent  = [];
z2spd_dependentR = [];
z2spd_dependentI = [];
homexdir = pwd;
cd /home/sergio/SPECTRA/FORTRANLINUX/SpeedDependentVoigt/
fipspeed = mktemp('fipspd.txt');
fopspeed = mktemp('fopspd.txt');

%error('driver_HTP_all.x not yet done, this allows all lines to be read in ....')
tic
fid = fopen(fipspeed,'w');
fprintf(fid,' %4i  %10.6f %10.6e \n',length(jr),mass_amu*1.0d0,temperature);
fprintf(fid,' %10.6f %10.6f %10.6e \n',min(outwave),max(ceil(outwave)),mean(abs(diff(outwave))));  
for iix = 1 : length(jr)
  linepars = [freqr(iix) strenrt(iix) brd(iix) brd(iix)/2.0 0.0 0.0 0.0 0.0];
  fprintf(fid,' %10.6f %10.6e %10.6f %10.6e %10.6f %10.6e %10.6f %10.6e \n',linepars);
end  
fclose(fid);
spder = ['!driver_HTP_all.x < ' fipspeed ' > ' fopspeed];
eval(spder);
toc

sedder = ['!sed -i -e 1,12d ' fopspeed];
eval(sedder)
gah = load(fopspeed);
z2frq_dependent = gah(:,2);
z2spd_dependentR = gah(:,3);
z2spd_dependentI = gah(:,4);    
rmer = ['!/bin/rm ' fipspeed ' ' fopspeed]; eval(rmer)
cder = ['cd ' homexdir]; eval(cder);
z2spd_dependentR = z2spd_dependentR';
z2spd_dependentI = z2spd_dependentI';
z2frq_dependent  = z2frq_dependent';

plot(outwave,outvect,'b',outwave,xoutvect,'r',frq_dependent,spd_dependentR,'k.-');
plot(outwave,outvect./xoutvect,outwave(1:end-2),outvect(1:end-2)./spd_dependentR(3:length(outwave)));
plot(outwave,outvect,'b.-',outwave,xoutvect,'r',...
     outwave(1:end-2),spd_dependentR(3:length(outwave)),'k.-',outwave(1:end-2),zspd_dependentR(3:length(outwave)),'g*-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fudger_and_threeregions;

elower        = elowerr;
[elower,jall] = efitter_co(jr,band,elowerr,elower,prb); 

[f,voi,fc,theRATIO] = ...  
    driver4um(temperature,ptotal,pself,layeramt,...   %gas layer params  
                  band,prb,LVF,IO,birn,...               %spectra computation  
                  boxcar,hiResFreq,outwave,theRATIO,...  %boxcar, wavenumbers  
                  path_length,elower,jall,...            %from initial call  
                  jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr,stuff); 
