%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  tdif_birn.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file computes the difference between 1st order mix transmission
%       with birnbaum and actual transmission.
%
%  same as tdif_birn_Oct99 except that it gets rid of the magic numbers
%  in the Cousin-Mixing blending
%  
%  also for some reason the Lorentz lineshape in doVmix.f is not working
%  so call the Voigt lineshape
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function final_diff=tdif_birn_goodparams(bstart,freqr,freqp,jp,f,flag4)

global beta beta_pure beta_for bsm duration frequency_shift fudge
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor
global trans_amplr populationr population_tr t_rawdata frequency_shift
global trans_amplp populationp population_tp strenpt w_forp w_selfp
global voi_back pick_no
global K_back strenrt w_forr w_selfr abscoef 
global sigsig_comp sigsigR sigsigP ymixR wwwR
global c0 cspan

figure(1)

beta = bstart(1);
duration = bstart(2);
frequency_shift = bstart(3);
c0 = bstart(4);
cspan = bstart(5);

band_strength_multiplier=flag4;
tau2=duration;
%%fprintf(1,'at beginning beta = %8.6f dur = %8.6f \n',beta,duration);

load strenths
for i=1:length(temperature)
  eval(['global W_plus_r' num2str(i)])
  eval(['global W_plus_p' num2str(i)])
  end

diff=[];

for iii=1:length(temperature)
  ind=find(f(:,iii)~=0);
  fnow=f(ind,iii);
  trd=t_rawdata(ind,iii);
  kback=K_back(ind,iii);
  
  subplot(211); plot(fnow,trd);
  subplot(212); plot(fnow,exp(-kback));
  title('input data and backgnd'); pause(0.1)

  K_scale_lor=density*pressure_self(iii)/pressure_ref*...
		temperature_ref*path_length(iii)/temperature(iii)/pi;
  abscoefP=zeros(length(fnow),1);  
  abscoefR=zeros(length(fnow),1);  
  abscoef=zeros(length(fnow),1);  

  %%%%%%%%%%%%% R-branch lines %%%%%%%%%%%%%%%%%%%%%%%%%%%
  freqr_shift=freqr+frequency_shift/100;
  no_lines=length(freqr_shift);
  Y_1stmix=zeros(1,no_lines);
  ZZ=(beta+(1-beta)*eye(no_lines,no_lines));
  eval(['W_plus2=W_plus_r' num2str(iii) '.*ZZ;']);clear ZZ

  ymix=y1sNEW(0,freqr_shift,0,0,W_plus2,trans_amplr(:,iii),0,beta);
  ymixR=ymix;

  %%ymix = zeros(size(ymix));
  %%ymixR = zeros(size(ymixR));

  w_tot=pressure_self(iii)*w_selfr(:,iii)+pressure_for(iii)*... 
                w_forr(:,iii)/pressure_ref; 

  clf
  plot(ymix); title('ymix'); pause(0.1)

  voivoi=v_lormix(iii,freqr,fnow',ymix,temperature(iii),w_forr(:,iii),... 
            w_selfr(:,iii),strenrt(:,iii),'V','b'); 
  clf
  plot(fnow,voivoi); title('try 1a'); pause(0.1)

  pah=find(voivoi < 0);
  voivoi(pah)=0.0;

  abscoefR=voivoi;
  clf
  plot(fnow,voivoi); title('try 1b'); pause(0.1)

  %new !!!!!!!! blend in cousin!!!!!!!!!!
  lorcousin=abscoefR;

  %%%%look at blender4um.m
  magicR=max(freqr_shift) + 40; 
  magicR=magicR-20.0;
  magicRR=magicR+20;

  disp('magicR and magicRR in tdif_birn')
  %[magicR magicRR]

  indblend=find((fnow > 0) & ((fnow >= magicR) & (fnow <= magicRR)));
  indcous=find((fnow > 0) & (fnow > magicRR));

  ymix0=zeros(size(freqr_shift));  
  voivoi=vcousin(iii,freqr,fnow',ymix0,temperature(iii),w_forr(:,iii),... 
            w_selfr(:,iii),strenrt(:,iii),'V','c',1.0); 

  lorcousin(indcous)=voivoi(indcous);    %no blending here ... just cous

  ffff=fnow(ind);

  %these are the start and stop points of current wavevector
  s1=ffff(1);   s2=ffff(length(ffff));
  %endfull tells where the end of full mixing is, and blend is the
  %width of the blending region
  %%%%%[endfull,blend]=blender4um(band,prb,freqr,jr,kk);
  endfull=magicR;
  blend=20.0;

  %%%%%%%%%%%%%%   FULL       |     BLEND        |       COUS
  %%%%%%%%%%%%%               fa                fb
  fa=endfull; fb=endfull+blend;          
  voi=(fb-ffff(indblend)')/(fb-fa).*abscoefR(indblend) + ...
       (ffff(indblend)'-fa)/(fb-fa).*voivoi(indblend);

  lorcousin(indblend)=voi;            %blending here

  abscoefR=lorcousin;
  sigsigR=abscoefR;

  %%%%%%%%%%%%% P-branch lines %%%%%%%%%%%%%%%%%%%%%%%%%%%
  freqp_shift=freqp+frequency_shift/100;
  no_lines=length(freqp_shift);  
  Y_1stmix=zeros(1,no_lines);
  ZZ=(beta+(1-beta)*eye(no_lines,no_lines));
  eval(['W_plus2=W_plus_p' num2str(iii) '.*ZZ;']);clear ZZ

  ymix=y1sNEW(0,freqp_shift,0,0,W_plus2,trans_amplp(:,iii),0,beta);

  w_tot=pressure_self(iii)*w_selfp(:,iii)+pressure_for(iii)*... 
                w_forp(:,iii)/pressure_ref; 

  voivoi=v_lormix(iii,freqp,fnow',ymix,temperature(iii),w_forp(:,iii),... 
            w_selfp(:,iii),strenpt(:,iii),'V','b'); 

  abscoefP=voivoi;

  lorcousin=abscoefP;

  magicP=max(freqp_shift)+40; 

  magicP=max(freqp_shift) + 40; 
  magicP=magicP-20.0;

  ind=find((fnow > 0) & (fnow > magicP));

  ymix0=zeros(size(freqp_shift));  
  voivoi=vcousin(iii,freqp,fnow',ymix0,temperature(iii),w_forp(:,iii),... 
            w_selfp(:,iii),strenpt(:,iii),'V','c',1.0); 
  lorcousin(ind)=voivoi(ind);

  abscoefP=lorcousin;
  sigsigP=abscoefP;

  %%%%%%%%%%%%% add P and R branch lines %%%%%%%%%%%%%%%%%%%%%%%%%%%

  abscoef=abscoefP+abscoefR;

  subplot(211); 
  plot(fnow,abscoefP,fnow,abscoefR,'r'); 
  subplot(212);
  plot(fnow,exp(-abscoefP),fnow,exp(-abscoefR),'r'); 
  title('sigsig P and R'); pause(0.1)

  abscoef=abscoef';
  sigsig_comp=abscoef;

  subplot(211); 
  plot(fnow,abscoef); 
  subplot(212);
  plot(fnow,kback); 
  title('P+R vs backgnd'); pause(0.1)

  abscoef=abscoef*bsm+kback*bsm;
  error=trd-exp(-abscoef);

  %figure(iii);clg;
  %plot(fnow,trd,'r',fnow,exp(-abscoef),'g',fnow,error,'b');grid;drawnow

  figure(2)	
  subplot(length(temperature),1,iii)
  plot(fnow,trd,'b',fnow,exp(-abscoef),'r',fnow,error,'g');grid;drawnow
  legend('data','sims','error');
  axis([min(fnow) max(fnow) -0.2 1.2])
  pause(1.0)
  fss=frequency_shift;
%eval(['save fileGL' num2str(pick_no(iii)) ' fnow trd abscoef tau2 beta frequency_shift'])
%eval(['save fileIND' num2str(pick_no(iii)) ' fnow trd abscoef tau2 beta frequency_shift'])
eval(['save file' num2str(pick_no(iii)) ' fnow trd abscoef tau2 beta fss kback sigsigR sigsigP fudge c0 cspan'])

  diff=[diff;error];
  end

nbeta=num2str(beta);nbsm=num2str(bsm);
ndur=num2str(duration);nfs=num2str(frequency_shift);
disp(['bsm=' nbsm '   beta=' nbeta '   tau2=' ndur '   fs=' nfs])
fprintf(1,'beta= %10.8f dofc = %10.8f \n',beta,duration);

final_diff=diff;

