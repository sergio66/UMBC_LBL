%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  tdif_birn.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file computes the difference between 1st order mixing transmission
%       with birnbaum and actual transmission.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function final_diff=tdif_birn(bstart,freqr,freqp,jp,f,flag4)

global beta beta_pure beta_for bsm duration frequency_shift fudge
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor
global trans_amplr populationr population_tr t_rawdata frequency_shift
global trans_amplp populationp population_tp strenpt w_forp w_selfp
global voi_back pick_no
global K_back strenrt w_forr w_selfr abscoef 
global sigsig_comp sigsigR sigsigP ymixR wwwR

load strenths
fprintf(1,'at beginning beta = %8.6f dur = %8.6f \n',beta,duration);

for i=1:length(temperature)
  eval(['global W_plus_r' num2str(i)])
  eval(['global W_plus_p' num2str(i)])
  end

i=0;
if flag4(1)==1
  i=i+1;bsm=bstart(i);
  end
if flag4(2)==1
  i=i+1;beta=bstart(i);
  end
if flag4(3)==1
  i=i+1;duration=bstart(i);
  end
if flag4(4)==1
  i=i+1;frequency_shift=bstart(i);
  end

band_strength_multiplier=bsm;

if (beta > 0)
  beta=min([ beta  1]);
else
  beta=0.95;
  end

if (duration < 0)
  duration=5.0e-3;
  end
if (duration > 1.05e-2)
  duration=5e-3;
  end

tau2=duration;

diff=[];

for iii=1:length(temperature)
  ind=find(f(:,iii)~=0);
  fnow=f(ind,iii);
  trd=t_rawdata(ind,iii);
  kback=K_back(ind,iii);

  K_scale_lor=density*pressure_self(iii)/pressure_ref*...
		temperature_ref*path_length(iii)/temperature(iii)/pi;
  abscoefP=zeros(length(fnow),1);  
  abscoefR=zeros(length(fnow),1);  
  abscoef=zeros(length(fnow),1);  

  %%%%%%%%%%%%% R-branch lines %%%%%%%%%%%%%%%%%%%%%%%%%%%
%  fprintf(1,'the pressures are \n');
%  [pressure_self(iii) pressure_for(iii)]
%  fprintf(1,'the strengths are \n');
%  strenrt(:,iii)
%  fprintf(1,'the w forr are \n');
%  w_forr(:,iii)
%  fprintf(1,'the w selfr are \n');
%  w_selfr(:,iii)

  freqr_shift=freqr+frequency_shift/100;
  no_lines=length(freqr_shift);
  Y_1stmix=zeros(1,no_lines);
  ZZ=(beta+(1-beta)*eye(no_lines,no_lines));
  eval(['W_plus2=W_plus_r' num2str(iii) '.*ZZ;']);clear ZZ

  ymix=y1sNEW(0,freqr_shift,0,0,W_plus2,trans_amplr(:,iii),0,beta);
  ymixR=ymix;

  w_tot=pressure_self(iii)*w_selfr(:,iii)+pressure_for(iii)*... 
                w_forr(:,iii)/pressure_ref; 

  voivoi=v_lormix(iii,freqr,fnow',ymix,temperature(iii),w_forr(:,iii),... 
            w_selfr(:,iii),strenrt(:,iii),'V','b'); 

  pah=find(voivoi < 0);
  voivoi(pah)=0.0;

  abscoefR=voivoi;

%new !!!!!!!! blend in cousin!!!!!!!!!!
  lorcousin=abscoefR;
  ind=find((fnow > 0) & (fnow > 2417)); 
  ymix0=zeros(size(freqr_shift));  
  voivoi=vcousin(iii,freqr,fnow',ymix0,temperature(iii),w_forr(:,iii),... 
            w_selfr(:,iii),strenrt(:,iii),'V','c',1.0); 
  ffff=fnow(ind);
  slpe=(ffff(length(ffff))-ffff(1)); 
  s1=ffff(1);   s2=ffff(length(ffff)); 

  slpe=(2437-ffff(1)); 
  s1=ffff(1);   s2=2437; 

  voi=((voivoi(ind).*(ffff'-s1) + abscoefR(ind).*(s2-ffff')))/slpe; 
  lorcousin(ind)=voi;

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
  ind=find((fnow > 0) & (fnow > max(freqp_shift)+40)); 
  ymix0=zeros(size(freqp_shift));  
  voivoi=vcousin(iii,freqp,fnow',ymix0,temperature(iii),w_forp(:,iii),... 
            w_selfp(:,iii),strenpt(:,iii),'V','c',1.0); 
  lorcousin(ind)=voivoi(ind);

  abscoefP=lorcousin;
  sigsigP=abscoefP;

  %%%%%%%%%%%%% add P and R branch lines %%%%%%%%%%%%%%%%%%%%%%%%%%%

  abscoef=abscoefP+abscoefR;

  plot(fnow,exp(-abscoefP),fnow,exp(-abscoefR)); 
  title('sigsig P and R'); pause(0.1)

  abscoef=abscoef';
  sigsig_comp=abscoef;

  abscoef=abscoef*bsm+kback*bsm;
  error=trd-exp(-abscoef);

  %figure(iii);clg;
  %plot(fnow,trd,'r',fnow,exp(-abscoef),'g',fnow,error,'b');grid;drawnow
	
  subplot(length(temperature),1,iii)
  plot(fnow,trd,'r',fnow,exp(-abscoef),'g',fnow,error,'b');grid;drawnow
  axis([min(fnow) max(fnow) -0.2 1.2])
  fss=frequency_shift;
%eval(['save fileGL' num2str(pick_no(iii)) ' fnow trd abscoef tau2 beta frequency_shift'])
%eval(['save fileIND' num2str(pick_no(iii)) ' fnow trd abscoef tau2 beta frequency_shift'])
eval(['save file' num2str(pick_no(iii)) ' fnow trd abscoef tau2 beta fss kback sigsigR sigsigP fudge'])

  diff=[diff;error];
  end

nbeta=num2str(beta);nbsm=num2str(bsm);
ndur=num2str(duration);nfs=num2str(frequency_shift);
disp(['bsm=' nbsm '   beta=' nbeta '   tau2=' ndur '   fs=' nfs])
fprintf(1,'beta= %10.8f dofc = %10.8f \n',beta,duration);

final_diff=diff;

