function [f,voi,fc]=...
     yrun(temperature,band,ptotal,pself,layeramt,prb,LVF,IO,birn,frequser)
%[f,y]=yrun(temperature,band,ptotal,pself,layeramt,frequser,prb,LVF,IO,birn)
%this function computes CO2 R-line mixing due to the 2350 band
%disp('bands: 2350')

%all this is put into structure "stuff"
%global bsm pressure_self pressure_for btz B0 temperature_ref;  
%global density Boltzmann mass_CO2 speed_light path_length pressure_ref; 
%global K_scale_mixing K_scale_voigt K_scale_lor; 
%global population population_t t_rawdata voi_back voi_pr; 
%global beta_delt_self beta_delt_air beta_pi_self beta_pi_air;  

if (band ~= 2350)
  error('need 2350 band!!!!!!!!!')
  end

if (prb == 'p')
  prb='P';
  end
if (prb == 'r')
  prb='R';
  end

cc=cputime;

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

f=frequser;

[jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr,stuff]=...
    loader(temperature,band,path_length,ptotal,pself,prb);
elower=elowerr;
[elower,jall]=efitter(jr,band,elowerr,elower,prb); 

%beta=(pself*stuff.beta_pure+(ptotal-pself)*stuff.beta_for)./ptotal;
[W_co2for,W_co2self]=wfunco2er(jr,elower,elowerr,w_selfr,w_forr,band,jall,...
                               temperature,stuff);

[trans_ampl,population_t]=...
      trans_pop(temperature,freqr,jr,elowerr,strenr,stuff);

W_plus=(pself*W_co2self+(ptotal-pself)*W_co2for)/stuff.pressure_ref;

% multiply off diagonals of W by beta 
% this is new, and conforms to the code in 
%/salsify/scratch4/Strow/Tobin_home/tobin/Co2q/B_sigpie/Run 
%and is done because of the last line in wfun1co2.m is 
%   dif=K+diag(wq_tape);     and not 
%   dif=beta*K1+diag(wq_tape); 
beta=stuff.beta;
W_plus=W_plus.*(beta+(1-beta)*eye(length(jr)));

rd=population_t.*trans_ampl;
len=length(rd);
sss=W_plus.*(ones(len,1)*rd').*(trans_ampl*ones(1,len));
sss=triu(sss,1);
denom=diag(W_plus).*rd.*trans_ampl;
denom=sum(denom);
sss=1+2*sum(sum(sss))/denom

%[lor,lormix,lor_birn,lormix_birn]=...
%         klormix_birn(f,freqr,W_plus,jr,temperature,beta,w_selfr,w_forr,...
%      trans_ampl,population_t,stuff,strenrt,ymix);
%fullmix4=full_mix4(freqr,f,W_plus,jr,w_selfr,w_forr,temperature,...
%      trans_ampl,population_t,stuff,birn);

%%%%this is a real full mixing calculation : does not matter what IO is
if ((LVF == 'F') | (LVF == 'f')) 
  voi=full2(freqr,f,W_plus,jr,w_selfr,w_forr,temperature,trans_ampl,...
                population_t,stuff,birn);  
elseif ((LVF == 'B') | (LVF == 'b'))  
  %do the fullmix first 
  voif=full2(freqr,f,W_plus,jr,w_selfr,w_forr,temperature,trans_ampl,... 
                population_t,stuff,birn);   
  %do the first mix second 
  ymix=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta); 
  voi1=voigtmix2(freqr,f,ymix,jr,...  
               temperature,w_forr,w_selfr,strenrt,stuff,layeramt,'V',birn);   
  voi=(ptotal-stuff.p1)/(stuff.p2-stuff.p1)*(voif-voi1) + voi1; 
else  
%%%%%%%%this is first order line mixing calculation using voigt/lor lineshape 
  if (IO == '1')
    ymix=y1s(jr,freqr,elowerr,strenr,W_plus,trans_ampl,stuff,beta);
  elseif (IO == '0')
%%%%%%%%this is zeroth order line mixing calculation using voigt/lor lineshape 
    ymix=zeros(size(freqr)); 
  else
    error('need IO to be 0 or 1 ');
    end
  voi=voigtmix2(freqr,f,ymix,jr,... 
               temperature,w_forr,w_selfr,strenrt,stuff,layeramt,LVF,birn);  
  end 
 

fc=freqr;

cc=cputime-cc
