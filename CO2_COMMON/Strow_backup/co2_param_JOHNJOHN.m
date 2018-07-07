function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
         co2_param(band,ptotal,pself)
%function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
%         co2_param(band,ptotal,pself)
%this file sets the beta parameters etc

%b2s,b2f=0 for all bands except the PR_deltpi bands
b2s=0.0;
b2f=0.0;

%use this d.of.c parameter everywhere!!!!!!!
%this is from the files 8,9,10 in PR_sigsig/YIPPEE
duration_self = 9.244381131571670e-03;

%this is from the files 2,3,4 in PR_sigsig/YIPPEE
%duration=(pself*duration_pure+(ptotal-pself)*duration_for)/ptotal; 
%for files 2,3,4 
ptotal= 1.002026315789474; pself=0.05005263157894737;
%%duration     = 3.872568524276932e-03;   %****** from ALLFREQS ********
%%duration_for = 3.590130719741353e-03;   %****** from ALLFREQS ********
%%%%******use this!!!!!!!*********
duration     = 5.423029791802164e-03;   %****** from SOMEFREQ ********
duration_for = 5.222111747621331e-03;   %****** from SOMEFREQ ********
%%%%******use this!!!!!!!*********
%%duration     = 6.349696153099001e-03;   %****** from SOMEFREQ3 ********
%%duration_for = 6.197500138474120e-03;   %****** from SOMEFREQ3 ********


%%%%%%%ZZZ
duration     = 5.600e-3;
duration_for = (ptotal*duration - pself*duration_self)/(ptotal-pself);

if ((band == 618)|(band == 648)|(band == 662)|(band == 667)|(band == 720)|...
              (band == 791)|(band == 2080))
  %this is for the 15 um band
  %sigpi band   
  %for Q branch, duration is irrelevant, so this is used by PR_sigpi/loader.m
  %used for PQR branches
  beta_self = 0.5358;            %self broadened        %%%%%used to be 0.5
  beta_for  = 0.6002;            %foreign broadened     %%%%%used to be 0.62
  
elseif  ((band == 668)|(band == 740)|(band ==2093))
  %this is for the 15 um band
  %deltpi band    for Q branch, duration is irrelevant
  beta_pi_self   = 0.5407505;  % fitted from pure Pi-Sigma data 
  beta_pi_air    = 0.5995593;  % fitted from air-broadened Pi-Sigma data 
  beta_delt_self = 0.4880067;  % fitted from pure Pi-Delta data 
  beta_delt_air  = 0.72172644; % fitted from air-broadened Pi-Delta data 

  beta_self     = beta_pi_self;
  beta_for      = beta_pi_air;
  b2s           = beta_delt_self;
  b2f           = beta_delt_air;

  %for PR branch, we DO NOT do complete line mixing ie we just use k/klor=0.5
  %and then use the d. of .c parameter there
  %so this file is called by PR_deltpi but just uses d. of. c
  %this is for the 4 um band
elseif  ((band == 2310)|(band == 2311)  |  ...                      %deltdelt
         (band == 2320)|(band == 2321)|(band == 2322) | ...         %pipi
         (band== 2350)|(band== 2351)|(band == 2352)|(band ==2353)| ...
         (band==2354))                                              %isgsig
  %this is from the files 8,9,10 in PR_sigsig/YIPPEE
  beta_self = 7.178170098503445e-01; 
  %this is from the files 2,3,4 in PR_sigsig/YIPPEE
  %beta=(pself*beta_pure+(ptotal-pself)*beta_for)/ptotal; 
  %for files 2,3,4 ptotal= 1.002026315789474, pself=0.05005263157894737
%%  beta     = 9.815991792047214e-01;  %****** from ALLFREQS ********
%%  beta_for = 9.954682513847358e-01;  %****** from ALLFREQS ********
%%%%******use this!!!!!!!*********
  beta     = 9.368267525307101e-01;  %****** from SOMEFREQ ********       
  beta_for = 9.483417913856668e-01;  %****** from SOMEFREQ ********       
%%%%******use this!!!!!!!*********
  %%beta     = 9.160267964817331e-01;  %****** from SOMEFREQ3 ********       
  %%beta_for = 9.264482205086274e-01;  %****** from SOMEFREQ3 ********       
  %%%%%%%ZZZ
  %%beta_for=9.5e-01;

else
  error('none of bands in this file matches YOUR band!!')
  end
