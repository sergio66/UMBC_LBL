function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
         co2_param(band)
%function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
%         co2_param(band)
%this file sets the beta parameters etc using RAL data

%b2s,b2f=0 for all bands except the PR_deltpi bands
b2s=0.0;
b2f=0.0;

%use this d.of.c parameter everywhere!!!!!!!
duration_self = 7.475618919000517e-03;        %%before using cousin
duration_self = 7.240460331135358e-03;        %%after using cousin
duration_self = 6.785957156237349e-03;        %%redoing backgnd, using cousin
duration_self = 9.965661865269277e-03;

%duration=(pself*duration_pure+(ptotal-pself)*duration_for)/ptotal; 
ptotal= 961.70/1013.25; pself=26.80/1013.25;
duration     = 5.141075920725852e-03;         %%before using cousin
duration     = 5.607214396015617e-03;         %%after using cousin
duration     = 5.039772980195424e-03;         %%redoing backgnd, using cousin
duration     = 3.937305689748352e-03;

duration_for = (ptotal*duration - pself*duration_self)/(ptotal-pself);
duration     = 3.937305689748352e-03;

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
  beta_self = 8.551805493192096e-01;            %before using cousin
  beta_self = 8.758325975627885e-01;            %after using cousin
  beta_self = 8.719961873350909e-01;            %using cousin, redoing backgnd
  beta_self = 6.828529651649395e-01;

  %beta=(pself*beta_pure+(ptotal-pself)*beta_for)/ptotal; 
  beta     = 1.0;  
  beta_for = (ptotal*beta - pself*beta_self)/(ptotal-pself);
  beta_for = 1.009091390024152e+00;
  
else
  error('none of bands in this file matches YOUR band!!')
  end

