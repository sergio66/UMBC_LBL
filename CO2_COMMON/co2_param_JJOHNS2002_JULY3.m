function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
         co2_param(band,ptotal,pself)
%function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
%         co2_param(band,ptotal,pself)

%this file sets the beta parameters etc using JJOHNS data in the 4 um region
%from the AIRS adta, as well as CAMEX1 and WINTEX and CLAMS, we see that the
%obs - calcs in the 4.3 um bandhead, says that our linemixing gives atmospheres
%that are too transperent ... so we need to "reduce" beta and thus increase abs
%inbetween the lines
%Looking at the JJOHNS data, this is eveident. So i did some more fits

%b2s,b2f=0 for all bands except the PR_deltpi bands
b2s=0.0;
b2f=0.0;

%%% first four from RAL data sets 1-4
%%% last  two are from JJOHNS data set 1 and avg JJOHNS data sets (2-4)
%%% note all were done at T = 296 K
pt=[2.6800e+01 1.6140e+02 5.6030e+02 9.6170e+02 1.0033e+03 1.0157e+03];  
ps=[2.6800e+01 2.6800e+01 2.6800e+01 2.6800e+01 1.0132e+00 5.0764e+01];  
pt=pt/1013.25; 
ps=ps/1013.25;
pf=pt-ps;
ratio=pf./pt;      %%%%%now we need to sort this!!!!!!
[Y,I] = sort(ratio);
ratio = ratio(I);

%use this d.of.c parameter everywhere!!!!!!!  
%%% first four from RAL data sets 1-4
%%% last  two are from JJOHNS data set 1 and avg JJOHNS data sets (2-4)
%dur=[10.46423654878965 6.573033770334661 5.614703762280506 5.251501801911941];
%dur=[dur 23.05402968282328 7.150923071809694];
dur=[10.46423654878965 6.573033770334661 5.614703762280506 5.251501801911941];
dur=[dur 23.05402968282328 5.650000000000000];
dur=dur(I)*1e-3;             %use sort index I
duration_self = dur(1);

if ((ptotal-pself)/ptotal < ratio(3))
  duration     = interp1(ratio,dur,(ptotal-pself)/ptotal);
else
  duration     = spline(ratio,dur,(ptotal-pself)/ptotal);
  end

if ((ptotal-pself)/ptotal > 0.00001)
  duration_for = (ptotal*duration - pself*duration_self)/(ptotal-pself);
else
  duration_for = duration_self;
  end


%%%%%%%%%%%%%%%%%%%% set up beta parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (band < 2200)    %%%%this sets up the P/Q/R branch SigPi, DeltPi mixing
  par_beta_self =   0.5358;            %self broadened     %%%%%used to be 0.5
  par_beta_for  =   0.6002;            %for  broadened     %%%%%used to be 0.62
  end 

if ((band == 618)|(band == 648)|(band == 662)|(band == 667)|(band == 720)|...
              (band == 791)|(band == 1932)|(band == 2080)|(band == 2129))
  %this is for the 15 um band
  %sigpi band   
  %for Q branch, duration is irrelevant, so this is used by PR_sigpi/loader.m
  %used for PQR branches
  beta_self = par_beta_self;    %self broadened  
  beta_for  = par_beta_for;     %foreign broadened
  
elseif  ((band == 668)|(band == 740)|(band ==2093))
  %this is for the 15 um band
  %deltpi band    for Q branch, duration is irrelevant

  %these parameters were gotten from the various Tobin files! no idea when!
  %however, in all the obs-calcs done prior to Sept 99, these were the 
  %parameters that were used. seem to be the best!!!!!!
  %version1
  beta_pi_self   = 0.5407505;  % fitted from pure Pi-Sigma data 
  beta_pi_air    = 0.5995593;  % fitted from air-broadened Pi-Sigma data 
  beta_delt_self = 0.4880067;  % fitted from pure Pi-Delta data 
  beta_delt_air  = 0.72172644; % fitted from air-broadened Pi-Delta data 

  %%%%%%%%since things are not working  Aug 2, 2000
  beta_delt_self=0.51; beta_delt_air=0.89; 
  beta_pi_self=0.54; beta_pi_air=0.60; % k89
  
  %beta_delt_self=0.9; beta_delt_air=0.9; beta_pi_self=0.9; 
  %beta_pi_air=0.9; % k89new
  
  %beta_delt_self=0.60; beta_delt_air=0.95; beta_pi_self=0.54; 
  %beta_pi_air=0.60; % k95
  
  %beta_delt_self=0.60; beta_delt_air=0.60; beta_pi_self=0.54; 
  %beta_pi_air=0.60; % k60
  
  %beta_delt_self=0.00; beta_delt_air=0.00; beta_pi_self=0.00; 
  %beta_pi_air=0.00; % k00

  beta_self     = beta_pi_self;
  beta_for      = beta_pi_air;
  b2s           = beta_delt_self;
  b2f           = beta_delt_air;

  %for PR branch, we DO NOT do complete line mixing ie we just use k/klor=0.5
  %and then use the d. of .c parameter there
  %so this file is called by PR_deltpi but just uses d. of. c
  %this is for the 4 um band
elseif  ((band == 2310)|(band == 2311)  |  ...                  %deltdelt
         (band == 2320)|(band == 2321)|(band == 2322) | ...     %pipi
         (band == 2350)|(band == 2351)|(band == 2352) | ...
         (band == 2353)|(band == 2354))                          %isgsig

  %%% first four from RAL data sets 1-4
  %%% last  two are from JJOHNS data set 1 and avg JJOHNS data sets (2-4)
  besp = [6.831266673446708 8.942573853443219 ...
          9.369380338865858 9.502017378026286];
  besp = [besp 6.428171629898014 9.217988039544647];
  besp = [6.831266673446708 8.942573853443219 ...
          9.369380338865858 9.502017378026286];
  besp = [besp 6.428171629898014 9.150000000000000];
  besp=besp(I)*0.1;               %use sort index I
  beta_self = besp(1);

  %beta=(pself*beta_pure+(ptotal-pself)*beta_for)/ptotal; 
  %if ((ptotal-pself)/ptotal < ratio(4))
  %  beta     = interp1(ratio,besp,(ptotal-pself)/ptotal);
  %else
  %  beta = besp(4);
  %  end
  if ((ptotal-pself)/ptotal < ratio(3))
    beta = interp1(ratio,besp,(ptotal-pself)/ptotal);
  else
    beta = spline(ratio,besp,(ptotal-pself)/ptotal);
    end

  if ((ptotal-pself)/ptotal > 0.00001)
    beta_for = (ptotal*beta - pself*beta_self)/(ptotal-pself);
  else
    beta_for = beta_self;
    end

else
  error('none of the bands in this file matches YOUR band!!')
  end

