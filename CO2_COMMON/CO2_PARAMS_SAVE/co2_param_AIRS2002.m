function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
         co2_param(band,ptotal,pself)
%function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
%         co2_param(band,ptotal,pself)

%%%see history of beta, doc in SPECTRA/CO2_COMMON/co2_param.m.JULY20_2002

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT   IMPORTANT %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% this file sets the general parameters for all bands. ie in year 1999-2002,
%%%% we fit the R branch of the RAL data, and got linemixing,doc parameters
%%%% which kept us very happy, especially in the 15 um region (the CAMEX and
%%%% WINTEX data certainly said he 2380-2405 region was improved over the
%%%% Genln2/Cousin lineshape)

%%%% then came AIRS data in mid 2002 which said
%%%%%  a) 15 um region is great!
%%%%%  b) there are still some problems in the 2380-2405 region (as suspected)
%%%%%  c) there is a big problem, high in the atmosphere, at 2280 cm-1

%%%% so we used JJOHNS data to fix the 2380-2405 region, by using and blending
%%%%   fits from 2380-2390, 2390-2400,2400-2405 regions
%%%% and used Cousin lineshape to fix the 2280 region ... get linemix and doc 
%%%%   parameters for the next strong isotope (2351) that mimic Cousin here
%%%% and used combination of 2351 (cousin) and 2350 (R) to fix 2200-2230 cm-1

%%% this file still sets general parameters for all bands, in 15 um and 4 um
%%% regions. However, it modifies the parameters for the 2351 band
%%% After that, PR_sigsig/yrun_sigsig.m further calls co2_param_JJOHNS2002.m 
%%% to fine tune the parameters in the 
%%%    2380-2390, 2390-2400,2400-2405 regions for the main (2350) isotope
%%%    2270-2290  regions for the nextr strongest (2351) isotope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT   IMPORTANT %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this file sets the beta parameters etc using RAL data in the 4 um region
% however OBS-CALCS imply that we have too transperent an atmosphere in
% the 2385-2400 region, so mebbe we should slightly INCREASE kabs
%%%%%%%%  thus we need to slightly REDUCE effect of linemixing beta %%%%%%%
%%%%%%%%  thus we need to slightly REDUCE effect of linemixing beta %%%%%%%
%%%%%%%%  thus we need to slightly REDUCE effect of linemixing beta %%%%%%%

%b2s,b2f=0 for all bands except the PR_deltpi bands
b2s=0.0;
b2f=0.0;

%%from RAL fits : note all were done at T = 296 K
pt=[2.6800e+01   1.6140e+02   5.6030e+02   9.6170e+02];  pt=pt/1013.25;
ps=[2.6800e+01   2.6800e+01   2.6800e+01   2.6800e+01];  ps=ps/1013.25;
pf=pt-ps;
ratio=pf./pt;

%use this d.of.c parameter everywhere!!!!!!!  from RAL fits
%%%%%%%%%%%%%%%%%%% BEST BEST BEST FOR RAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
dur=[10.46423654878965 6.573033770334661 5.614703762280506 5.251501801911941];
%%%%%%%%%%%%%%%%%%% BEST BEST BEST FOR RAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%AIRS FUDGE FACTOR
airs_fudge_dur = 1.000;  dur = dur * airs_fudge_dur;
dur=dur*1e-3;
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

%%%%%%%%%%%%%%%%%%% BEST BEST BEST FOR RAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
  besp=[6.831266673446708 8.942573853443219 ...
        9.369380338865858 9.502017378026286];
%%%%%%%%%%%%%%%%%%% BEST BEST BEST FOR RAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%AIRS FUDGE FACTOR
  airs_fudge_beta = 1.000;  besp = besp * airs_fudge_beta;
  besp=besp*0.1;
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
  %fprintf(1,'in co2_param, bs,bf = %8.6f %8.6f \n',beta_self,beta_for)

else
  error('none of bands in this file matches YOUR band!!')
  end

%%% to maintain compatibility with PR_sigsig/driver4um when doing datafits,
%%% have this parameter
dofudge = -1;  
if ((band == 2351) & (dofudge == -1))

  %%%%% these comes from the CO2_COUSIN_FITS; see the ReadMe file Aug 05
  %%%%% they sure seemed to improve the 2351 lineshapes, right at the
  %%%%% P-R bandhead at 2270-2290 cm-1

  P =  [-1.9267e-01   5.3168e-01  -3.1755e-02   1.9469e-02];
  mod_beta = polyval(P,ptotal);
  mod_doc = 1.000;

  %%% aug 5, 2002
  %%% at 2200 - 2250 cm-1, 
  %%% the RAL data seems to imply using the 2350 linemixing parameters 
  %%%   gives transmissions that are too large (t_obs - t_calcs < 0)
  %%% the RAL data seems to imply using the 2351 Cousin linemixing parameters
  %%%   gives transmissions that are too small (t_obs - t_calcs > 0)
  %%% so let beta --> (betaORIG + betaCOUSIN)/2
  %%%        beta --> (betaORIG + betaORIG * mod_beta)/2
  %%%        beta --> betaORIG (1 + mod_beta)/2

  %%% aug 9, 2002
  %%% this improves things (OBS-CALCS for RAL DATA)
  %%% but it still hurts the AIRS data!!!
  %%% so comment it out, and go back to orig RAL parameters from R2350
 % disp('modifying (basic) beta and doc for 2351 band, 636 CO2 isotope')
 % disp('using the AVERAGE of the 2350 global results, and 2351 COUSIN center')
 % beta_for  = beta_for  * (1 + mod_beta)/2;
 % beta_self = beta_self * mod_beta;
 % beta=(pself*beta_self+(ptotal-pself)*beta_for)/ptotal; 

  duration_for  = duration_for  * mod_doc;
  duration_self = duration_self * mod_doc;
  end

