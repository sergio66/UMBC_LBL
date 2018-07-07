function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
         co2_param(band,ptotal,pself)
%function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
%         co2_param(band,ptotal,pself)

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
%%this was used upto April 2000
%dur=[9.999997661522524 6.762294060993856 5.759869628913939 5.381980420978643];
%%%%%% this is new 4.00pm on April 13, 2000 - fudge factor == 0.9975 --> 0.96
%dur=[10.46423715226885 6.573035392098870 5.614703720544624 5.251501811591357];
%%%%%% this is new 4.00pm on April 14, 2000 - used newN2 and co2 inputs
%%%%%% as yesterday we had been using Temp=296K,unavg pressures instead of 
%%%%%% the correct avg temps and pressures; fudge factor == 0.9975 --> 0.96

%%%%%% this was fudged, using jmax = 60 on July 17, 2002
%%%%%% fudge = [1.00 0.99567 0.99959 1.00156]
dur=[10.46423654878965 6.544572534109112 5.612401733737971  5.259694144722923];
%%%%%%%%%%%%%%%%%%% BEST BEST BEST FOR RAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% this was used upto July 15, 2002
%%%%%% global pr2350_jmax = 600 in run7co2.m
%%%%%% dofudge = -1 in PR_sigsig/driver4um.m
%%%%%% dofudge = +1 in PR_sigsig/yrun_sigsig.m
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

  %%%%%this was before April 2000  - fudge factor = 1.0
%  besp=[6.858492724649634 8.895720470576258 ... 
%         9.320465583739526 9.453734030764889];
  %%%%%% this is new 3.00pm on April 13, 2000
%  besp=[6.831610402702522  8.942621789321642 ...
%        9.369734582258997 9.502437339211491];
  %%%%%% this is new 4.00pm on April 13, 2000 - fudge factor == 0.9975 --> 0.96
%  besp=[6.831266545391199 8.942573481127287 ...
%        9.369380351382158 9.502017375162639];
%%%%%% this is new 4.00pm on April 14, 2000 - used newN2 and co2 inputs
%%%%%% as yesterday we had been using Temp=296K instead of the correct avg 
%%%%%% temps and pressures; fudge factor == 0.9975 --> 0.96

%%%%%% this was fudged, using jmax = 60 on July 17, 2002
%%%%%% fudge = [1.00 0.993073 0.995393 0.987198]
  besp=[6.831266673446708 8.88062864436041 ...
        9.326215603644702 9.380372551552794];
%%%%%%%%%%%%%%%%%%% BEST BEST BEST FOR RAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% this was used upto July 15, 2002
%%%%%% global pr2350_jmax = 600 in run7co2.m
%%%%%% dofudge = -1 in PR_sigsig/driver4um.m
%%%%%% dofudge = +1 in PR_sigsig/yrun_sigsig.m
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

