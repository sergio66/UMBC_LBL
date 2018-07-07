function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
         co2_param_JJOHNS2002(band,ptotal,pself,region)
%function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
%         co2_param_JJOHNS2002(band,ptotal,pself,region)

%%% where you see %%% JJ      ==> fits the JJOHNS data very well
%%%               %%% JJ,RAL  ==> fits the JJOHNS + RAL quite very well

%%%%% to look at birnbaum factors %%%%%%
%path(path,'/home/sergio/SPECTRA/FORTRANLINUX');
%doc1=0.001;  y1 =birnbaumWORKS2(000:0.25:500,250,0.1,300,doc1);
%doc5=0.005;  y5 =birnbaumWORKS2(000:0.25:500,250,0.1,300,doc5);
%doc10=0.010; y10=birnbaumWORKS2(000:0.25:500,250,0.1,300,doc10);
%doc20=0.020; y20=birnbaumWORKS2(000:0.25:500,250,0.1,300,doc20);
%plot((000:0.25:500)-250,[y1; y5; y10; y20])

%%%%%%%%%%%%%%%%%%%% for 2350 band %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% use standard co2_param parameters everywhere EXCEPT
%%%   region == 1 ==> 2380-2395
%%%   region == 2 ==> 2390-2395
%%%   region == 3 ==> 2395-2405
%%%%%%%%%%%%%%%%%%%% for 2350 band %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% for 2351 band %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% use standard co2_param parameters everywhere EXCEPT
%%%   region == 1 ==> 2265-2290
%%%%%%%%%%%%%%%%%%%% for 2350 band %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this file sets the beta parameters etc using JJOHNS data in the 4 um region
%from the AIRS data, as well as CAMEX1 and WINTEX and CLAMS, we see that the
%obs - calcs in the 4.3 um bandhead, says that our linemixing gives atmospheres
%that are too transperent ... so we need to "reduce" beta and thus increase abs
%inbetween the lines
%Looking at the JJOHNS data, this is evident. So i did some more fits

%b2s,b2f=0 for all bands except the PR_deltpi bands
b2s=0.0;
b2f=0.0;

%%% first one from RAL data sets 1
%%% last  two are from JJOHNS data set 1 and avg JJOHNS data sets (2-4)
%%% note all were done at T = 296 K
pt=[2.6800e+01 1.0153e+03 1.0033e+03];  
ps=[2.6800e+01 5.0716e+01 1.0132e+00];  
pt=pt/1013.25; 
ps=ps/1013.25;
pf=pt-ps;
ratio=pf./pt;      %%%%%now we need to sort this!!!!!!
[Y,I] = sort(ratio);
ratio = ratio(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% basically, this is  2350 stuff ................................
%%%% but right at the end, we do some 2351 stuff at the 2280-2290 bandhole
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%use this d.of.c parameter everywhere!!!!!!!  
if region == 3
  %this is great for JJOHNS data, from 2395-2405 cm-1
  dur=[10.46423654878965 5.500000000000 22.97981366176344];   %%% JJ
  dur=[10.46423654878965 5.500000000000  5.00000000000000];   %%% JJ,RAL
elseif region == 2
  %this is great for JJOHNS data, from 2390-2395 cm-1
  dur=[10.46423654878965 5.5000000000000 22.97981366176344];   %%% JJ
  dur=[10.46423654878965 5.5000000000000  4.50000000000000];   %%% JJ,RAL
elseif region == 1
  %this is great for JJOHNS data, from 2380-2394 cm-1
  dur=[10.46423654878965 24.0000000000000 22.97981366176344];  %%% JJ,RAL
  end
dur=dur(I)*1e-3;             %use sort index I
duration_self = dur(1);

duration     = interp1(ratio,dur,(ptotal-pself)/ptotal,'linear','extrap');
if ((ptotal-pself)/ptotal > 0.00001)
  duration_for = (ptotal*duration - pself*duration_self)/(ptotal-pself);
else
  duration_for = duration_self;
  end

%%%%%%%%%%%%%%%%%%%% set up beta parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (band < 2200)    %%%%this sets up the P/Q/R branch SigPi, DeltPi mixing
  par_beta_self =   0.5358;            %self broadened  %%%%%used to be 0.5
  par_beta_for  =   0.6002;            %for  broadened  %%%%%used to be 0.62
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

  beta_delt_self=0.51; beta_delt_air=0.89; 
  beta_pi_self=0.54; beta_pi_air=0.60; % k89
  
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
         (band == 2353)|(band == 2354))                         %isgsig

  %%% first one from RAL data sets 1
  %%% last two from JJOHNS data set 1 and avg JJOHNS data sets (2-4) july3
  if region == 3
    %this is great for JJOHNS data, from 2395-2405 cm-1
    besp = [6.831266673446708 9.350000000000000 7.027982578711146]; %% JJ
    besp = [6.831266673446708 9.350000000000000 9.750000000000000]; %% JJ,RAL
  elseif region == 2
    %this is great for JJOHNS data, from 2390-2395 cm-1
    besp = [6.831266673446708 9.300000000000000 7.027982578711146]; %% JJ
    besp = [6.831266673446708 9.300000000000000 9.600000000000000]; %% JJ,RAL
elseif region == 1
    %this is great for JJOHNS data, from 2380-2394 cm-1
    besp = [6.831266673446708 7.000000000000000 7.027982578711146]; %% JJ,RAL
    end
  besp=besp(I)*0.1;               %use sort index I
  beta_self = besp(1);

  %beta=(pself*beta_pure+(ptotal-pself)*beta_for)/ptotal; 
  %if ((ptotal-pself)/ptotal < ratio(4))
  %  beta     = interp1(ratio,besp,(ptotal-pself)/ptotal);
  %else
  %  beta = besp(4);
  %  end
  beta = interp1(ratio,besp,(ptotal-pself)/ptotal,'linear','extrap');
  if ((ptotal-pself)/ptotal > 0.00001)
    beta_for = (ptotal*beta - pself*beta_self)/(ptotal-pself);
  else
    beta_for = beta_self;
    end
 %fprintf(1,'in co2_paramJJ, bs,bf = %8.6f %8.6f \n',beta_self,beta_for)

else
  error('none of the bands in this file matches YOUR band!!')
  end

if ((band == 2351) & (region == +1))
  %%%%%these comes from the CO2_COUSIN_FITS; see the ReadMe file Aug 05
  P =  [-1.9267e-01   5.3168e-01  -3.1755e-02   1.9469e-02];
  mod_beta = polyval(P,ptotal);
  mod_doc = 1.000;

  disp('modifying (region 1) beta and doc for 2351 band, 636 CO2 isotope')
  disp('using the 2351 COUSIN center')
  beta_for  = beta_for  * mod_beta;
  beta_self = beta_self * mod_beta;
  beta=(pself*beta_self+(ptotal-pself)*beta_for)/ptotal; 

  duration_for  = duration_for  * mod_doc;
  duration_self = duration_self * mod_doc;
  end

