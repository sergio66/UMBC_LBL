function [outwave,out_array,out_linemixarray]=run8co2_FULLlinemixUMBC(...
                                              gasID,fmin,fmax,profname,topts);

% same as run8co2_linemixUMBC.m except it ONLY DOES full linemixing
%

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% to create the line files per band eg CO2_MATFILES/H16/hit618.mat
%% make sure you have run driver_makeDAVEhitlin for latest HITRAN!!!
%% make sure you have run driver_makeDAVEhitlin for latest HITRAN!!!
%% make sure you have run driver_makeDAVEhitlin for latest HITRAN!!!
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%gasID = 2;
%fmin = 1905; fmax = 1930;
%profname = 'IPFILES/co2_std_layer70';

% this is essentially the run7 linemix code 
% 
% ********* also need lots of stuff from Global_Data_HITRAN2004 **************
% same as run7.m except it does "load mass.dat" and the rest of the isotope
% stuff, a little more smartly than does run7.m
% this meant I have to move the "load mass.dat command to AFTER hitread.m
%
% default database = '/asl/data/hitran/h12.by.gas';
% for all "efitter.m" code have changed   leastsq --> lsqnonlin
% for all "wfunco2er.m" code have changed leastsq --> lsqnonlin W/O JACOBIANS!
% ********* also need lots of stuff from Global_Data_HITRAN2004 **************

% same as run6co2 except there are a whole bunch of defaults
% also have introduced new parameter "band" and "PQRallowed" as well as
% parameter LVF now has the (G)ENLN2-92 lineshape,  which does    
%                             Q branch line mixing      w/o cousin
%                             no P/R branch line mixing w/ cousin
%      see history note (14) below
%
% [outwave,out_array]=run4co2(gasID,fmin,fmax,ffin,fmed,fcor,fstep,...
%       xnear,xmed,xfar,nbox,strength_far,strength_near,LVF,IO,birn,profname);
% becomes [outwave,out_array]=run7co2(gasID,fmin,profname,{topts})
% if you do 0.0005 res with 5 pt averaging to give 0.0025 
%   -- > [fr,l]=run6co2(2,705,730,0.0005,0.1,0.5,1,1,2.0,250,5,0.0,0.0,...
%              'F','1','b',BANDS,PQRallowed,profname); 
% becomes [fr,k]=run7co2(gasID,fmin,fmax,profname);
%
% so for example, to mimic GENLN2 in the 15um, do one of the following
% (a) turn OFF P/R line mixing :
%       topts.PQRallowed = [01 02 -14 +14 -15 +15 -11 +11 -12 +12 -13 +13];
%       topts.PQRallowed = [01 02                 -11 +11 -12 +12 -13 +13];
%     turn off birnbaum 
%       topts.birn = 'c'
%     turn off fist order linemixing (!!!!what a waste)
%       topts.IO = '0'
%     turn off full linemixing (!!!!what a waste)
%       topts.LVF = 'V'
% (b) turn on GENLN2-92 lineshape and cousin 
%       topts.LVF = 'G92'
%       topts.birn = 'c'
%       topts.HITRAN='/asl/data/hitran/h92.by.gas';
%      see history note (14) below

% another example, to turn off ALL linemixing and do run7.m cacluation,
%   toptsVoigt.LVF = 'V';   toptsVoigt.IO = '0'; toptsVoigt.birn = 'n';
%   toptsVoigt.band = [];  toptsVoigt.xfar = 25.0;

%this function computes the lines spectra for gasID
%
%  TYPE     VAR           DESCRIPTION              TYPICAL VALUE
%----------------------------------------------------------------------
%            these need to be passed in by the user
%------------------------------------------------------------------------
%integer   gasID          HITRAN gas ID                   2
%
%integer   fmin           minimum freq (cm-1)            705
%integer   fmax           maximum freq (cm-1)            730
%
%matrix   profname       this contains N x 5 matrix : layernum,p,pp,t,x
%                         where column 1 = layer number 
%                               column 2 = layer pressure (atm)
%                               column 3 = layer partial pressure (atm)
%                               column 4 = layer temperature (k)
%                               column 5 = layer gas amt (kmoles/cm2)
%
%------------------------------------------------------------------------
%  TYPE     VAR           DESCRIPTION                   DEFAULT VALUE
%            these are default values; can be reset using structure {topts}
%------------------------------------------------------------------------
%real      ffin           fine point spacing (cm-1)           0.0005
%real      fmed           medium point spacing (cm-1)         0.1
%real      fcor           coarse point spacing (cm-1)         0.5
%
%real      fstep          wide mesh width size (cm-1)         1.0
%real      xnear          near wing distance(cm-1)            1.0
%real      xmed           med wing distance(cm-1)             2.0
%real      xfar           far wing distance(cm-1)             250.0
%
%integer   nbox           boxcar sum size (odd integer)           5
%
%real      str_far        min line strength for far wing lines   0.0
%real      str_near       min line strength for near wing lines  0.0
%
%char      LVF            'L','l' for lorentz                     'F'
%                         'V','v' for voigt/vaan Huber
%                         'F','f' for full line mixin
%                         'G','g' for GENLN2-92 lineshape : 
%                            Q branch linemixing in the Q deltpi,Qsigpi bands
%                            all other bands/branches - only cousin lineshape
%      see history note (14) below
%
%string    HITRAN         path to HITRAN database   /asl/data/hitran/h12.by.gas
%
%char      IO             '1' for first order line mixing          '1'
%                         '0' for no line mixing
% note L1,L0; V1,V0 are different while F1,F0 are the same
%
%char     birn            'b','B' for include birnbaum              'b'
%                         'c','C' for include cousin
%                         'n','N' for no birnbaum/cousin 
%                         (always off for Q-sigpi,Q-deltpi)
%

%array    band           tells which bands to use for line mixing  -1 = all
% else choose from the following
%   CO2q_sigpi     = [618 648 662 667 720 791 1932 2080 2129];
%   CO2q_delpi     = [668 740 2093];
%   CO2pr_sigsig   = [2350 2351 2352 2353 2354];
%   CO2pr_deltdelt = [2310 2311];
%   CO2pr_pipi     = [2320 2321 2322];
%   CO2pr_sigpi    = [618 648 662 667 720 791];         
%   CO2pr_delpi    = [668 740];   %note k/klor for PR_delpi=0.5
% WARNING : bands that are NOT selected have their lines left "untouched" 
% ie they will be included in the "mainloop" as voigt/lorentz lines

%array   PQRallowed      tells which P,Q,R branches to use        -1 = all
%PQR is a code that tells us which bands to add in  01 = Q delt pi
%                                                   02 = Q sig  pi
%                                                  -14 = P sig  pi
%                                                  -15 = P delt pi
%                                                  +14 = R sig  pi
%                                                  +15 = R delt pi
%                                                  -11 = P sig  sig
%                                                  -12 = P delt delt
%                                                  -13 = P pi   pi
%                                                  +11 = R sig  sig
%                                                  +12 = R delt delt
%                                                  +13 = R pi   pi
%
%real      stren_mult     multiplies ALL strengths             1.0
%real      width_mult     multiplies ALL widths                1.0
%real      tsp_mult       multiplies ALL pressure shifts       1.0
% integer   iLayDo      if multiple layers read in (ie profname is a textfile with many layers)
%                          then if iLayDo == -1 do all lays, else select which ONE to do
%                          default                               -1
%
%integer   mainloop       execute main loop? (1 yes, -1 no)     1
%integer   linemixloop    execute linemix loop? (1 yes, -1 no)  1

%------------------------------------------------------------------------
%------------------------------------------------------------------------

%restrictions :
%(1) xnear <= xmed <= xfar        
%(2) xnear >= fstep
%(3) xmed/fmed  xnear/fmin   fstep/fmed   fstep/ffin       are integers
%(4)fstep/(nbox*ffin)        fcor/ffin                     are integers
%(5)(fmax-fmin)/fstep        (fmax-fmin)/fcor              are integers
%
%thus ffin,fmed,fcor are all 10/2^n   n <= 5

% HISTORY : 
% 23) v12p1 = +1 uses Matlab v12.1 leastsq fitting routines for wfunco2
%           = -1 uses Matlab v14+  fsolve or nonlinleastsq
%           =  0 uses Matlab v14+  lsqcurvefit
% 22) "which_isotope" is now a topts option to send in all isotopes (0) or
%      a list of isotopes to be used                           May 16, 2007
% 21) iFudgeType = -1 (default) for no fudge, so just use linemixing (with 
%                     the R2350 band being subdivided into 3 regions)
%                = +2 use my fudges on top of the 3 region linemixing at 4 um,
%                     derived from the Sept 2002 AIRS data
%                = +3 use Scott's fudges on top of the 3 region linemixing, 
%                     at 4 um, derived from the 2003-2004 AIRS data
% 20) added new parameter to do the following RemoveSomeBands (default = -1)
%       default : weak background lines, use voigt
%                 strong mixing line eg 720, 791, 2350 use linemixing
%                 so DO NOT remove linemixing bands
%       special : weak background lines, use voigt
%                 some strong mixing lines eg 720, 791, 2350 use linemixing
%                 other strong mixing lines eg 667,721, 2320 IGNORE
%                 this is used for building the NLTE background database in
%      /home/sergio/KCARTA/SRC/NONLTE2/sergio/NONLTE_MATLAB/loop_background.m
%      So in this case, RemoveSomeBands = +1
%      and if for example, bands has been set to 
%          topts.band = [2353 2354 2322];   %%so do these here, as "weak" bands
%      then the "removeCO2lines" is called twice, once to remove ALL linemix 
%      lines, once to remove only lines from [2353 2354 2322]; this way what
%      is done is ... all weak lines and only 2353,2354, 2322 lines; the rest 
%      eg 2350,2351,2352,2310,2320 are discarded

% 19) PR_sigsig/loader.m keeps f0 --> f0 instead of 
%        f0 --> f0 + P*dv(P) 
%        ie has NO pressure shift of linecenters
%     the shifts for the pipi, deltdelt, sigpi and deltpi lines are left at 0,
%     when looking at PR_XY/loader.m  
% 18) PR_sigsig yrun_sigsig and driver4um.m now have a "dofudge = -1/+1" 
%     parameter so that we can use directly run7co2 to do data fits
% 17) can now tell run7co2 to turn on/off the main loop or linemix loop
% 16) making PR_sigsig/yrun_sigsig more smart in that the (beta,d of c) indices
%     are now wavenumber dependent. This is so that AIRS, and RAL and JJOHNS 
%     data agree      
% 15) introduced parameters stren_mult,width_mult, tsp_mult which
%     multiply ALL strengths, foreign/self widths              Apr 03, 2002
% 14) parameter LVF now allows the (G)enln2-92 lineshape where there is
%                             Q branch line mixing      w/o cousin
%                             no P/R branch line mixing w/ cousin
% 13) added two new parameters : BANDS and PQR 
%    bands tell you where to do linemixing eg 720, 2093
%    PQR tells you to do P(-1), Q(0), R(+1) mixing
%  so if you are smart, you can do 720 Q branch mixing ONLY! *All* other lines
%    will be done with lorentz lineshape (ie 667, 791, 741 lines etc)
%    band = 720, PQR = 0. 
%  this means procedure removeCO2lines has been updated to removeCO2linesNEW
% **************************
%    WARNING : if you are trying to fit RAL or JOHN-JOHNS data, you need to 
%    carefully read CO2_RAL_FITS/Readme and mess around with 
%    removeco2linesNEW() vs removeco2lines(),  as well as symbolic links for 
%    removeCO2lines.m. also need to play with main loop on/off, linemix loop
%    on/off. As i said, carefully read CO2_RAL_FITS/Readme
% **************************
% 12) changed it to run7(gasID,v1,v2,profile,{topts})          Feb 25,2001
% 11) replaced UseToth   with HITRAN                           Feb 23,2001
% 10) can only do cousin if there is NO line mixing. Thus LVF must not be 'F'
%    and IO must not be '1' for this; instead LVF must be 'G'
% 9) the code internally can set LVF = 'B' to blend full/first order mixing
% 8) many of the loops are now MEXed
% 7) slowly putting in the PQR branches, as well as allowing user LVF option
% 6) this is JUST for CO2
% 5) allows the user to choose ALL isotopes of a gas (0) or one of the
%    isotopes
% 4) this differs from run3.m in that the code computes the optical depths
%    of all lines in EACH layer (k*length), and then does a UNION of all the
%    lines, so that lines do not turn on and off, driving the SVD crazy
% 3) we have parameter xnear xmed xfar to define fine,med,coarse meshes
% 2) better version of run2.m (vectorised as much as possible)
% 1) same as run1WORKS except that this computes the spectra for all 100
%    layers (run1WORKS is mainly a debugger for ONE layer)

% differs from run2a.m in that we have redefined what FINE,MEDIUM,COARSE
% lines to use. take [f3 f4] as the limits of the current wide mesh
% (a) fine lines are always within  [f3-xnear  f4+xnear]     xnear ~ 1 cm-1
%     ie these lines are always within +/- 1 cm-1 of current wide mesh
% (b) medium lines are always within 
%             [f3-xmed f3-xnear] U [f4+xnear f4+xmed]      xmed ~ 2 cm-1
%     ie these lines are always 1-2 cm-1 away from current wide mesh
% (c) far lines are always within 
%             [f3-xfar f3-xmed] U [f4+xmed f4+xfar]        xfar ~ 25 cm-1
%     ie these lines are always 2-25 cm-1 away from current wide mesh
% in addition, this means that the layer loop over coarse freq grid has been
% removed, so that there is only ONE loop over layers

%%% the bands are identified from v_l and v_u
%%% see /salsify/packages/Genln2/Genln2/new_linmix.f for details for the Qbands
%                  1   618 cm-1  3  2, 626 (Q50 - Q 2,  615 - 619 cm-1)
%                  2   648 cm-1  2  1, 636 (Q 2 - Q50,  648 - 652 cm-1)
%                  3   662 cm-1  2  1, 628 (Q 1 - Q50,  662 - 665 cm-1)
%                  4   667 cm-1  2  1, 626 (Q 2 - Q50,  667 - 671 cm-1)
%                  5   668 cm-1  4  2, 626 (Q 2 - Q81,  667 - 675 cm-1)
%                  6   720 cm-1  5  2, 626 (Q50 - Q 2,  718 - 721 cm-1)
%                  7   741 cm-1  8  4, 626 (Q81 - Q 2,  733 - 742 cm-1)
%                  8   791 cm-1  8  3, 626 (Q 2 - Q50,  791 - 794 cm-1)
%                  9  1932 cm-1  6  1, 626 (Q 2 - Q50, 1932 -1937 cm-1)
%                 10  2076 cm-1  8  1, 626 (Q 2 - Q50, 2076 -2080 cm-1)
%                 11  2093 cm-1 14  2, 626 (Q 2 - Q70, 2093 -2095 cm-1)
%                 12  2129 cm-1 15  2, 626 (Q50 - Q 2, 2128 -2130 cm-1)
%  MGQL - IS THE AFGL GLOBAL QUANTA ID OF THE LOWER STATE FOR EACH SET
%  MGQU - IS THE AFGL GLOBAL QUANTA ID OF THE UPPER STATE FOR EACH SET
%     DATA (MGQL(ISET),ISET=1,12)/2,1,1,1,2,2,4,3,1,1,2, 2/
%     DATA (MGQU(ISET),ISET=1,12)/3,2,2,2,4,5,8,8,6,8,14,15/
%%% under /asl/www/pub/strow/mixing/mixing_coef_tables.txt are the
%%% v1_v2_l_v3 quantum IDs
% sigpi
%    Band 10002  01101;  Isotope index 1;  freq = 618  cm-1; 
%    Band 01101  00001;  Isotope index 2;  freq = 648  cm-1;
%    Band xxxxx  yyyyy;  Isotope index 3;  freq = 662  cm-1;
%    Band 01101  00001;  Isotope index 1;  freq = 667  cm-1;  
%    Band 10001  01101;  Isotope index 1;  freq = 720  cm-1; 
%    Band 11101  10002;  Isotope index 1;  freq = 791  cm-1; 
%    Band 11102--00001;  Isotope index 1;  freq = 1932  cm-1; 
%    Band 11101--00001;  Isotope index 1;  freq = 2080  cm-1; 
%    Band 20001  01101;  Isotope index 1;  freq = 2129  cm-1;
%  deltpi
%    Band 02201  01101;  Isotope index 1;  freq = 668  cm-1; 
%    Band 11101  02201;  Isotope index 1;  freq = 740  cm-1; 
%    Band 12201  01101;  Isotope index 1;  freq = 2093  cm-1;  
% \hline 
%   $\Sigma \Sigma$ &   2350         &   Main            & 2230-2390 \\ 
%   $\Sigma \Sigma$ &   2280         &   C13             & 2180-2330 \\ 
%   $\Sigma \Sigma$ &   2320         &   C14             & 2220-2380 \\ 
%   $\Sigma \Sigma$ &   2325         &   Main            & 2230-2380 \\ 
%   $\Sigma \Sigma$ &   2330         &   Main            & 2230-2360 \\ 
% \hline 
%   $\Delta \Delta$ &   2320         &   Main            & 2230-2380 \\ 
%   $\Delta \Delta$ &   2260         &   C13             & 2180-2320 \\ 
%   $\Delta \Delta$ &   2310         &   C14             & 2230-2360 \\ 
% \hline 
%   $\Pi \Pi$ &   2310         &   Main            & 2230-2370 \\ 
%   $\Pi \Pi$ &   2250         &   C13             & 2180-2310 \\ 

%  Q_deltpi        =  [668 740 2093];
%  Q_deltpi_lower  =  [2   4   2   ];
%  Q_deltpi_upper  =  [4   8   14  ];
%  Q_isotope       =  [1   1   1   ];

%  Q_sigpi        =  [618 648 662 667 720 791 1932 2080 2129];
%  Q_sigpi_lower  =  [2   1   1   1   2   3   1    1    2 ];
%  Q_sigpi_upper  =  [3   2   2   2   5   8   6    8    15];
%  Q_isotope      =  [1   2   3   1   1   1   1    1    1 ];

%  P_sigsig        =  [2350 2351 2352 2353 2354];
%  P_sigsig_lower  =  [1    1    1     3    5];
%  P_sigsig_upper  =  [9    9    9     23   25];
%  P_isotope      =   [1    2    3     1    1];

%  P_deltdelt        =  [2310 2311];
%  P_deltdelt_lower  =  [4    4];
%  P_deltdelt_upper  =  [24   24];
%  P_isotope         =  [1    2];

%  P_pipi        =  [2320 2321 2322];
%  P_pipi_lower  =  [2    2    2];
%  P_pipi_upper  =  [16   16   16];
%  P_isotope      = [1    2    3];

% all the voigt, vhh fcns have been MEXed
%if brd > 1e-2 can use the Martin/Puerta voigt fcn (very fast and accurate)
%if brd < 1e-2 use the Tobin's voigt fcn (slow)

%        path(path,'/salsify/users/sergio/KCARTA/SPECTRA/Q_deltpi');
%        path(path,'Q_deltpi');

%%%%homepath='/salsify/users/sergio/KCARTA/SPECTRA/';

rand('seed', sum(1000*clock))

warning off
currentdir=what;            %get current directory
currentdir=currentdir.path;

homepath=[currentdir '/'];

str1=[homepath 'Q_sigpi'];      path(path,str1);
str1=[homepath 'Q_deltpi'];     path(path,str1);
str1=[homepath 'PR_sigpi'];     path(path,str1);
str1=[homepath 'PR_deltpi'];    path(path,str1);
str1=[homepath 'PR_sigsig'];    path(path,str1);
str1=[homepath 'PR_deltdelt'];  path(path,str1);
str1=[homepath 'PR_pipi'];      path(path,str1);

% wherever code says %%%%%%%%%%%%%% USER_WARN the user can reset some params

global quiet p2311_21_jmax p2350_jmax r2350_jmax pr2351_jmax stretch v12p1

stretch = 300;          %% this is magic number so that "hitread " looks for 
                        %% lines that are +/- stretch away from band center
%%%%%%%%   now find the 60 lowest j lines for the P_sigsig 2351 !!!!!!!!!
%%%%%%%%   or P_pipi 2321 or P_deltdelt 2311
%%%%%%%% these are the strongest isotope bands and give problems on the 
%%%%%%%% left side (high j's)
p2311_21_jmax = 600; %% this is magic number for PR_sigsig mixing, P branch
                     %% in case if we use too many j's, abs < 0 and so 
                     %% removeneg makes everything almost 0 (to left of the 
                     %% bandhead). but looking at RAL data indicates using 
                     %%all J's is fine!!!
r2350_jmax =  600;   %% this is a magic number for strong R_mixing
                     %% in case if we use too many j's, abs < 0
                     %% also leave dofudge = +1 in PR_sigsig/yrun_sigsig.m
                     %% also leave dofudge = +1 in PR_sigsig/driver4um.m
p2350_jmax =  80;    %% this is a magic number for strong P_mixing
                     %% in case if we use too many j's, abs < 0
                     %% also leave dofudge = +1 in PR_sigsig/yrun_sigsig.m
                     %% also leave dofudge = +1 in PR_sigsig/driver4um.m
pr2351_jmax = 600;   %% this is a magic number for next strongest isotope,
                     %% PR_sigsig mixing
                     %% in case if we use too many j's, abs < 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% optimal parameters as of July 20, 2002
%%% these seem to make AIRS data better in the 2380 and 2280 chunks
%%% but in 2255 chunk, they make knew/kold <= 1 (a little too small)
%%% this is probably because we stop P2350 linemixing at j == 80 ....
%%% while this makes the RAL lab data better at 2220, it seems to make 
%%% AIRS data bad at 2260
p2311_21_jmax = 600; 
r2350_jmax    = 600;
p2350_jmax    = 600;  %%%did this after FORCING 2351 stuff!!!! 
pr2351_jmax   = 600; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

quiet=-1;               %%avoid printing out debug statements

if (gasID ~= 2)
  error('this is a specialised code for CO2. use run8 for other gases');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set up defaults
%real      ffin           fine point spacing (cm-1)      0.0005
%real      fmed           medium point spacing (cm-1)    0.1
%real      fcor           coarse point spacing (cm-1)    0.5
%
%real      fstep          wide mesh width size (cm-1)      1.0
%real      xnear          near wing distance(cm-1)         1.0
%real      xmed           med wing distance(cm-1)          2.0
%real      xfar           far wing distance(cm-1)          250.0
%
%integer   nbox           boxcar sum size (odd integer)    5
%
%real      str_far        min line strength for far wing lines     0.0
%real      str_near       min line strength for near wing lines    0.0
%
%char      LVF            'L','l' for lorentz                     'F'
%                         'V','v' for voigt/vaan Huber
%                         'F','f' for full line mixing
%                         'G92','g92' for GENLN2-92 lineshape : 
%                            Q branch linemixing in the Q deltpi,Qsigpi bands
%                            all other bands/branches - only cousin lineshape
%
%string    HITRAN         path to HITRAN database   /asl/data/hitran/h12.by.gas
%char      IO             '1' for first order line mixing          '1'
%                         '0' for no line mixing
% note L1,L0; V1,V0 are different while F1,F0 are the same
%
%char     birn            'b','B' for include birnbaum              'b'
%                         'c','C' for include cousin
%                         'n','N' for no birnbaum/cousin 
%                         (always off for Q-sigpi,Q-deltpi)
%
%real      stren_mult     multiplies ALL strengths             1.0
%real      width_mult     multiplies ALL widths                1.0
%real      tsp_mult       multiplies ALL pressure shifts       1.0
% integer   iLayDo      which of the PROFNAME lays to do     -1
%
%integer   mainloop       execute main loop? (1 yes, -1 no)     1
%integer   linemixloop    execute linemix loop? (1 yes, -1 no)  1
%
% integer  v12p1          use Matlab v12.1 leastsq fitting routines  +1
%                         use current fsolve or nonlinlssq           -1

%%function [outwave,out_array]=...
%      run6co2(gasID,fmin,fmax,ffin,fmed,fcor,fstep,...
%      xnear,xmed,xfar,nbox,strength_far,strength_near,LVF,IO,birn,profname,
%      {topts});

v12p1         = +1;

ffin          = 0.0005;
fmed          = 0.1;
fcor          = 0.5;
fstep         = 1.0;
xnear         = 1.0;
xmed          = 2.0;
xfar          = 250.0;          %%%% <----------------- note this
if fmin > 2830 | fmin < 605
  newparams = runXtopts_params_smart(fmin);
  ffin  = newparams.ffin;
  fmed  = newparams.fmed;
  fcor  = newparams.fcor;
  fstep = newparams.fstep;
  xnear = newparams.xnear;
  xmed  = newparams.xmed;
  xfar  = newparams.xfar;
  xfar  = 250.0;          %%%% <----------------- note this
end

nbox          = 5;
strength_far  = 0.0;
strength_near = 0.0;
LVF           = 'F';
HITRAN        = '/asl/data/hitran/h12.by.gas';
HITRAN        = '/asl/data/hitran/h16.by.gas';
HITRAN        = '/asl/data/hitran/h20.by.gas';
IO            = '1';
birn          = 'b';
RemoveSomeBands = -1;            %%%% <--- assume we do not care about NLTE
band          = -1;              %%%% <---------------- do all bands
PQRallowed    = -1;              %%%% <---------------- do all branches
stren_mult    = 1.0;
width_mult    = 1.0;
tsp_mult      = 1.0;
mainloop      = 1;
linemixloop   = 1;
iLayDo        = -1;
iFudgeType    = -1;    %%% do not use /KCARTADATA/General/ChiFile/co2_4um_fudg*
which_isotope = 0;                             %default use all isotopes

%%%%%%%%%%%%%%% REMOVE THE BAND DATA FROM THE BACKGROUND %%%%%%%%%%%%%%%
%   CO2q_sigpi     = [618 648 662 667 720 791 1932 2080 2129];
%   CO2q_delpi     = [668 740 2093];
%   CO2pr_sigsig   = [2350 2351 2352 2353 2354];
%   CO2pr_deltdelt = [2310 2311];
%   CO2pr_pipi     = [2320 2321 2322];
%   CO2pr_sigpi    = [618 648 662 667 720 791];         
%   CO2pr_delpi    = [668 740];   %note k/klor for PR_delpi=0.5

CO2q_sigpi = [618 648 662 667 720 791 1932 2080 2129];
CO2q_delpi = [668 740 2093];
bandQ=union(CO2q_sigpi,CO2q_delpi);

%%%%%%%%%%%%%%% 2355 CO2pr_sigsig   = [2350 2351 2352 2353 2354 2355];
CO2pr_sigsig   = [2350 2351 2352 2353 2354];
CO2pr_deltdelt = [2310 2311];
CO2pr_pipi     = [2320 2321 2322];

CO2pr_sigpi    = [618 648 662 667 720 791];         
CO2pr_delpi    = [668 740];   %note k/klor for PR_delpi=0.5, so we reset to : 
CO2pr_delpi    = [];          %note k/klor for PR_delpi=0.5
%note we don't have to worry abput including PR mixing for sigpi,deltpi
%here, as function "removeCO2lines" takes care of things!!!!!!!!

bandPR=union(CO2pr_sigsig,CO2pr_deltdelt);
bandPR=union(bandPR,CO2pr_pipi);
bandPR=union(bandPR,CO2pr_sigpi);

band=union(bandQ,bandPR);
band0 = band;

%array   PQRallowed      tells which P,Q,R branches to use        -1 = all
%PQR is a code that tells us which bands to add in  01 = Q delt pi
%                                                   02 = Q sig  pi
%                                                  -14 = P sig  pi
%                                                  -15 = P delt pi
%                                                  +14 = R sig  pi
%                                                  +15 = R delt pi
%                                                  -11 = P sig  sig
%                                                  -12 = P delt delt
%                                                  -13 = P pi   pi
%                                                  +11 = R sig  sig
%                                                  +12 = R delt delt
%                                                  +13 = R pi   pi
PQRallowed = [01 02 -14 +14 -15 +15 -11 +11 -12 +12 -13 +13];
PQRallowed0 = PQRallowed;

allowedparams = [{'ffin'},     {'fmed'},        {'fcor'},     {'fstep'},  ...
                 {'xnear'},    {'xmed'},         {'xfar'},                ...
                 {'nbox'},     {'strength_far'}, {'strength_near'},       ...
                 {'LVF'},                        {'HITRAN'},              ...
                 {'IO'},       {'birn'},         {'band'},                ...
                 {'PQRallowed'},   {'RemoveSomeBands'}                    ...
                 {'stren_mult'},   {'width_mult'}, {'tsp_mult'},          ...
                 {'mainloop'},     {'linemixloop'}, {'iFudgeType'},       ...
                 {'which_isotope'},{'v12p1'},{'iLayDo'}];

%read in {topts}
if nargin == 5
  optvar = fieldnames(topts);
  for i = 1 : length(optvar)
   if (length(intersect(allowedparams,optvar{i})) == 1)
     eval(sprintf('%s = topts.%s;', optvar{i}, optvar{i}));
   else
     fprintf(1,'topts param not in allowed list ... %s \n',optvar{i});
     error('quitting run8co2');
   end
end
end

%v12p1 = -1; %%use v14   lsqnonlin
%v12p1 = +1; %%use v12.1 leastsq
if v12p1 > 0
  disp('doing v12.1 leastsq fits for wfunco2/efitter')
  addpath /home/sergio/SPECTRA/OPTIM12.1
  addpath /home/sergio/SPECTRA/OPTIM12.1/private
else
  disp('doing v14 lsqnonlin fits for wfunco2/efitter')
  rmpath /home/sergio/SPECTRA/OPTIM12.1
  rmpath /home/sergio/SPECTRA/OPTIM12.1/private
end

if band == -1
  band = band0;
end
if (PQRallowed == -1)
  PQRallowed = PQRallowed0;
end

if (abs(mainloop) ~= 1)
  error('need mainloop parameter to be -1/+1');
end
if (abs(linemixloop) ~= 1)
  error('need linemixloop parameter to be -1/+1');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% this uses MEX files for 
%%%%%%%% (1) boxint2      (boxint.f is the slow MATLAB file)
%%%%%%%% (2) vhh          (vhhF is fast file which is faster than tobin's)

%see that all parameters makes sense
z=checkpar(fstep,ffin,fmed,nbox,fcor,fmax,fmin,xnear,xmed,xfar);
if (z ~= 1) 
  error('there is an error in your parameters')
end
allowedLVF={'l','L','v','V','f','F','G92','g92'};
if (isempty(intersect(LVF,allowedLVF)))
  error('there is an error in your LVF parameter')
end
allowedbirn=['n','N','b','B','c','C'];
if (isempty(intersect(birn,allowedbirn)))
  error('there is an error in your birn parameter')
end
allowedIO=['1','0'];
if (isempty(intersect(IO,allowedIO)))
  error('there is an error in your IO (first/none) parameter')
end

if ((IO == '1') & ((birn=='C')|(birn=='c')))
  error('cannot have line mixing AND cousin turned on')
end

if ( ((LVF == 'f')|(LVF == 'F')) & ((birn=='C')|(birn=='c')) )
  error('cannot have line mixing AND cousin turned on; try LVF == ''g92'' ')
end

if ( ((LVF == 'f')|(LVF == 'F')) )
  IO='1';                    %this is a mixing calculation
end

%get all lines above strengthM value
if (strength_far > strength_near)
  strengthM=strength_near;
else
  strengthM=strength_far;
end

%%%%%%%%%%%%%% USER_WARN
% set these parameters as they determine what mininum Optical Depth 
% should be used  === k(line strengh)*max(VHH)
OptDepth_close = 1e-10;
OptDepth_close = 0.0;    %%%new, as of July 2007, to make ALL lines be used
OptDepth_far   = OptDepth_close*1000;

%%%%%%%%%%%%%% USER_WARN
% adjust parameter "mult" to adjust how far away from fmin,fmax we
% want to consider lines to include from HITRAN
mult=1;
low= fmin-mult*xfar;
high=fmax+mult*xfar;
low= fmin-stretch;
high=fmax+stretch;

%%%%%%%%%%%%%%%%%%% LOAD IN ISOTOPE MASSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
load mass.dat        %get the mass of the isotopes  
[x,y]=size(mass);
%first 32 are the numbers of the isotopes  <iGasID  NumIsotopes  0.0>
isotope_num=mass(1:32,2:2);
%mext bunch are the masses of the isotopes  <iGasID  MassIsotopes Abundance>
themass=mass(33:x,2:2);
%now pick up the gas masses we need
dummy1=sum(isotope_num(1:gasID-1));
dummy2=sum(isotope_num(1:gasID));
dummy=dummy1+1:dummy2;
mass_iso=themass(dummy);
liso=length(mass_iso);

%%%%%%%%%%%% IS THIS INTERACTIVE SESSION OR CLUNK THRU PROFILE %%%%%%%%%
useruser=-1;
if (useruser > 0)
  which_isotope=input('Enter which isotope : ');
end

%%%%%%%%%%%%%%%%%%% LOAD IN GAS PROFILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_gas_profile

if min(iLayDo) > 0 & max(iLayDo) > length(GasAmt)
  [iLayDo length(GasAmt)]
  error(' you have iLayDo > 0 but also greater than length(profile)!!!')
elseif min(iLayDo) > 0 & max(iLayDo) <= length(GasAmt)
  fprintf(1,'read in info for %4i layers, but only going to do layer %2i \n',length(GasAmt),iLayDo)
  GasAmt = GasAmt(iLayDo);
  press  = press(iLayDo);
  partpress = partpress(iLayDo);
  temperature = temperature(iLayDo);
end

pflag=quiet;

%which_isotope=0;                             %default use all isotopes
MinLayer=1; MaxLayer=length(GasAmt); Step=1; %default step thru ALL layers
NumLayers=(MaxLayer-MinLayer+Step)/Step;
TheLayers=MinLayer:Step:MaxLayer;

%%%%%%%%%%%%%%%%%% DEFINE OUTPUT ARRAYS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define output freq array, output 100 layer matrix
%output point spacing grid
output_spacing=nbox*ffin;
%eg 605:0.0025:629.9975;
fmaxUSE    = fmax-output_spacing;
outwave    = fmin:output_spacing:fmaxUSE; 
lennnn = (1:length(outwave))-1; outwave=fmin+lennnn*output_spacing;    %<-----
out_array  = zeros(NumLayers,length(outwave));

%%%%%%% READ IN RELEVANT DATA FROM HITRAN DATABASE  %%%%%%%%%%%%%%%%%%%%%%%%
%% fnamePRE='/salsify/scratch4/h96.by.gas/g';        %H96 -- old
%% fnamePRE='/salsify/scratch4/h98.by.gas/g';        %H98 -- KCARTA database 
%% fnamePRE='/salsify/scratch4/h2k.by.gas/g';        %H98 -- KCARTA database 
if (HITRAN(length(HITRAN)) == '/')
  fnamePRE = [HITRAN 'g' ];
else
  fnamePRE = [HITRAN '/g'];
end
fnamePOST='.dat';
fnameIN=int2str(gasID);
fname=[fnamePRE fnameIN fnamePOST]

%%%ZZZ
%%%%%%%%following line if we want to use latest Hitran* database
[lineORIG,hitran_version,hlist_qtips] = ...
  hitread(low,high,strengthM,gasID,fname,+1);
lineORIG.stren            = lineORIG.stren  * stren_mult;
lineORIG.tsp              = lineORIG.tsp    * tsp_mult;
lineORIG.abroad           = lineORIG.abroad * width_mult;
lineORIG.sbroad           = lineORIG.sbroad * width_mult;

%% subset for isotopes; not really needed in all its gory detail
subset_for_isotopes; 

%%%%%%%%

lineORIG = should_I_translate2oldHITparams(lineORIG,gasID,hitran_version);

exchangelinecenters = +1;  %%% run this ONCE with exchangelinecenters = +1 to 
                           %%% put in Hartman's line centers
exchangelinecenters = -1;  %%% run with old stuff "from HITRAN" 
if exchangelinecenters == 1
  disp('looking to replace HITRAN linecenters with Hartman linecenters')
  disp('by clearing out the CO2_MATFILES/hit*.mat files');
  lineIN = lineORIG;
  wn1    = fmin;
  wn2    = fmax;
  %% exchange HITRAN linecenters for Hartmann linecenters (in lineOUT); 
  [hartmann_bands,lineX]=co2lines_jmhartmann(wn1,wn2,lineIN,which_isotope,+1);
  lineORIG = lineX;
  rmer = ['!/bin/rm CO2_MATFILES/hit*.mat']; eval(rmer);
end

%%%%%%%%%%%%%%%%%%% LOAD IN ISOTOPE MASSES and do QTIPS %%%%%%%%%%%%%%%%%%%%%%
load_mass_iso_dat;
initializeABCDG_qtips;
if length(intersect(hitran_version,{'h92','h96','h98'})) == 1
  if gasID == 19
    lineORIG = reduce19(lineORIG);
  end
end
%%%%%%%%%%%%%%%%%%% LOAD IN ISOTOPE MASSES and do QTIPS %%%%%%%%%%%%%%%%%%%%%%

fc_orig=lineORIG.wnum;
fprintf(1,'\n read in hittomat \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PQR=[];
bandtype=[];
num_band=0;
band = sort(band);

PQR0=[];
bandtype0=[];
num_band0=0;

% useless dummy, to tell which version of the HITRAN database to use
if (gasID == 1)
  vers=12;
else
  vers=10;
end

fprintf(1,'\n finding start/stop wavenumbers of bands \n');
[hh,ll] = band_end_find(gasID,homepath,band,vers,strengthM,fname,hitran_version);
fprintf(1,'band               start                stop \n');
fprintf(1,'-------------------------------------------------\n');
pprriinntt = [band; ll; hh;];
fprintf(1,'%12.5f     %12.5f    %12.5f \n',pprriinntt);
fprintf(1,'-------------------------------------------------\n');

fprintf(1,'now going to "remove" BAND lines .. might take a few secs \n');
fprintf(1,'initial number of HITRAN lines = %3i \n',lineORIG.linct);

lineORIG0 = lineORIG;

if RemoveSomeBands > 0
  display(' ')
  band
  band0
  fprintf(1,'-->>>>>> start RemoveSomeBands == +1                  <<<<<--\n')
  fprintf(1,'--> WARNING : may have some 4um "Q" branches with no mixing \n')
  fprintf(1,'--> removing some of the 4um/15 um bands (eg for NLTE)!!! \n')
  fprintf(1,'--> finding start/stop wavenumbers of bands \n');
  [hh0,ll0] = band_end_find(gasID,homepath,band0,vers,strengthM,fname,hitran_version);
  fprintf(1,'band               start                stop \n');
  fprintf(1,'-------------------------------------------------\n');
  pprriinntt = [band0; ll0; hh0;];
  fprintf(1,'%12.5f     %12.5f    %12.5f \n',pprriinntt);
  fprintf(1,'-------------------------------------------------\n');

  line0 = lineORIG;
  lx0 = 0;
  for ii = 1:length(band0)
    if ((low <= hh0(ii)) & (high >= ll0(ii)))
      if quiet > 0
        fprintf(1,'ii band0(ii) = %3i  %4i \n',ii,band0(ii))
      end
      %latest,greatest
      [line0,num_band0,PQR0,bandtype0,lx] = ...
        removeCO2linesNEW(band0(ii),line0,low,high,...
                          num_band0,PQR0,bandtype0,...
                          vers,strengthM,homepath,...
                          fname,xfar,PQRallowed,exchangelinecenters);  
      lx0 = lx0 + lx;
      if quiet > 0
        fprintf(1,'\n');
      end
    end
  end
  display(' ')
  fprintf(1,'--> initial number of HITRAN lines = %3i \n',lineORIG.linct);
  fprintf(1,'--> removed %3i lines \n',lx0);
  fprintf(1,'--> final num of "background" lines = %3i \n \n',line0.linct);
  fprintf(1,'-->>>>>> end RemoveSomeBands == +1                  <<<<<--\n')
end

fprintf(1,'>>>--> WARNING : may have some 4um "Q" branches with no mixing \n')
fprintf(1,'>>>--> initial num of lines = %3i \n',lineORIG.linct);
lx0 = 0;
for ii=1:length(band)
  if ((low <= hh(ii)) & (high >= ll(ii)))
    if quiet > 0
      fprintf(1,'ii band(ii) = %3i  %4i \n',ii,band(ii))
    end
    [lineORIG,num_band,PQR,bandtype,lx] = ...
       removeCO2linesNEW(band(ii),lineORIG,low,high,num_band,PQR,bandtype,...
          vers,strengthM,homepath,fname,xfar,PQRallowed,exchangelinecenters);
       %latest,greatest
    lx0 = lx0 + lx;
    if quiet > 0
      fprintf(1,'\n');
    end
  end
end
fprintf(1,'>>>--> removed %6i lines \n \n',lx0);
fprintf(1,'>>>--> final num of "background" lines = %3i \n',lineORIG.linct);

%% for NLTE
%% [2310 2320 2350 2351 2352]
%%   save /carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/otherlines.mat lineORIG 
%% [2310 2311 2320 2321 2350 2351 2352]
%%   save /carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/otherlines2.mat lineORIG
%% [2310 2311 2320 2321 2350 2351 2352 2355]
%%   save /carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/otherlines3.mat lineORIG
%% [2310 2320 2350 2351 2353 2354]
%%   save /carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/otherlines4.mat lineORIG
%% [          2350               ]
%%   save /carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/otherlines5.mat lineORIG
%%error('oops')

if (pflag > 0)
  semilogy(lineORIG0.wnum,lineORIG0.stren,lineORIG.wnum,lineORIG.stren,'r')
  legend('before','after'); pause(0.1);
end

if RemoveSomeBands > 0
  lineORIG = line0;
  semilogy(lineORIG0.wnum,lineORIG0.stren,lineORIG.wnum,lineORIG.stren,'r')
  legend('before','after'); pause(0.1);
end

fprintf(1,'\n  final number of background lines = %3i \n',lineORIG.linct);
%%%%%%%%%%%%% GET READY TO LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n read in hittomat \n\n'); 
number_of_lines = lineORIG.linct; 
fprintf(1,'orig number of lines in HITRAN file = %8i \n',lineORIG.linct);  

nwide=(fmax-fmin)/fstep;    %number of wide meshes

%number of points in wide mesh (output resolution)
%eg if fstep=1.0, nbox*ffin=0.0025, then we really should have 
% fstep/(nbox*ffin)+1 = 401 pts
%but we are letting the "last" point be the "first" point of the next wide 
%mesh, so we only have 400 pts == fstep/(nbox*ffin) 
chunklen = find_chunklen(fstep,nbox,ffin); 
chunk=1:chunklen;

cpptotal=cputime; 

%%now do a UNION of all lines so that lines do not turn on and off, turning 
%%both SVD and Fast Models crazy. Also here include only the isotope the
%%user specified
%if ((lineORIG.linct > 0) & (NumLayers > 1))
%  line = findUnionNew(lineORIG,GasAmt,temperature,press,partpress,...
%                      A,B,C,D,G,hitran_version,mass_info,...
%                      gasID,mass_iso,OptDepth_close,which_isotope,TheLayers);
%  number_of_lines = line.linct;  
%else 
%  line = lineORIG; 
%  line = findUnionNew_oneLayer(lineORIG,GasAmt,temperature,press,partpress,...
%                      A,B,C,D,G,hitran_version,mass_info,...
%                      gasID,mass_iso,OptDepth_close,which_isotope,TheLayers);
%  number_of_lines = line.linct;  
%end 
%if (line.linct == 0)
%  disp('found NO lines in the HITRAN database!!!');
%end

%now do a UNION of all lines so that lines do not turn on and off, turning  
%both SVD and Fast Models crazy. Also here include only the isotope the 
%user specified 
if ((NumLayers > 1) & (number_of_lines > 0))  
  line = findUnionNew(hlist_qtips,...
                      lineORIG,GasAmt,temperature,press,partpress,...  
                      A,B,C,D,G,hitran_version,mass_info,...  
                      gasID,mass_iso,OptDepth_close,which_isotope,TheLayers);  
  number_of_lines = line.linct;    
elseif ((NumLayers == 1) & (number_of_lines > 0))  
  line=lineORIG;  
  line = findUnionNew_oneLayer(lineORIG,GasAmt,temperature,press,partpress,...
                      A,B,C,D,G,hitran_version,mass_info,...  
                      gasID,mass_iso,OptDepth_close,which_isotope,TheLayers);  
  number_of_lines = line.linct;    
end  
 
if (line.linct == 0)  
  disp('found NO lines in the HITRAN database!!!!'); 
else  
  fprintf(1,'found %8i lines in the HITRAN file ... proceeding\n',line.linct);
end 

if strengthM < 1e-50
  if (lineORIG.linct ~= line.linct)
    error('findUnionNew somehow missed some lines (strengthM == 0)')
  end
end 

if ((birn ~= 'c') & (birn ~= 'C')) 
  chichi = -1;
else
  chichi = +1;
end 

if ((LVF == 'L')|(birn == 'l')) 
  linetype = -1;
else
  linetype = +1;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% do background calc  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
domesh = +1;   %this does stuff the fast way!!!!!!! ie Genln2 way
domesh = -1;   %this does stuff the slow way!!!!!!! add ALL lines everywhere

domesh = +1;   %this does stuff the fast way!!!!!!! ie Genln2 way

%%%%%%%%%%%%%%%%% slow background code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%% this is the slow way to do the background 
%%%%%%%%%%%%%%%%% slow background code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

f1=fmin - ffin*(nbox-1)/2;    
%stop freq of current wide mesh  ... note we have  the -ffin*nbox because 
% we do not want to include the "last" point here 
f2=fmax - ffin*nbox + ffin*(nbox-1)/2;  
fine_freq=f1:ffin:f2; 
lennnn = (1:length(fine_freq))-1; fine_freq=f1+lennnn*ffin;        %<-----

[aaaa,bbbb] = size(mass_iso);
if aaaa > bbbb
  mass_iso = mass_iso';
end

if domesh < 0
  fprintf(1,'\n doing BACKGROUND using SLOOOOOOOOOOOOW method ... \n');
  for jj=MinLayer:Step:MaxLayer  

    if line.linct > 0 
      %INNER LOOP OVER LAYERS 1..100 = bottom -> top  
      nn    = (jj-MinLayer)/Step+1;
      fprintf(1,'layer number nn = %4i \n',nn);
      tempr = temperature(jj);  
      p     = press(jj);  
      ps    = partpress(jj);  
      qisotopes = getq_oldVSnew(hlist_qtips,hitran_version,A,B,C,D,G,...
                                mass_info,gasID,tempr);

      qfcn     = qoldVSqnew_fast(hlist_qtips,hitran_version,A,B,C,D,G,...
                                 qisotopes,mass_info,gasID,line,tempr);
      pwr      = line.abcoef; 
      for_brd  = line.abroad; 
      self_brd = line.sbroad; 
      freq     = line.wnum+press(jj)*line.tsp;%freq shift 
      energy   = line.els; 
      s0       = line.stren; 
      brd      = broad(p,ps,1.0,for_brd,self_brd,pwr,tempr,gasID); 
      strength = find_stren(qfcn,freq,tempr,energy,s0,GasAmt(jj)); 

      outvect  = loopco2(line.iso,mass_iso,brd,strength,freq,fine_freq,... 
                   tempr,line.linct,length(fine_freq),chichi,linetype,...
                   p-ps,ps); 

      scum            = boxint2(outvect,nbox); 
      out_array(nn,:) = out_array(nn,:)+scum; 
    end
  end
end  %if domesh < 0

if (pflag > 0)
  %fprintf(1,'\n hit return to continue ..... \n');
  plot(outwave,exp(-out_array(1,:))); pause(0.1);
end

%%%%%%%%%%%%%%%%% fast background code %%%% main loop %%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%% this is the fast way to do the background 
%%%%%%%%%% this is the way run6.m and GENLN2 do things
%%%%%%%%%%%%%%%%% fast background code %%%% main loop %%%%%%%%%%%%%%%%%%%%%%   

noloops = -1;
if mainloop > 0
  noloops = nwide;
end

% disp('for debigging in run8co2_linemixUMBC.m have set noloops == -1')
% disp('ret to continue');
% pause
% noloops = -1

if ((domesh > 0) & (noloops > 0))
 fprintf(1,'\n doing BACKGROUND using FAST method... \n');
%%%%%%these are the non mixing background weaker CO2 lines
%%%%%%%ZZZ
  for  ii=1:noloops    %OUTER LOOP OVER WIDE MESH 
 
    %%%%%%%% define the fine,medium, coarse resolution wavenumber grids %%%
    %this is defining the temporary fine freq resolution array 
    %start freq of current wide mesh 
    f1 = fmin+(ii-1)*fstep - ffin*(nbox-1)/2;    
    %stop freq of current wide mesh  ... we have  the -ffin*nbox because 
    % we do not want to include the "last" point here 
    f2 = fmin+ii*fstep - (ffin*nbox) + ffin*(nbox-1)/2;  
    fine_freq = f1:ffin:f2; 
    lennnn = (1:length(fine_freq))-1; fine_freq=f1+lennnn*ffin;        %<----- 

    %this is defining the temporary medium freq resolution array 
    f3 = fmin+(ii-1)*fstep;     %start freq of current wide mesh 
    f4 = fmin+(ii)*fstep;         %stop freq of current wide mesh 
    med_freq = f3:fmed:f4; 
    lennnn = (1:length(med_freq))-1; med_freq=f3+lennnn*fmed;        %<-----
 
    %this is defining the temporary coarse freq resolution array 
    f5 = fmin+(ii-1)*fstep;     %start freq of current wide mesh 
    f6 = fmin+(ii)*fstep;         %stop freq of current wide mesh 
    coarse_freq = f3:fcor:f4; 
    lennnn = (1:length(coarse_freq))-1; coarse_freq=f3+lennnn*fcor;    %<-----
 
    %%%%%%%%% define the current output resolution wavenumber grids %%%%%%% 
    %at the output wavenumber resolution, this is current wavenumber array 
    output_freq = f3:nbox*ffin:f4-nbox*ffin; 
    lennnn = (1:length(output_freq))-1; output_freq=f3+lennnn*nbox*ffin; %<----
    %occupying the following columns in the output array outwave,out_array 
    thechunk = chunk+(ii-1)*chunklen; 
 
    %%%%%%%% sort the lines into very near, medium near, far %%%%%%%%%%%%% 
    [very_near]   = sortbins(line,f3,f4,-1,xnear); 
    [medium_near] = sortbins(line,f3,f4,xnear,xmed); 
    [far_wing]    = sortbins(line,f3,f4,xmed,xfar); 
 
    fprintf(1,'\n mesh %3i start stop = %8.5f %8.5f',ii,f3,f4); 
 
    for jj=MinLayer:Step:MaxLayer 
      %INNER LOOP OVER LAYERS 1..100 =bottom->top 
      nn=(jj-MinLayer)/Step+1;

      tempr = temperature(jj); 
      p     = press(jj); 
      ps    = partpress(jj); 
      qisotopes = getq_oldVSnew(hlist_qtips,hitran_version,A,B,C,D,G,...
                                mass_info,gasID,tempr);

      %%%%%%%%%%%%%%% do the FINE GRID first (very near points) 
      %fprintf(1,'\n number of VERY NEAR lines = %3i',very_near.linct); 
 
      if very_near.linct > 0 
        qfcn     = qoldVSqnew_fast(hlist_qtips,hitran_version,A,B,C,D,G,...
                                   qisotopes,mass_info,gasID,very_near,tempr);
        pwr      = very_near.abcoef; 
        for_brd  = very_near.abroad; 
        self_brd = very_near.sbroad; 
        freq     = very_near.wnum+press(jj)*very_near.tsp;  %freq shift 
        energy   = very_near.els; 
        s0       = very_near.stren; 
        brd      = broad(p,ps,1.0,for_brd,self_brd,pwr,tempr,gasID); 
        strength = find_stren(qfcn,freq,tempr,energy,s0,GasAmt(jj)); 
      
        outvect = loopco2(very_near.iso,mass_iso,brd,strength,freq,...
                       fine_freq,tempr,very_near.linct,length(fine_freq),...
                       chichi,linetype,p-ps,ps); 
        scum = boxint2(outvect,nbox); 
        out_array(nn,thechunk) = out_array(nn,thechunk)+scum; 
 
      end  
 
      %%%%%%%%%%%%%%% do the MEDIUM GRID second  (medium near points) 
      %fprintf(1,'\n number of MEDIUM lines = %3i',medium_near.linct); 
 
      clear y 
      z = zeros(size(med_freq)); 
 
      if medium_near.linct > 0 
        qfcn     = qoldVSqnew_fast(hlist_qtips,hitran_version,A,B,C,D,G,...
                                  qisotopes,mass_info,gasID,medium_near,tempr);
        pwr      = medium_near.abcoef; 
        for_brd  = medium_near.abroad; 
        self_brd = medium_near.sbroad; 
        freq     = medium_near.wnum+press(jj)*medium_near.tsp;  %freq shift 
        energy   = medium_near.els; 
        s0       = medium_near.stren; 
        brd      = broad(p,ps,1.0,for_brd,self_brd,pwr,tempr,gasID); 
        strength = find_stren(qfcn,freq,tempr,energy,s0,GasAmt(jj)); 
 
        z = loopco2(medium_near.iso,mass_iso,brd,strength,freq,med_freq,... 
               tempr,medium_near.linct,length(med_freq),chichi,linetype,...
               p-ps,ps); 
        %do a spline interpolation of med_freq onto output_freq  
        xspline = spline(med_freq,z,output_freq);    %linear not good enough 
        out_array(nn,thechunk) = out_array(nn,thechunk)+xspline; 
 
      end 
 
      %%%%%%%%%%%%%%% do the COARSE GRID third  (farpoints) 
      %fprintf(1,'\n number of FAR lines = %3i',far_wing.linct); 
  
      clear y 
      z=zeros(size(coarse_freq)); 
 
      if far_wing.linct > 0 
        qfcn     = qoldVSqnew_fast(hlist_qtips,hitran_version,A,B,C,D,G,...
                                   qisotopes,mass_info,gasID,far_wing,tempr);
        pwr      = far_wing.abcoef; 
        for_brd  = far_wing.abroad; 
        self_brd = far_wing.sbroad; 
        freq     = far_wing.wnum+press(jj)*far_wing.tsp;  %freq shift 
        energy   = far_wing.els; 
        s0       = far_wing.stren; 
        brd      = broad(p,ps,1.0,for_brd,self_brd,pwr,tempr,gasID); 
        strength = find_stren(qfcn,freq,tempr,energy,s0,GasAmt(jj)); 
 
        z = loopco2(far_wing.iso,mass_iso,brd,strength,freq,coarse_freq,... 
               tempr,far_wing.linct,length(coarse_freq),chichi,linetype,...
               p-ps,ps); 
        %do a spline interpolation of coarse_freq onto output_freq  
        if (length(coarse_freq) == 2)   
          xspline = interp1(coarse_freq,z,output_freq);  %linear good enough  
        elseif (length(coarse_freq) == 3)   
          xspline = interp_quad(coarse_freq,z,output_freq);%quad good enough  
        elseif (length(coarse_freq) > 3)   
          xspline = spline(coarse_freq,z,output_freq);  %cubic good enough  
        end  
        out_array(nn,thechunk) = out_array(nn,thechunk)+xspline; 
      end 
    end           %inner loop over layers 
  end             %outer loop over wide meshes (frequencies) 
end               %if domesh > 0

out_backgnd = out_array;

%%%%%%%%%%%%%%%%%%%%%%%%%%% GUTS OF CODE == LINEMIXING %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% now we have to add on the PQR branches!!!!
%%%%%% do everything at FINE resolution for now
%%%%%%%%%%%%%%%%%%%%%%%%%%% GUTS OF CODE == LINEMIXING %%%%%%%%%%%%%%%%%%%%%%%

fc_bands=[];

LVF_user  = LVF;
birn_user = birn;
IO_user   = IO;

%linemixloop = -1;

noloops = -1;
if linemixloop > 0
  noloops = num_band;
end

fprintf(1,' \n');
fprintf(1,'running the linemixing portion with v12p1 = %3i \n',v12p1);
if ((num_band > 0) & (noloops > 0))

  fprintf(1,'\n doing PQR bands ... \n');
  f1 = fmin - ffin*(nbox-1)/2;    
  %stop freq of current mesh  ... note we have  the -ffin*nbox because 
  % we do not want to include the "last" point here 
  f2 = fmax - ffin*nbox + ffin*(nbox-1)/2;  
  fine_freq = f1:ffin:f2; 
  lennnn = (1:length(fine_freq))-1; fine_freq=f1+lennnn*ffin;        %<-----

  pqrfreq = fine_freq; 
 
  raRatio = zeros(1,NumLayers);
  rRatio  = 0.0;    %dummy value .. will change!!!!! kmix/klor=raRatio

%%%%%%%% this is the linemix loop; to turn off do for iii=1:-1
%%%%%%%%%ZZZ
  for iii = 1:noloops
    fprintf(1,' --> %6i : P(-ve)QR(+ve) type = %4i \n',bandtype(iii),PQR(iii));
 
    for jj = MinLayer:Step:MaxLayer  
      %INNER LOOP OVER LAYERS 1..100 = bottom -> top  
      nn = (jj-MinLayer)/Step+1;

      if quiet > 0
        fprintf(1,'layer number nn = %4i \n',nn);
      end

      tempr = temperature(jj);  
      p     = press(jj);  
      ps    = partpress(jj);  

%PQR is a code that tells us which bands to add in  01 = Q delt pi  
%                                                   02 = Q delt sig  
 
%                                                  -11 = P sig  sig  
%                                                  -12 = P delt delt  
%                                                  -13 = P pi   pi  
%                                                  -14 = P sig  pi  
 
%                                                  +11 = R sig  sig  
%                                                  +12 = R delt delt  
%                                                  +13 = R pi   pi  
%                                                  +14 = R sig  pi  
%bandtype tells you which band to do eg 740 etc  

      %worry about the transition from full to first order mixing
      %airslevels(50) = 1.6050e+02                ===== p2 (mb)
      %airslevels(50)/1013.25 = 1.5840e-01 (atm)
      %airslevels(75)/1013.25 = 2.0649e-02 (atm)  ===== p1

      p2 = 1.5840e-01;  %% atm = 160.45 mb
      p1 = 2.0649e-02;  %% atm = 20.91  mb

      %default set these to orig settings
      LVF  = LVF_user;
      birn = birn_user;
      IO   = IO_user; 
      if ((LVF_user == 'g92') | (LVF_user == 'G92'))
        if ((PQR(iii) == 01) | (PQR(iii) == 02)) 
          %%the Q deltpi,sigpi bands need first order line mixing, no cousin
          LVF  = 'V';
          birn = 'n';
          IO   = '1';
        elseif ((abs(PQR(iii))==11)|(abs(PQR(iii))==12)|(abs(PQR(iii))==13)) 
          %%the PR pipi,sigsig,deltdelt bands need NO line mixing, but cousin
          LVF  = 'V';
          birn = 'c';
          IO   = '0';
        elseif ((abs(PQR(iii))==14)|(abs(PQR(iii))==15))
          %%the PR deltpi,sigpi bands need NO line mixing, but cousin
          LVF  = 'V';
          birn = 'c';
          IO   = '0';
        end
      end

      if ((LVF == 'F') | (LVF == 'f'))
        disp('does only full line mixing')
        disp('does only full line mixing')
        disp('does only full line mixing')
        if (p > 0)   % do full mixing at any finite pressure
          if quiet > 0 
            fprintf(1,'doing only FULL mixing \n');
          end
          rRatio = raRatio(jj);
          [fhaha,yfull,fc,rRatio] = chooserCO2(tempr,bandtype(iii),p,ps,...
                GasAmt(jj),PQR(iii),'F',IO,birn,nbox,pqrfreq,outwave,...
                rRatio,homepath); 
          fc_bands = union(fc,fc_bands); 
          %plot(fc_bands);pause(0.1) 
          raRatio(jj) = rRatio;
          scum        = yfull;
        end
      else      %%%%%%%don't worry
        rRatio = raRatio(jj);
        [fhaha,y,fc,rRatio] = chooserCO2(tempr,bandtype(iii),p,ps,...
            GasAmt(jj),PQR(iii),LVF,IO,birn,nbox,pqrfreq,outwave,...
            rRatio,homepath); 
        fc_bands = union(fc,fc_bands); 
        %plot(fc_bands);pause(0.1) 
        scum        = y;
        raRatio(jj) = rRatio;
      end

      outold(nn,:)    = out_array(nn,:);
      out_array(nn,:) = out_array(nn,:) + scum; 

      pflag0 = pflag;
      %pflag = 1;
      if (pflag > 0)
    fprintf(1,'band = %6i P(-ve)QR(+ve) type = %4i \n',bandtype(iii),PQR(iii));
        subplot(311);
        semilogy(fhaha,scum,'b',outwave,out_array(nn,:),'r'); 
        title('current line mixing (blue) and total (red) so far');
        subplot(312);
        plot(fhaha,exp(-outold(nn,:))./exp(-out_array(nn,:)),'r'); 
        plot(fhaha,exp(-scum),'b',outwave,exp(-out_array(nn,:)),'r'); 
        subplot(313);
        semilogy(fhaha,scum)
        fprintf(1,'ret to continue ....');
        pause;
        disp(' ')
      end  
      pflag = pflag0;

    end 
  end 
end 

if iFudgeType > 0
  chifcn = makechi(outwave,iFudgeType);
  [iii,jjj] = size(out_array);
  out_array = out_array .* (ones(iii,1)*chifcn);
end

cpptotal = cputime-cpptotal

eval(['cd ' currentdir]);
warning on

fprintf(1,'ran the linemixing portion with v12p1 = %3i \n',v12p1);
out_linemixarray = out_array - out_backgnd;
