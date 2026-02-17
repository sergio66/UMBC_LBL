function [outwave,out_array]=run8co2(gasID,fmin,fmax,profname,topts);

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% to create the line files per band eg CO2_MATFILES/H20/hit618.mat
%% make sure you have run driver_makeDAVEhitlin for latest HITRAN!!!
%% make sure you have run driver_makeDAVEhitlin for latest HITRAN!!!
%% make sure you have run driver_makeDAVEhitlin for latest HITRAN!!!
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% load /carrot/s1/sergio/IASIPRODUCTS_JACOBIANS/STD/g2_jac.mat
% semilogy(fout,sum(abs(jout'))); look for bands where the sum > 10^0 (==1)
% new default (topts.hartmann_linemix = +1)
%              is to do JMHARTMAN code from  555  to  855 cm-1
%                                            905  to 1105 cm-1
%                                      from 1830  to 2680 cm-1
%                      run8      code everywhere else
% if topts.hartmann_linemix = 0 then just call run8!!!!!!
% so most of the defaults are IF you call topts.hartmann_linemix = -1, which
%    is when the code goes to call run8co2_linemixUMBC.m
%
% ********* also need lots of stuff from Global_Data_HITRAN2004 **************
% same as run7.m except it does "load mass.dat" and the rest of the isotope
% stuff, a little more smartly than does run7.m
% this meant I have to move the "load mass.dat command to AFTER hitread.m
%
% default database = '/asl/data/hitran/h16.by.gas';
% default database = '/asl/data/hitran/h20.by.gas';
% for all "efitter.m" code have changed   leastsq --> lsqnonlin
% for all "wfunco2er.m" code have changed leastsq --> lsqnonlin W/O JACOBIANS!
% ********* also need lots of stuff from Global_Data_HITRAN2004 **************

% same as run6co2 except there are a whole bunch of defaults
% also have introduced new parameter "band" and "PQRallowed" as well as
% parameter LVF now has the (G)ENLN2-92 lineshape,  which does    
%                             Q branch line mixing      w/o cousin
%                             no P/R branch line mixing w/ cousin
%char      NIF              (V)oigt,F(I)rst Order,(F)ull linemix         'F'
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
%                         'G92','g92' for GENLN2-92 lineshape : 
%                            Q branch linemixing in the Q deltpi,Qsigpi bands
%                            all other bands/branches - only cousin lineshape
%      see history note (14) below
%char      NIF            'V','v' for nothing (voigt)
%                         'I','i' for first order linemix
%                         'F','f' for full line mixin
%
%string    HITRAN          path to HITRAN database   /asl/data/hitran/h16.by.gas
%                          path to GEISA  database   /asl/data/geisa/g15.by.gas
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
%
%integer   mainloop       execute main loop? (1 yes, -1 no)     1
%integer   linemixloop    execute linemix loop? (1 yes, -1 no)  1
%
%integer   v12p1          use Matlab v12.1 for fitting          +1
%                         use nonlinsq or fsolve for fitting    -1
% integer   iLayDo      if multiple layers read in (ie profname is a textfile with many layers)
%                          then if iLayDo == -1 do all lays, else select which ONE to do
%                          default                               -1
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
%  15) introduced "iLayDo" so user can choose which of the input lays he/she wants to process
% 24) Aug 2010 : iVersHartmann = 2010 or 2007 uses J-M Hartmann code supplied in 
%     year 20XX. Also allows user to include watervapor partial pressure 
%     (for the 2010 code version)
% 23) hartmann_linemix = +1 does hartman linemixing, which
%                           JMHARTMAN code from  555  to  855 cm-1 
%                                                905  to 1105 cm-1 
%                                          from 1830  to 2680 cm-1 
%                           run8                     code everywhere else 
%                      =  0 uses run8                code everywhere
%                      = -1 uses run8co2_linemixUMBC code everywhere
% 22A) can now use which_isotopes to subset for isotopes
%     if which_isotopes = 0            use all isotopes
%     if which_isotopes = [-1 iaList]  throw away isotopes in iaList, keep rest
%     if which_isotopes = [iaList]     keep isotopes in iaList, throw rest
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
%      is done is ... all wek lines and only 2353,2354, 2322 lines; the rest eg
%       2350,2351,2352,2310,2320 are discarded
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
%     char      NIF          (V)oigt,F(I)rstOrder,(F)ull linemix 'F'
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
%sigpi
% Band 10002  01101;  Isotope index 1;  freq = 618  cm-1; 
% Band 01101  00001;  Isotope index 2;  freq = 648  cm-1;
% Band xxxxx  yyyyy;  Isotope index 2;  freq = 662  cm-1;
% Band 01101  00001;  Isotope index 1;  freq = 667  cm-1;  
% Band 10001  01101;  Isotope index 1;  freq = 720  cm-1; 
% Band 11101  10002;  Isotope index 1;  freq = 791  cm-1; 
% Band 11102--00001;  Isotope index 1;  freq = 1932  cm-1; 
% Band 11101--00001;  Isotope index 1;  freq = 2080  cm-1; 
% Band 20001  01101;  Isotope index 1;  freq = 2129  cm-1;
% deltpi
% Band 02201  01101;  Isotope index 1;  freq = 668  cm-1; 
% Band 11101  02201;  Isotope index 1;  freq = 740  cm-1; 
% Band 12201  01101;  Isotope index 1;  freq = 2093  cm-1;  

%$\Sigma \Sigma$ &   2350         &   Main            & 2230-2390 \\ 
%$\Sigma \Sigma$ &   2280         &   C13             & 2180-2330 \\ 
%$\Sigma \Sigma$ &   2320         &   C14             & 2220-2380 \\ 
%$\Sigma \Sigma$ &   2325         &   Main            & 2230-2380 \\ 
%$\Sigma \Sigma$ &   2330         &   Main            & 2230-2360 \\ 
%\hline 
%$\Delta \Delta$ &   2320         &   Main            & 2230-2380 \\ 
%$\Delta \Delta$ &   2260         &   C13             & 2180-2320 \\ 
%$\Delta \Delta$ &   2310         &   C14             & 2230-2360 \\ 
%\hline 
%$\Pi \Pi$ &   2310         &   Main            & 2230-2370 \\ 
%$\Pi \Pi$ &   2250         &   C13             & 2180-2310 \\ 

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
%char      NIF            (V)oigt,F(I)rstOrder,(F)ull linemix 'F'
%
%string    HITRAN         path to HITRAN database   /asl/data/hitran/h16.by.gas
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
which_isotope = 0;                             %default use all isotopes
hartmann_linemix = -1;
iVersHartmann = 2010;                          % use his newest code
%%%%% most of this is for run8co2_linemixUMBC
LVF           = 'F';
NIF           = 'F';

do_HITRAN_vers                                        %% << set whether to use H96,H2k,H04,H08,H12,H16,H20,H24 >>
%% do_GEISA_vers; HITRANpathNyear = GEISApathNyear;   %% << if you want to use GEISA dbase, uncomment this to use G15 >>

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
iFudgeType    = -1;    %%% do not use /KCARTADATA/General/ChiFile/co2_4um_fudg*
iLayDo        = -1;

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
                 {'LVF'},      {'NIF'},          {'HITRAN'},              ...
                 {'IO'},       {'birn'},         {'band'},                ...
                 {'PQRallowed'},   {'RemoveSomeBands'}                    ...
                 {'stren_mult'},   {'width_mult'}, {'tsp_mult'},          ...
                 {'mainloop'},     {'linemixloop'}, {'iFudgeType'},       ...
                 {'which_isotope'},{'v12p1'},                             ...
                 {'hartmann_linemix'},{'iVersHartmann'},{'iLayDo'},       ...
                 {'WV_partialpressure'}];

%read in {topts}
if nargin == 5
  %% first check to see if earlier driver codes are using topts.HITRAN instead of topts.HITRANpathNyear
  topts = do_translate_params_HITRAN_2_HITRANpathNyear(topts);
  %% then proceed as usual      
  optvar = fieldnames(topts);
  for i = 1 : length(optvar)
    if (length(intersect(allowedparams,optvar{i})) == 1)
      eval(sprintf('%s = topts.%s;', optvar{i}, optvar{i}));
    else
      fprintf(1,'topts param not in run8co2 allowed list ... %s \n',optvar{i});
      error('quitting run8co2');
    end
  end
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

allowedNIF=['v','V','i','I','f','F'];
if (isempty(intersect(NIF,allowedNIF)))
  error('there is an error in your NIF parameter')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len_estimate = (fmax-fmin)/25*10000; %% in IR, 25 cm-1 = 10000 pts AFTER boxcar
if len_estimate * nbox > 200000
  fprintf(1,'estimate 25 cm-1 span = 10000 points after %2i boxcar integration \n',nbox)
  fprintf(1,'so for ffin = finemode, you need to x %5i the number of estimated output points \n',nbox)
  fprintf(1,'currently have parameter(MaxLen=200010) \n')
  error('please reduce fmax-fmin')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hartmann_linemix == -1
  disp('----------> branching to run8co2_linemixUMBC <----------------');
  if nargin == 4
    [outwave,out_array] = run8co2_linemixUMBC(gasID,fmin,fmax,profname);
  elseif nargin == 5
    toptsL = topts;
    if isfield(topts,'hartmann_linemix')
      toptsL = rmfield(toptsL,'hartmann_linemix');
    end
    [outwave,out_array] = run8co2_linemixUMBC(gasID,fmin,fmax,profname,toptsL);
  end
elseif hartmann_linemix == 0
  disp('----------> branching to run8 <----------------');
  if nargin == 4
    [outwave,out_array] = run8(gasID,fmin,fmax,profname);
  elseif nargin == 5
    toptsN = topts;
    if isfield(topts,'hartmann_linemix')
      toptsN = rmfield(toptsN,'hartmann_linemix');
    end
    if isfield(toptsN,'LVF')
      toptsN = rmfield(toptsN,'LVF');
      toptsN.LVG = topts.LVF;
    end
    if isfield(toptsN,'xfar') == 0
      disp('setting xfar == 250 cm-1')
      toptsN.xfar = 250;
    end
    [outwave,out_array] = run8(gasID,fmin,fmax,profname,toptsN);
  end
elseif hartmann_linemix == +1
  disp('----------> branching to run8co2_hartmann <----------------');
  if nargin == 4
    [outwave,out_array] = run8co2_hartmann(gasID,fmin,fmax,profname);
  elseif nargin == 5
    toptsH = topts;
    if isfield(toptsH,'hartmann_linemix')
      toptsH = rmfield(toptsH,'hartmann_linemix');
    end
    if isfield(toptsH,'LVF')
      toptsH = rmfield(toptsH,'LVF');
      toptsH.LVG = topts.LVF;
    end
    [outwave,out_array] = run8co2_hartmann(gasID,fmin,fmax,profname,toptsH);
  end
end

