function [outwave,out_array] = run8(gasID,fmin,fmax,profname,topts)

%tic
%tic0 = tic;

% ********* also need lots of stuff from Global_Data_HITRAN2004 **************
% same as run7.m except it does "load mass.dat" and the rest of the isotope
% stuff, a little more smartly than does run7.m
% this meant I have to move the "load mass.dat command to AFTER hitread.m
%
% default database = '/asl/data/hitran/h04.by.gas';
% default database = '/asl/data/hitran/h08.by.gas';
% default database = '/asl/data/hitran/h12.by.gas';
% default database = '/asl/data/hitran/h16.by.gas;'
% default database = '/asl/data/hitran/h20.by.gas;'
% ********* also need lots of stuff from Global_Data_HITRAN2004 **************

% same as run6.m, except parameters come in thru {param}; others are defaulted
% have also added on a new parameter : HITRAN version (1996,1998,2000)
% function [outwave,out_array]=run6(gasID,fmin,fmax,ffin,fmed,fcor,fstep,...
%       xnear,xmed,xfar,nbox,strength_far,strength_near,LVG,CKD,profname);
% function [outwave,out_array]=run7(gasID,fmin,fmax,profname,{topts})
% function [outwave,out_array]=run8(gasID,fmin,fmax,profname,{topts})
%
% if you straight away want 0.0025 res with no 5 pt averaging  
% also note the coarse freq point spacing is 0.5 ==> can do quad interp  
% also note the far wing coarse inclusion is 100 cm  
% [fr,k] = run6(3,705,755,0.0025,0.1,0.5,1,1,2,100,1,0.0,0.0,LVG,-1,profname);
%     or {topts.nbox = 1, topts.xfar = 100.0}
%     [fr,k] = run7(3,705,755,profname,topts);
% if you do 0.0005 res with 5 pt averaging to give 0.0025  
% also note the coarse freq point spacing is 1.0 ==> can do linear interp  
% [fr,k] = run6(3,705,730,0.0005,0.1,0.5,1,1,2,25,5,0.0,0.0,LVG,-1,profname);
%     [fr,k] = run7(3,705,755,profname);
%
% this function computes the lines spectra for gasID
% 
%   TYPE     VAR           DESCRIPTION              TYPICAL VALUE
% ----------------------------------------------------------------------
%             these need to be passed in by the user
% ------------------------------------------------------------------------
% integer   gasID          HITRAN gas ID                   3
% 
% integer   fmin           minimum freq (cm-1)            605
% integer   fmax           maximum freq (cm-1)            630
% 
% matrix   profname       this contains N x 5 matrix : layernum,p,pp,t,x
%                          where column 1 = layer number 
%                                column 2 = layer pressure (atm)
%                                column 3 = layer partial pressure (atm)
%                                column 4 = layer temperature (k)
%                                column 5 = layer gas amt (kilomoles/cm2)
%
% topts fields
% ------------------------------------------------------------------------
%   TYPE     VAR           DESCRIPTION                   DEFAULT VALUE
%             these are default values; can be reset using structure {topts}
% ------------------------------------------------------------------------
% real      ffin           fine point spacing (cm-1)             0.0005
% real      fmed           medium point spacing (cm-1)           0.1
% real      fcor           coarse point spacing (cm-1)           0.5
% 
% real      fstep          wide mesh width size (cm-1)           1.0
% real      xnear          near wing distance(cm-1)              1.0
% real      xmed           med wing distance(cm-1)               2.0
% real      xfar           far wing distance(cm-1)               25.0
% 
% integer   nbox           boxcar sum size (odd integer)            5
% 
% real      str_far        min line strength for far wing lines     0.0
% real      str_near       min line strength for near wing lines    0.0
% 
% char      LVG            (W)ithout basement == what LBLRTM has for "weak molecules" 
%                            ie other than CO2 or WV
%                          these are mostly Voivec = GENLN2 voigt
%                          (L)orentz,Voi(G)t,(V)anHuber             'V'
%                          (S)peed dependent voigt
%                          these are mostly Humlicek voigt                          
%                          Voi(GH)t,(VH)anHuber
%
% integer   CKD            continumm no (-1)                         -1
%                          yes water : (0,21,23) 
%                          yes O2,N2 : (+1) 
% string    HITRAN         path to HITRAN database   /asl/data/hitran/h16.by.gas
%                          path to GEISA  database   /asl/data/geisa/g15.by.gas
% 
% real      stren_mult     multiplies ALL strengths             1.0
% real      width_mult     multiplies ALL widths                1.0
% real      tsp_mult       multiplies ALL pressure shifts       1.0
%                          if negative, multiply main iso by abs(tsp_mult)
%                                        multiply other isos by 0
% integer   iLayDo      if multiple layers read in (ie profname is a textfile with many layers)
%                          then if iLayDo == -1 do all lays, else select which ONE to do
%                          default                               -1
%
% string   str_unc      uncertainty values for                   []
%                       wavenumber, intensity, air- and self-
%                       broadened halfwidths, temperature-dependence, and pressure shift
%                       can be + - R X which are max,min(-max),random,nothing
%                       so eg '++++++' would reset all the abpve line params to be X + dXmax
%                       so eg '------' would reset all the abpve line params to be X - dXmax
%                       so eg 'RRRRRR' would reset all the abpve line params to be X + random between
%                                  (dXmin and dXmax)
%                       so eg 'XXXXXX' would do nothing
%
% integer strongestline  picks out the strongest line so you can do your calc for one line only -1
%
% restrictions :
% (1) xnear <= xmed <= xfar        
% (2) xnear >= fstep
% (3) xmed/fmed  xnear/ffin   fstep/fmed   fstep/ffin       are integers
% (4)fstep/(nbox*ffin)        fcor/ffin                     are integers
% (5)(fmax-fmin)/fstep        (fmax-fmin)/fcor              are integers
% 
% thus ffin,fmed,fcor are all 10/2^n   n <= 5

%  HISTORY : 
%  16) introduced str_unc so can use the uncertainty codes to re-adjust parameters
%  15) introduced "iLayDo" so user can choose which of the input lays he/she wants to process
%  14) can now use which_isotopes to subset for isotopes
%      if which_isotopes = 0            use all isotopes
%      if which_isotopes = [-1 iaList]  throw away isotopes in iaList, keep rest
%      if which_isotopes = [iaList]     keep isotopes in iaList, throw rest
%  13) "which_isotope" is now a topts option to send in all isotopes (0) or
%       a list of isotopes to be used                           May 16, 2007
%  12) introduced parameters stren_mult,width_mult, tsp_mult which multiplies
%       ALL strengths, foreign/self widths, pressure shifts     Apr 03, 2002
%  11) changed it to run7(gasID,v1,v2,profile,{topts})          Feb 25, 2002
%  10) added HITRAN parameter                                   Feb 25, 2002
%  9) profile file name, which contains 5 column atmosphere profile
%  8) continuum for water, O2,N2
%  7) allows 3 different line shapes : lorentz, voigt or van huber
%  6) differsn from run4.m in that the loop over lines is now MEXED!!!
%  5) allows the user to choose ALL isotopes of a gas (0) or one of the
%     isotopes
%  4) this differs from run3.m in that the code computes the optical depths
%     of all lines in EACH layer (k*length), and then does a UNION of all the
%     lines, so that lines do not turn on and off, driving the SVD crazy
%  3) we have parameter xnear xmed xfar to define fine,med,coarse meshes
%  2) better version of run2.m (vectorised as much as possible)
%  1) same as run1WORKS except this should compute the spectra for all 100
%     layers (run1WORKS is mainly a debugger for ONE layer)

%  differs from run2a.m in that we have redefined what FINE,MEDIUM,COARSE
%  lines to use. take [f3 f4] as the limits of the current wide mesh
%  (a) fine lines are always within  [f3-xnear  f4+xnear]     xnear ~ 1 cm-1
%      ie these lines are always within +/- 1 cm-1 of current wide mesh
%  (b) medium lines are always within 
%              [f3-xmed f3-xnear] U [f4+xnear f4+xmed]      xmed ~ 2 cm-1
%      ie these lines are always 1-2 cm-1 away from current wide mesh
%  (c) far lines are always within 
%              [f3-xfar f3-xmed] U [f4+xmed f4+xfar]        xfar ~ 25 cm-1
%      ie these lines are always 2-25 cm-1 away from current wide mesh
%  in addition, this means that the layer loop over coarse freq grid has been
%  removed, so that there is only ONE loop over layers

% 
%  all the voigt, vhh fcns have been MEXed
% if brd > 1e-2 can use the Martin/Puerta voigt fcn (very fast and accurate)
% if brd < 1e-2 use the Tobin's voigt fcn (slow)

%  wherever code says %%%%%%%%%%%%%% USER_WARN the user can reset some params

%%%%%%%%%% this uses MEX files for 
%%%%%%%% (1) boxint2        (boxint.f is the slow MATLAB file)
%%%%%%%% (2) vhh,vhh1       (vhh is faster than tobin's vhh1)

% real      ffin           fine point spacing (cm-1)      0.0005
% real      fmed           medium point spacing (cm-1)    0.1
% real      fcor           coarse point spacing (cm-1)    0.5
% 
% real      fstep          wide mesh width size (cm-1)      1.0
% real      xnear          near wing distance(cm-1)         1.0
% real      xmed           med wing distance(cm-1)          2.0
% real      xfar           far wing distance(cm-1)          25.0
% 
% integer   nbox           boxcar sum size (odd integer)    5
% 
% real      str_far        min line strength for far wing lines     0.0
% real      str_near       min line strength for near wing lines    0.0
% 
% char      LVG            (W)ithout basement == what LBLRTM has for "weak molecules" 
%                            ie other than CO2 or WV
%                          these are mostly Voivec = GENLN2 voigt
%                          (L)orentz,Voi(G)t,(V)anHuber             'V'
%                          (S)peed dependent voigt
%                          these are mostly Humlicek voigt                          
%                          Voi(GH)t,(VH)anHuber
%
% integer   CKD            continumm no (-1)               -1
%                          yes water : (0,21,23) 
%                          yes O2,N2 : (+1) 
% string    HITRAN         path to HITRAN database   /asl/data/hitran/h08.by.gas
%                          path to GEISA  database   /asl/data/geisa/g15.by.gas
%
% real      stren_mult     multiplies ALL strengths             1.0
% real      width_mult     multiplies ALL widths                1.0
% real      tsp_mult       multiplies ALL pressure shifts       1.0
%
% integer   iLayDo         which of the PROFNAME lays to do     -1
% string    str_unc        uncertainty corrections              []
%
rand('seed', sum(1000*clock))

%  function [outwave,out_array] = run6(gasID,fmin,fmax,ffin,fmed,fcor,fstep,...
%   xnear,xmed,xfar,nbox,strength_far,strength_near,LVG,CKD,HITRAN,profname);

% %  default parameters for the 605-2830 cm-1 chunk
ffin          = 0.0005;
fmed          = 0.1;
fcor          = 0.5;
fstep         = 1.0;
xnear         = 1.0;
xmed          = 2.0;
xfar          = 25.0;
strongestline = -1;

%% xfar = max(xfar,25.0);    %%% <<<< new
if gasID == 2
  fprintf(1,' asking for run8 to do CO2, so set xfar = 500 cm-1 \n')
  fprintf(1,' asking for run8 to do CO2, so set xfar = 200 cm-1 \n')
  xfar = 500;
  xfar = 200;
end  

if (fmin > 2830 | fmin < 605)
  newparams = runXtopts_params_smart(fmin);
  ffin  = newparams.ffin;
  fmed  = newparams.fmed;
  fcor  = newparams.fcor;
  fstep = newparams.fstep;
  xnear = newparams.xnear;
  xmed  = newparams.xmed;
  xfar  = newparams.xfar;

  %% xfar = max(xfar,25.0);    %%% <<<< new
  if gasID == 2
    fprintf(1,' asking for run8 to do CO2, so set xfar = 500 cm-1 \n')
    fprintf(1,' asking for run8 to do CO2, so set xfar = 200 cm-1 \n')
    xfar = 500;
    xfar = 200;
  end  
end

nbox          = 5;
strength_far  = 0.0;
strength_near = 0.0;
LVG           = 'V';
if (gasID == 7 | gasID == 22)
  CKD  = +1;
elseif (gasID == 1)
  CKD = 4
else
  CKD = -1;
end
HITRAN        = '/asl/data/hitran/h16.by.gas';
HITRAN        = '/asl/rta/hitran/h16.by.gas';
HITRAN        = '/asl/rta/hitran/h20.by.gas';
stren_mult    = 1.0;
width_mult    = 1.0;
tsp_mult      = 1.0;
iLayDo        = -1;
which_isotope = 0;            % default use all isotopes
O2O3N2continuumVers = 3;      % used this is making H2008 database
O2O3N2continuumVers = 5;      % used this is making H2016 database together with MTCKD3.2
  %% allowed versions are [1,3,4,5] and [0 for only continuum (NO lines)]
  %% see CKDLINUX/calcong.F
str_unc       = [];

allowedparams = [{'ffin'},     {'fmed'},         {'fcor'},     {'fstep'}, ...
                 {'xnear'},    {'xmed'},         {'xfar'},                ...
                 {'nbox'},     {'strength_far'}, {'strength_near'},       ...
                 {'LVG'},      {'CKD'},                       {'HITRAN'},            ...
                 {'stren_mult'},   {'width_mult'},            {'tsp_mult'},          ...
                 {'which_isotope'}, {'O2O3N2continuumVers'},  {'iLayDo'},            ...
                 {'str_unc'},       {'strongestline'}];

%read in {topts} 
if nargin == 5
  optvar = fieldnames(topts);
  for i = 1 : length(optvar)
   if (length(intersect(allowedparams,optvar{i})) == 1)
     eval(sprintf('%s = topts.%s;', optvar{i}, optvar{i}));
   else
     fprintf(1,'topts param not in allowed list ... %s \n',optvar{i});
     error('quitting run8');
   end
 end
end

% Howard's latest suggestion
% option to override the above defaults with topts fields
% (other options set with kopts are passed along to rtcalc)
%if nargin == 5
%  optvar = fieldnames(topts);
%  for i = 1 : length(optvar)
%    vname = optvar{i};
%    if exist(vname, 'var')
%      eval(sprintf('%s = topts.%s;', vname, vname));
%    else
%      fprintf(1,'topts param %s not in allowed list ... \n', vname);
%      error('quitting run8');    
%  end
%end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (gasID == 1)         %water
  if (CKD >= 0)
    docontinuum = -1;
    fprintf(1,'\n this code WILL NOT DO water continuum ... \n');
    fprintf(1,'\n use run8water instead ... \n');
    error('exiting out of run8');
  end
end

if (gasID == 7 | gasID == 22)
  if (CKD ~= 1)
    disp('-------> warning : you did not turn on N2 or O2 continuum <----')
    disp('-------> warning : you did not turn on N2 or O2 continuum <----')
    disp('-------> warning : you did not turn on N2 or O2 continuum <----')
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see FORTRANLINUX/loop.F
%      = -3   voigt                  voigt1RI        HUMLICEK = used in HITTRAN LM 2018 
% rlvg = -2 : vhh without pedestal vhh2RI_basement VOIVEC = GENLN2
%        -1 : lerentz
%         0 : vanhuber               vhh2RI          VOIVEC = GENLN2    <<< DEFAULT >>>>
%         1 : voigt                  voigt2RI        VOIVEC = GENLN2
%         2 : speed dependent
%      = +3   vanhuber               vhh1RI           HUMLICEK = used in HITTRAN LM 2018

allowedline = {'w','W','l','L','v','V','g','G','s','S'};                                         %% these are basically GENLN2 voivec
allowedline = {'w','W','l','L','v','V','g','G','s','S','gh','GH','gH','Gh','vh','VH','vH','Vh'}; %% also allow Humlicek
okdokey = intersect(allowedline,LVG);
if (length(okdokey) ~= 1)
  error('lineshape must be W or L or V or G or S');
end

%% first five use GENLN2 Voivec
if ((LVG == 'w') | (LVG == 'W'))      %% without basement, just like WV
  lvgNum = -2;
elseif ((LVG == 'l') | (LVG == 'L'))
  lvgNum = -1;
elseif ((LVG == 'v') | (LVG == 'V'))
  lvgNum = 0;
elseif ((LVG == 'g') | (LVG == 'G'))
  lvgNum = 1;
elseif ((LVG == 's') | (LVG == 'S'))
  lvgNum = 2;
%% last two use Humlicek
elseif ((LVG == 'vh') | (LVG == 'VH') | (LVG == 'vH') | (LVG == 'Vh'))
  lvgNum = 3;
elseif ((LVG == 'gh') | (LVG == 'GH') | (LVG == 'gH') | (LVG == 'Gh'))
  lvgNum = -3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%see that all parameters makes sense
z = checkpar(fstep,ffin,fmed,nbox,fcor,fmax,fmin,xnear,xmed,xfar);
if (z ~= 1) 
  error('there is an error in your parameters')
end

%get all lines above strengthM value
if (strength_far > strength_near)
  strengthM = strength_near;
else
  strengthM = strength_far;
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
mult = 1;
low =  fmin-mult*xfar;
high = fmax+mult*xfar;

%%%%%%%%%%%%%%%%%%% IS THIS INTERACTIVE SESSION OR CLUNK THRU PROFILE %%%%%%%%%
useruser = -1;
if (useruser > 0)
  which_isotope = input('Enter which isotope : ');

  do_load   = 0;
  MinLayer  = 1; MaxLayer = 1; Step = 1;        %use only ONE "layer"
  NumLayers = 1;
  TheLayers = MinLayer:Step:MaxLayer;

  MGC = 8.314674269981136  ;  
  press       = input('Enter total pressure (in atm) : ');
  partpress   = input('Enter gas partial pressure (in atm) : ');
  temperature = input('Enter temperature (in K) : ');
  GasAmt      = input('Enter path cell length (in cm) ');
  %change to kmoles cm-2 
  GasAmt = GasAmt*101325*partpress/1e9/MGC/temperature; %change to kmoles/cm2 
end

%%%%%%%%%%%%%%%%%%% LOAD IN GAS PROFILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_gas_profile   %%% load profname

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

MinLayer  = 1; MaxLayer = length(GasAmt); 
Step      = 1; %default step thru ALL layers
NumLayers = (MaxLayer-MinLayer+Step)/Step;
TheLayers = MinLayer:Step:MaxLayer;

%%%%%%%%%%%%%%%%%% DEFINE OUTPUT ARRAYS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define output freq array, output 100 layer matrix
%output point spacing grid
output_spacing = nbox*ffin;
%eg 605:0.0025:629.9975;
fmaxUSE = fmax-output_spacing;
outwave = fmin:output_spacing:fmaxUSE;  
lennnn  = (1:length(outwave))-1; 
outwave = fmin+lennnn*output_spacing;   %<-----

out_array = zeros(NumLayers,length(outwave));

%%%%%%% READ IN RELEVANT DATA FROM HITRAN DATABASE  %%%%%%%%%%%%%%%%%%%%%%%%
%% fnamePRE = '/salsify/scratch4/h96.by.gas/g';        %H96 -- old
%% fnamePRE = '/salsify/scratch4/h98.by.gas/g';        %H98 -- KCARTA database 
%% fnamePRE = '/salsify/scratch4/h2k.by.gas/g';        %H98 -- KCARTA database 
if (HITRAN(length(HITRAN)) == '/')
  fnamePRE = [HITRAN 'g' ];
else
  fnamePRE = [HITRAN '/g'];
end
fnamePOST = '.dat';
fnameIN   = int2str(gasID);
fname     = [fnamePRE fnameIN fnamePOST];

%%%%%%%%%%ZZZ
%%%%%%%%following line if we want to use latest Hitran* database
[lineORIG,hitran_version,hlist_qtips] = ...
   hitread(low,high,strengthM,gasID,fname,-1);

addpath /home/sergio/MATLABCODE
lineX = lineORIG;
lineORIG = adjust_for_uncertainty(lineORIG,str_unc,gasID,fmin);

lineORIG.stren            = lineORIG.stren  * stren_mult;

if tsp_mult >= 0
  lineORIG.tsp              = lineORIG.tsp    * tsp_mult;
elseif tsp_mult < 0
  %% do what LBLRTM does namely only do tsp mult for MAIN isotope, and no line shift for isotopes
  junk = find(lineORIG.iso == 1);
    lineORIG.tsp(junk)  = lineORIG.tsp(junk) * abs(tsp_mult);
  junk = find(lineORIG.iso > 1);
    lineORIG.tsp(junk)  = lineORIG.tsp(junk) * abs(tsp_mult) * 0;
end

lineORIG.abroad           = lineORIG.abroad * width_mult;
lineORIG.sbroad           = lineORIG.sbroad * width_mult;
%%%%%%%%

%%%%%%%%%%%%%%%%%%% LOAD IN ISOTOPE MASSES and do QTIPS %%%%%%%%%%%%%%%%%%%%%%
%% subset for isotopes; not really needed in all its gory detail
%subset_for_isotopes
[lineORIG,which_isotope] = subset_for_isotopes(lineORIG,which_isotope);

lineORIG = should_I_translate2oldHITparams(lineORIG,gasID,hitran_version);

load_mass_iso_dat;
initializeABCDG_qtips;
if length(intersect(hitran_version,{'h92','h96','h98'})) == 1
  if gasID == 19
    lineORIG = reduce19(lineORIG);
  end
end

if gasID == 27 & ~strfind(hitran_version,'h12')
  %% had problems with number of isotopes (2) and qtips08.m
  if nargin == 5
    [outwave,out_array] = run8_g27(gasID,fmin,fmax,profname,topts);
  else
    [outwave,out_array] = run8_g27(gasID,fmin,fmax,profname);
  end
  return
end

%%%%%%%%%%%%%%%%%%% LOAD IN ISOTOPE MASSES and do QTIPS %%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'\n read in hittomat \n\n');
number_of_lines = lineORIG.linct;
fprintf(1,'orig number of lines in HITRAN file = %8i \n',lineORIG.linct);  

nwide = (fmax-fmin)/fstep;    %number of wide meshes

%number of points in wide mesh (output resolution)
%eg if fstep = 1.0, nbox*ffin = 0.0025, then we really should have 
% fstep/(nbox*ffin)+1  =  401 pts
%but we are letting the "last" point be the "first" point of the next wide 
%mesh, so we only have 400 pts == fstep/(nbox*ffin) 
chunklen = find_chunklen(fstep,nbox,ffin); 
chunk=1:chunklen;

cpptotal = cputime; 

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
  line = lineORIG;
  line = findUnionNew_oneLayer(lineORIG,GasAmt,temperature,press,partpress,...
                      A,B,C,D,G,hitran_version,mass_info,...
                      gasID,mass_iso,OptDepth_close,which_isotope,TheLayers);
  number_of_lines = line.linct;  
end

if strongestline == +1
  boo = find(line.stren == max(line.stren),1);
  semilogy(line.wnum,line.stren,line.wnum(boo),line.stren(boo),'ro');
  thefields = fieldnames(line);
  for ii = 1 : length(thefields)
    str = ['junk = line.' thefields{ii} ';'];
    eval(str);
    [mmjunk,nnjunk] = size(junk);
    if mmjunk == 1 & nnjunk == 1
      str = ['newjunk.' thefields{ii} ' = line.' thefields{ii} ';'];
      eval(str);
    elseif mmjunk == 1 & nnjunk > 1
      str = ['newjunk.' thefields{ii} ' = line.' thefields{ii} '(boo);'];
      eval(str);
    elseif mmjunk > 1 & nnjunk > 1
      str = ['newjunk.' thefields{ii} ' = line.' thefields{ii} '(boo,:);'];
      eval(str);
    end
  end
  newjunk.linct = 1;
  line = newjunk;
  number_of_lines = line.linct;  
  fprintf(1,'rounding linecenter from %12.5f to %12.5f \n',line.wnum,round(line.wnum))
  line.wnum = round(line.wnum);
  disp('found strongest line')
  line
end

if (number_of_lines <= 0)
  fprintf(1,'found NO lines in the HITRAN file \n');
else
  fprintf(1,'found %8i lines in the HITRAN file ... proceeding\n',line.linct);
end

%% had to put in a check for number_of_lines > 0 since for 0 < f < 320, 
%% there are NO lines for nitrogen (gid 22), but there is continuum!!!
%% so we need code to skip thru the lines calc, but do the continuum calc,
%% without barfing

if strengthM < 1e-50 & number_of_lines > 0 & which_isotope == 0
  if (lineORIG.linct ~= line.linct)
    disp('WARNING .... findUnionNew somehow missed some lines (strengthM == 0)')
  end
end 

[aaaa,bbbb] = size(mass_iso);
if aaaa > bbbb
  mass_iso = mass_iso';
end

iTestIso = -1;
if iTestIso > 0
  clf
  donk = find(line.wnum >726.23 & line.wnum < 726.24);
  donk = find(line.wnum >726.354 & line.wnum < 726.360);
  donk = find(line.iso >= 2);
  semilogy(line.wnum,line.stren,...
         line.wnum(donk),line.stren(donk),'ro'); grid
  %axis([726.20 726.45 1e-30 1e-22]); grid on
    index = donk;
    line.linct  = length(index);
    line.igas   = line.igas; 
    line.iso    = line.iso(index);  
    line.wnum   = line.wnum(index); 
    line.stren  = line.stren(index); 
    line.tprob  = line.tprob(index); 
    line.abroad = line.abroad(index); 
    line.sbroad = line.sbroad(index); 
    line.els    = line.els(index); 
    line.abcoef = line.abcoef(index); 
    line.tsp    = line.tsp(index); 
    line.iusgq  = line.iusgq(:,index); 
    line.ilsgq  = line.ilsgq(:,index);
    line.gasid  = line.gasid(index);
end

%do near wing calc
%%%%%%%%%%%%%%%%%% ZZZZZZ
if (number_of_lines > 0)
  for  ii = 1:nwide    %OUTER LOOP OVER WIDE MESH
  %%for  ii = 1:-1    %OUTER LOOP OVER WIDE MESH

    %%%%%%%%% define the fine,medium, coarse resolution wavenumber grids %%%%
    %this is defining the temporary fine freq resolution array
    %start freq of current wide mesh
    f1 = fmin+(ii-1)*fstep - ffin*(nbox-1)/2;   
    %stop freq of current wide mesh  ... note we have  the -ffin*nbox because
    % we do not want to include the "last" point here
    f2 = fmin+ii*fstep - (ffin*nbox) + ffin*(nbox-1)/2; 
    fine_freq = f1:ffin:f2;
    lennnn = (1:length(fine_freq))-1; fine_freq = f1+lennnn*ffin;        %<-----

    %this is defining the temporary medium freq resolution array
    f3 = fmin+(ii-1)*fstep;   %start freq of current wide mesh
    f4 = fmin+(ii)*fstep;     %stop freq of current wide mesh
    med_freq = f3:fmed:f4;
    lennnn = (1:length(med_freq))-1; med_freq = f3+lennnn*fmed;          %<-----

    %this is defining the temporary coarse freq resolution array
    f5 = fmin+(ii-1)*fstep;   %start freq of current wide mesh
    f6 = fmin+(ii)*fstep;     %stop freq of current wide mesh
    coarse_freq = f3:fcor:f4;
    lennnn = (1:length(coarse_freq))-1; coarse_freq = f3+lennnn*fcor;    %<-----

    %%%%%%%%% define the current output resolution wavenumber grids %%%%%%%
    %at the output wavenumber resolution, this is the current wavenumber array
    output_freq = f3:nbox*ffin:f4-nbox*ffin;
    lennnn  =  (1:length(output_freq))-1; coarse_freq = f3+lennnn*nbox*ffin;  %<--
    %occupying the following columns in the output array outwave,out_array
    thechunk = chunk+(ii-1)*chunklen;

    %%%%%%%% sort the lines into very near, medium near, far %%%%%%%%%%%%%%%
    [very_near] = sortbins(line,f3,f4,-1,xnear);
    [medium_near] = sortbins(line,f3,f4,xnear,xmed);
    [far_wing] = sortbins(line,f3,f4,xmed,xfar);

    fprintf(1,'gasid %3i mesh %3i start stop freq = %8.5f %8.5f num of lines VeryNear/MediumNear/FarWing = %5i %5i %5i ',gasID,ii,f3,f4,very_near.linct,medium_near.linct,far_wing.linct);
%    fprintf(1,'gasid %3i mesh %3i start stop freq = %8.5f %8.5f \n',gasID,ii,f3,f4);
%    fprintf(1,'  number of fine/medium/coarse lines = %5i %5i %5i \n',very_near.linct,medium_near.linct,far_wing.linct)
 
    for jj = MinLayer:Step:MaxLayer 
      %INNER LOOP OVER LAYERS 1..100 = bottom -> top
      nn = (jj-MinLayer)/Step+1;

      if mod(jj,10) == 0
        fprintf(1,'.')
      end

      tempr     = temperature(jj);
      p         = press(jj);
      ps        = partpress(jj);
      qisotopes = getq_oldVSnew(hlist_qtips,hitran_version,A,B,C,D,G,...
                                mass_info,gasID,tempr);

%qfcn = qoldVSqnew_fast(hlist_qtips,hitran_version,A,B,C,D,G,...
%            qisotopes,mass_info,gasID,line,tempr);
%save co_lines.mat line qfcn
%error('jfs')

      %%%%%%%%%%%%%%% do the FINE GRID first (very near points)
      %fprintf(1,'\n number of VERY NEAR lines = %3i',very_near.linct)
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

%[freq; s0]
%very_near.iso
%mass_iso
%donk = find(line.wnum > min(freq)-1 & line.wnum < max(freq)+1);
%semilogy(line.wnum,line.stren,...
%         line.wnum(donk),line.stren(donk),'bo',freq,s0,'rx'); grid

%tempr
%very_near.linct
%length(fine_freq)
%lvgNum

%disp('this is testing MATLAB version of loop.mexa64')
%        xoutvect = xloop(very_near.iso,mass_iso,brd,strength,freq,fine_freq,...
%                   tempr,very_near.linct,length(fine_freq),lvgNum);
%error('kjsgkljsg')

        outvect = loop(very_near.iso,mass_iso,brd,strength,freq,fine_freq,...
                   tempr,very_near.linct,length(fine_freq),lvgNum);

% semilogy(fine_freq,outvect,'o-'); hold on
% pause(0.1)

        tempjunk    = boxint2(outvect,nbox);
        out_array(nn,thechunk) = out_array(nn,thechunk)+tempjunk;
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

        z = loop(medium_near.iso,mass_iso,brd,strength,freq,med_freq,...
                    tempr,medium_near.linct,length(med_freq),lvgNum);

        %do a spline interpolation of med_freq onto output_freq 
        xspline = spline(med_freq,z,output_freq);      %linear not good enough
        out_array(nn,thechunk) = out_array(nn,thechunk)+xspline;

      end

      %%%%%%%%%%%%%%% do the COARSE GRID third  (farpoints)
      %fprintf(1,'\n number of FAR lines = %3i',far_wing.linct);

      clear y
      z=zeros(size(coarse_freq));

      if far_wing.linct > 0
        qfcn      =  qoldVSqnew_fast(hlist_qtips,hitran_version,A,B,C,D,G,...
                                   qisotopes,mass_info,gasID,far_wing,tempr);
        pwr      = far_wing.abcoef;
        for_brd  = far_wing.abroad;
        self_brd = far_wing.sbroad;
        freq     = far_wing.wnum+press(jj)*far_wing.tsp;  %freq shift
        energy   = far_wing.els;
        s0       = far_wing.stren;
        brd      = broad(p,ps,1.0,for_brd,self_brd,pwr,tempr,gasID);
        strength = find_stren(qfcn,freq,tempr,energy,s0,GasAmt(jj));

        z = loop(far_wing.iso,mass_iso,brd,strength,freq,coarse_freq,...
                   tempr,far_wing.linct,length(coarse_freq),lvgNum);
        %do a spline interpolation of coarse_freq onto output_freq 

        if (length(coarse_freq) == 2)  
          xspline = interp1(coarse_freq,z,output_freq);    %linear good enough 
        elseif (length(coarse_freq) == 3)  
          xspline = interp_quad(coarse_freq,z,output_freq);  %quad good enough 
        elseif (length(coarse_freq) > 3)  
          xspline = spline(coarse_freq,z,output_freq);  %cubic good enough 
        end 

        out_array(nn,thechunk) = out_array(nn,thechunk)+xspline;
      end

    end           %inner loop over layers
    fprintf(1,'\n');
  end             %outer loop over wide meshes (frequencies)
end               %if number of lines > 0

%************************************************************************
%now do continuum .. it should be smooth enough that no need to do at
%fine mesh ... just do at output mesh
%       SUBROUTINE CALCON23( IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P, 
%     $    PARTP, AMNT, CON, CKD, whichlayer) 
%                    has been changed to
%      con = CALCON23( IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P, 
%     $    PARTP, AMNT, CKD, whichlayer)


docontinuum = -1;

if (gasID == 1)         %water
  if ((CKD == 0) | (CKD == 21) | (CKD == 23) | (CKD == 24))
    docontinuum = -1;
    fprintf(1,'\n this code WILL NOT DO water continuum ... \n');
    fprintf(1,'\n use run7water/ru76watercontinuum instead ... \n');
    error('exiting out of run7');
  end
end

if (gasID == 7)         %o2
  if (CKD == 1)
    docontinuum = 1;
    fprintf(1,'\n doing O2 continuum ... \n');
  end
end

if (gasID == 22)         %n2
  if (CKD == 1)
    docontinuum = 1;
    fprintf(1,'\n doing N2 continuum ... \n');
  end
end

%%%%%%%%%%%%%%%%% this is the orig code .... might blow up if input vectors
%%%%%%%%%%%%%%%%  are tooo large
%if (docontinuum > 0)
%    for jj=MinLayer:Step:MaxLayer %INNER LOOP 1..100 = bottom -> top
%      nn=(jj-MinLayer)/Step+1;
%      tempjunk=calcon(gasID,length(outwave),outwave,1,length(press),...
%                    temperature,press,partpress,GasAmt,CKD,jj);
%      out_array(nn,:)=out_array(nn,:)+tempjunk;
%    end
%end      

%******************************************************

%%%%from FORTRANFILES/max.inc
%c this is max length of arrays that can be used in the Mex Files  
%c this number came out of 
%c   200000 = max number of elements in mesh 
%c              eg (755-655)/0.0005 = 160000
%c        4 = number tacked on to arrays so boxint(y,5) can be done  
%      integer MaxLen
%      parameter(MaxLen=200010)
MaxLen = 200010;

%% docontinuum  =  docontinuum*2   %%% this way, only pass out CONTINUUM

dff = ffin*nbox;

%% to figurre out V1 V2 V3 V4 go to CKDLINUX/calcon.F
%% and uncomment the relevant iVers+includefile, compile with 
%% makewatermexLinux, then 
%%    mv  calcon.mexa64 calconVN.mexa64
%%

if (docontinuum == 1 & O2O3N2continuumVers < 0)
  fprintf(1,'GasID = %3i : needs continuum \n',gasID);
  disp('but you have set O2O3N2continuumVers = -1 so NO continuum added on!')

elseif (docontinuum == 1 & O2O3N2continuumVers > 0)
  disp('N2 or O2 : adding on continuum OD')
  [mm0,nn0] = size(out_array);
  if O2O3N2continuumVers == 5 & nn0 <= 10001
     homedirx = pwd;
     cd /home/sergio//SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/
     [fr5,od5] = wrapper_run8_gas7_gas22_continuum(gasID,fmin,fmax,0.0005,profname);
     out_array = out_array + od5;
     cder = ['cd ' homedirx];
     eval(cder);
  elseif O2O3N2continuumVers == 5 & nn0 > 10001
     addpath /home.sergio//SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/
     error('ooops need to increase array sizes for wrapper_run8_gas7_gas22_continuum MTCKD3.2')
     [fr5,od5] = wrapper_run8_gas7_gas22_continuum(gasID,fmin,fmax,0.0005,profname);
     out_array = out_array + od5;
  elseif (O2O3N2continuumVers ~= 5) & (length(outwave) <= MaxLen)          %can do it in one big chunk!
    for jj=MinLayer:Step:MaxLayer %INNER LOOP 1..100 = bottom -> top
      nn=(jj-MinLayer)/Step+1;
      switch O2O3N2continuumVers
        case 1
          tempjunk = calconV1(gasID,length(outwave),outwave,1,length(press),...
                  temperature,press,partpress,GasAmt,CKD,jj);
        case 3	  
          tempjunk = calconV3(gasID,length(outwave),outwave,1,length(press),...
                  temperature,press,partpress,GasAmt,CKD,jj);
        case 4
          tempjunk = calconV4(gasID,length(outwave),outwave,1,length(press),...
                  temperature,press,partpress,GasAmt,CKD,jj);
        case 5
          tempjunk = calconV_MTCKD3p2(gasID,length(outwave),outwave,1,length(press),...
                  temperature,press,partpress,GasAmt,CKD,jj);
        otherwise
          error('only have o2/n2/o3 continuum Vers 1,3,4,5');
      end
      out_array(nn,:) = out_array(nn,:)+tempjunk;
    end     %% for loop

  elseif O2O3N2continuumVers ~= 5 & (length(outwave) > MaxLen) %have to break it into frequency intervals
    for jj = MinLayer:Step:MaxLayer %LOOP 1..100 = bottom -> top
      nn = (jj-MinLayer)/Step+1;
      index0 = 1:10000;
      number = floor(length(outwave)/10000);

      for kk = 1:number             %LOOP OVER FREQUENCY INTERVALS
        index = index0+(kk-1)*10000;
        switch O2O3N2continuumVers
          case 1
            tempjunk = calconV1(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 3
            tempjunk = calconV3(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 4
            tempjunk = calconV4(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 5
            tempjunk = calconV_MTCKD3p2(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          otherwise
            error('only have o2/n2/o3 continuum Vers 1,3,4,5');
      end

        out_array(nn,index) = out_array(nn,index)+tempjunk;
      end

      %%see if anything left over
      if (index(length(index))+1 <= length(outwave)) 
        index = index(length(index))+1:length(outwave);
        switch O2O3N2continuumVers
          case 1
            tempjunk = calconV1(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 3
            tempjunk = calconV3(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 4
            tempjunk = calconV4(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 5
            tempjunk = calconV_MTCKD3p2(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          otherwise
            error('only have o2/n2/o3 continuum Vers 1,3,4,5');
      end
        out_array(nn,index) = out_array(nn,index)+tempjunk;
      end

    end
  end            %%%%if then else loop

elseif (docontinuum == 2)
  disp('N2 or O2 : replacing OD with continuum OD')
  error('why bother with this???????')
end

%******************************************************
cpptotal = cputime-cpptotal;
semilogy(outwave,out_array);

%toc
%toc0 = toc;
%fprintf(1,'toc-tic = %12.6f mins    cpptotal = %12.6f mins \n',(toc0-tic0)/60,cpptotal/60)
fprintf(1,'cpptotal = %12.6f mins \n',cpptotal/60)
