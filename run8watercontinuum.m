function [outwave,out_array] = ...
          run8watercontinuum(gasID,fmin,fmax,profname,topts);

%%% NEW : includes the multiplier files if asking only for self 
%%% or forn coefficients
%%%       so that the CKD is CKD .* mult

%%% CKD 0,21,23,24 are Tony Clough's continuums from LBLRTM
%%% CKD 1          is  Tobin/Clough latest water continuums from LBLRTM
%%% CKD 2          is  Tobin/Clough latest water continuums from LBLRTM
%%%                    with Scott's fixes (Apr 2003)
%%%               /asl/packages/sartaV104/Src_tuning/tunmlt_con1_14jul03.txt
%%%
%%% CKD 3          is  Tobin/Clough latest water continuums from LBLRTM
%%%                    with Scott's fixes (Jan 2004) 
%%%        /asl/data/sarta_database/Data_jan04untun/Coef/tunmlt_jan04deliv.txt
%%% CKD 4          is  Tobin/Clough latest water continuums from LBLRTM
%%%                    with Scott's fixes (Apr 2003) applied to SELF + FOREIGN
%%%                    as opposed to CKD2, which only applies to SELF
%%%               /asl/packages/sartaV104/Src_tuning/tunmlt_con1_14jul03.txt
%%%
%%% CKD 5          is  CKD1 with Scott's fixes (Apr 2003) applied to SELF
%%%           /home/sergio/SPECTRA/CKDLINUX/tunmlt_iasi_may09X.dat'; %% for v6
%%% CKD 6          is  CKD1 with Scott's fixes (Apr 2003) applied to SELF + FOREIGN
%%%                    as opposed to CKD5, which only applies to SELF
%%%           /home/sergio/SPECTRA/CKDLINUX/tunmlt_iasi_may09X.dat'; %% for v6
%%%  CKD 25,32,43 latest MT-CKD  (v2.5,v3.2,v4.3)

%%%%  same as run6watercontinuum, except we have {topts} = vargin
% [outwave,out_array]=run6watercontinuum(gasID,fmin,fmax,...
%       ffin,fmed,fcor,fstep,...
%       xnear,xmed,xfar,nbox,strength_far,strength_near,divide,CKD,...
%       selfmult,formult,usetoth,local,profname);
%%%% has been replaced by
% [fr,k] = run7watercontinuum(gasID,fmin,fmax,profname,{ffin,nbox,...
%                divide,CKD,selfmult,formult,local});
% [fr,k] = run7watercontinuum(gasID,fmin,fmax,profname,{topts})
%
% this is really a code to just compute the self and foreign parts of continuum
% it DOES NOT LOOP OVER lines in the database. but, depending on divide,
% it divides out (a) by 1.0 (leaves you entire continuum computation)
%                    q v tanh(c2 v/2T) (296/T) * (ps CS + pf CF)
%                (b) q v tanh(c2 v/2T) (296/T) * ps CS      or
%                    q v tanh(c2 v/2T) (296/T) * (p-ps) CF  or
%                    q v tanh(c2 v/2T) (296/T) 
%                    depending on selfmult,formult values
%
% divide   self    for    divide by                              result
% ----------------------------------------------------------------------
%   -1     ?       ?          1               qvtanh(c2v/2T)(296/T)(psCS+pfCF) = total ODcon
%    1     1       0     q v tanh(c2 v/2T) (296/T) * ps           CS
%    1     0       1     q v tanh(c2 v/2T) (296/T) * (p-ps)       CF
%    1     ?       ?     q v tanh(c2 v/2T) (296/T)                ps CS + pf CF
%
% for example, to plot out just the self coeffs for CKD2.1, do
% self.divide=1; self.selfmult=1; self.formult=0; self.CKD=21;
% [w,k21self]=run7watercontinuum(1,1580,1630,'IPFILES/co2one',self);

%
% recall output spacing = dff=ffin*nbox;
% might as well straight away use 0.0025 res with no 5 pt averaging  
% [fr,k]=run6watercontinuum(1,500,3000,0.0025,0.1,0.5,1,1,2,25,1,0.0,0.0,,...
%            divide,ckdvers,A,B,1,C,profname);
% [fr,k]=run6watercontinuum(1,500,3000,0.0025,0.1,0.5,1,1,2,25,1,0,0,...
%             -1,24,1,1,2000,0,'profname');
% becomes [fr,k] = run7watercontinuum(1,500,3000,profname)

% this function computes the lines spectra for gasID = 1 = WATER
%
%  TYPE     VAR           DESCRIPTION              TYPICAL VALUE
%----------------------------------------------------------------------
%            these need to be passed in by the user
%------------------------------------------------------------------------
% integer  useruser        interactive(+1) or file driven
%                          default = -1 = file driven      -1
% integer   gasID          HITRAN gas ID                   1
%
% integer   fmin           minimum freq (cm-1)            605
% integer   fmax           maximum freq (cm-1)            630
%
% matrix   profname        this contains N x 5 matrix : layernum,p,pp,t,x
%                          where column 1 = layer number 
%                                column 2 = layer pressure (atm)
%                                column 3 = layer partial pressure (atm)
%                                column 4 = layer temperature (k)
%                                column 5 = layer gas amnt (kilomolecules/cm2)
%
%------------------------------------------------------------------------
%  TYPE     VAR           DESCRIPTION                   DEFAULT VALUE
%            these are default values; can be reset using structure {topts}
%            most are useless values, from run7water.m
%------------------------------------------------------------------------
%
% real      ffin           fine point spacing (cm-1)      0.0005   USEFUL
% integer   nbox           boxcar sum size (odd integer)    5      USEFUL
%
%------------------------------------------------------------------------
% integer   divide         divide out by some values        -1 for no division
%                                                           +1 for yes division
%                          if -1 gives kc = v tanh(c2 v/2T)(296/T)(psCS + pfCF)
%                          if +1 can give you CS or CF, depending on values of
%                           selfmult,formult
% divide   self    for    divide by                              result
%   -1     ?       ?          1               qvtanh(c2v/2T)(296/T)(psCS+pfCF)
%    1     1       0     q v tanh(c2 v/2T) (296/T) * ps           CS
%    1     0       1     q v tanh(c2 v/2T) (296/T) * (p-ps)       CF
%    1     ?       ?     q v tanh(c2 v/2T) (296/T)                ps CS + pf CF

% integer   CKD            continumm no (-1)                      25
%                          yes water : (0,21,23,24,51,55,1,5,6,25,32 ....) 
%                          yes O2,N2 : (+1) 
%
% real     selfmult        multiplier for self part of contiuum  0<x<1
%          formult         multiplier for for  part of contiuum  0<x<1
%
% integer  local           use local lineshape                       0
%                                   +1 local defn with chi
%                                    0 local defn w/o chi
%                                   -1 to use run6 defn
% integer iLayDo          after reading in profile,             -1
%                         maybe you only want layer J           
%                         else -1 uses all layers
%
%------------------------------------------------------------------------
%
%restrictions :
%(1) xnear <= xmed <= xfar        
%(2) xnear >= fstep
%(3) xmed/fmed  xnear/fmin   fstep/fmed   fstep/ffin       are integers
%(4)fstep/(nbox*ffin)        fcor/ffin                     are integers
%(5)(fmax-fmin)/fstep        (fmax-fmin)/fcor              are integers
%
%thus ffin,fmed,fcor are all 10/2^n   n <= 5

% HISTORY : 
% 14) incluing multiplier files
% 13) really reduced the parameter list!
% 12) changed it to run7(gasID,v1,v2,profile,{topts})          Feb 25,2001
% 11) replaced UseToth   with HITRAN                           Feb 23,2001
% 10) included the local lineshape defn as another parameter
% 9) profile file name, which contains 5 column atmosphere profile
% 8) continuum for water
% 7) allows 3 different line shapes : lorentz, voigt or van huber
% 6) differsn from run4.m in that the loop over lines is now MEXED!!!
% 5) allows the user to choose ALL isotopes of a gas (0) or one of the
%    isotopes
% 4) this differs from run3.m in that the code computes the optical depths
%    of all lines in EACH layer (k*length), and then does a UNION of all the
%    lines, so that lines do not turn on and off, driving the SVD crazy
% 3) we have parameter xnear xmed xfar to define fine,med,coarse meshes
% 2) better version of run2.m (vectorised as much as possible)
% 1) same as run1WORKS except that this should compute the spectra for all 100
%    layers (run1WORKS is mainly a debugger for ONE layer)
%
% [fr,k] = run6watercontinuum(gasID,fmin,fmax,{ffin,fmed,fcor,fstep,...
%                xnear,xmed,xfar,nbox,strength_far,strength_near,,...
%                divide,CKD,selfmult,formult,HITRAN,local},profname);

%% before for 605-2830 cm-1, we have default ffin = 0.0005 cm-1, nbox = 5
%% and then ffin*nbox = 0.0025 cm-1;
%% since things are so smooth, forget the 5 point averaging and just do
useruser      = -1;
ffin          = 0.0025;           % <--- note this ---------->
nbox          = 1;                % <--- note this ---------->
iLayDo        = -1;               % do all layers
nbox0 = 5;

newparams = runXtopts_params_smart(fmin);  
  ffin  = newparams.ffin;  
  fmed  = newparams.fmed;  
  fcor  = newparams.fcor;  
  fstep = newparams.fstep;  
  xnear = newparams.xnear;  
  xmed  = newparams.xmed;  
  xfar  = newparams.xfar;  
ffin = nbox0*ffin;
nbox = 1;
          % <--- note this ---------->
divide        = -1;
CKD           = 1;        
selfmult      = 1.0;
formult       = 1.0;
local         = 0;
multfile      = 'NONE';    %% this generically replaces FNAMESELF && FNAMEFORN

allowedparams = [{'ffin'},     {'nbox'},         {'CKD'},       ...
                 {'divide'},   {'selfmult'},     {'formult'},   ...
                 {'local'},    {'multfile'},     {'iLayDo'},    ...
		 {'useruser'}];

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
     fprintf(1,'topts param not in allowed list ... \n');
     error('quitting run8watercontinuum');
   end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wherever code says %%%%%%%%%%%%%% USER_WARN the user can reset some params

%%%%%%%%%% this uses MEX files for 
%%%%%%%% (1) boxint2        (boxint.f is the slow MATLAB file)
%%%%%%%% (2) vhh,vhh1       (vhh is faster than tobin's vhh1)

if (abs(divide) ~= 1)
  error('to compute CS,CF, this code needs paramater divide = +/- 1 only!')
end

if (gasID ~= 1) 
  error('GasID must be 1 : this version is for water only!')
end

fprintf(1,'you have asked for CKD version %8.6f \n',CKD);

%% this was old allowed list
%% if (~ismember(CKD,[0 1 2 3 4 5 6 21 23 24 50 51 52 53 55 56 60]))
%%   disp('valid CKD vers = [1 2 3 4 5 6 - 0 21 23 24 -50 51 52 53 55 56 60]');
%%   error('invalid CKD version! please retry!')
%% end

%% this is new list
%% 0,21,23,24  are the 1990s CKD
%% 1,25,27,32  are new MT-CKD
%% 4 6         are derived from MT-CKD1 (cant remember how to derived 2,3,5)
%%             are derived from MT-CKD25

%origCKD = [0 21 23 24];
%MTCKD1  = [ [1 ] [4 6]];
%MTCKD25 = [ [25  27 32 43]];

valid_CKD_versions

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% see /home/sergio/CKDLINUX/ckd_lookupBIN_v4_ieee_le.m
%% see /home/sergio/CKDLINUX/ckd_lookupBIN_v5_ieee_le.m
%% see /home/sergio/CKDLINUX/ckd_lookupBIN_v6_ieee_le.m
FNAMESELF = 'NONE';
FNAMEFORN = 'NONE';
if CKD == 4
  %% for v4
  FNAMESELF = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_jan04deliv.dat';  
  FNAMEFORN = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_jan04deliv.dat';
elseif CKD == 6
  %% for v6 .. note from high spectral IASI
  FNAMESELF = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_iasi_may09X.dat'; 
  FNAMEFORN = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_iasi_may09X.dat';
end

if ~strcmp(multfile,'NONE')
  FNAMESELF = multfile;
  FNAMEFORN = multfile;
end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%%%%%%%%%%%%%%%%%%% IS THIS INTERACTIVE SESSION OR CLUNK THRU PROFILE %%%%%%%%%
%useruser=-1;
if (useruser > 0)
  which_isotope=input('Enter which isotope : ');

  do_load=0;
  MinLayer=1; MaxLayer=1; Step=1;        %use only ONE "layer"
  NumLayers=1;
  TheLayers=MinLayer:Step:MaxLayer;

  MGC=8.314674269981136  ;  
  press=input('Enter total pressure (in atm) : ');
  partpress=input('Enter gas partial pressure (in atm) : ');
  temperature=input('Enter temperature (in K) : ');
  GasAmt=input('Enter path cell length (in cm) ');
  %change to kilomolcles cm-2 
  GasAmt=GasAmt*101325*partpress/1e9/MGC/temperature; %change to kmolec/cm2 
end

%%%%%%%%%%%%%%%%%%% LOAD IN GAS PROFILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if useruser < 0
  load_gas_profile
end

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

which_isotope=0;                                %default use all isotopes
MinLayer=1; MaxLayer=length(GasAmt); Step=1;    %default step thru ALL layers
NumLayers=(MaxLayer-MinLayer+Step)/Step;
TheLayers=MinLayer:Step:MaxLayer;

%%%%%%%%%%%%%%%%%% DEFINE OUTPUT ARRAYS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define output freq array, output 100 layer matrix
%output point spacing grid
output_spacing=nbox*ffin;
%eg 605:0.0025:629.9975;
fmaxUSE   = fmax-output_spacing;
outwave   = fmin:output_spacing:fmaxUSE; 
lennnn = (1:length(outwave))-1; outwave=fmin+lennnn*output_spacing;  %<-----
out_array = zeros(NumLayers,length(outwave));

cpptotal=cputime; 

%%%%%%%%%%%%%%%%%%%%%%%%%% continuum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now do continuum .. it should be smooth enough that no need to do at
%fine mesh ... just do at output mesh

%%%%from FORTRANFILES/max.inc
%c this is max length of arrays that can be used in the Mex Files  
%c this number came out of 
%c   200000 = max number of elements in mesh 
%c              eg (755-655)/0.0005 = 160000
%c        4 = number tacked on to arrays so boxint(y,5) can be done  
%      integer MaxLen
%      parameter(MaxLen=200010)
MaxLen = 200010;

% divide   self    for    divide by                              result
%   -1     ?       ?          1               qvtanh(c2v/2T)(296/T)(psCS+pfCF)
%    1     1       0     q v tanh(c2 v/2T) (296/T) * ps           CS
%    1     0       1     q v tanh(c2 v/2T) (296/T) * (p-ps)       CF
%    1     ?       ?     q v tanh(c2 v/2T) (296/T)                ps CS + pf CF

%%%%%   divide out   [q v tanh(c2 v/2T) (296/T)]
c2 = 1.4387863;
AVOG = 6.022045E+26;
 
CKD_0 = CKD;

if (ismember(CKD_0,[50 51 52 53 55 56 60]))
  % use F77 MEX code to do the < 1300, > 1800 region, then Matlab code to 
  % do the 1300-1800 region
  CKD = 24;          
end

if (ismember(CKD_0,[1 5 6])) & ...
  ((divide == +1) & (min(outwave) <= 2830) & (max(outwave) >= 605))
  disp('  --> only doing CKD coeffs!!!')
  disp('  --> so use CKD1 and then multiply by relevant multiplier')
  CKD = 1;          
end

dff=ffin*nbox;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ((local == 0)|(local==1))          
  %local lineshape for water
  fprintf(1,'CKD = %3i CKD_0 = %3i local = %2i \n',CKD,CKD_0,local)
  out_array = do_local_lineshape_CKD(outwave,out_array,AVOG,c2,...
                   temperature,press,partpress,GasAmt,...
                   CKD,CKD_0,selfmult,formult,profname,...
                   local,divide,ffin,nbox,MaxLen,MinLayer,Step,MaxLayer);
elseif (local == -1)          
  %lorentz lineshape for water
  do_lorentz_lineshape_CKD
end      

[junkx,junky] = size(out_array);
for ii = 1 : junkx
  bad = find(out_array(ii,:) < eps);
  %if length(bad) > 0
  if length(bad) > 0 & CKD < 43    
    fprintf(1,'oh oh : found %6i bad points out of %6i in row %3i of out_array, CKD = %8.4f \n',length(bad),length(outwave),ii,CKD);
    good = setdiff(1:length(outwave),bad);
    out_array(ii,bad) = interp1(outwave(good),out_array(ii,good),outwave(bad),[],'extrap');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[junkx,junky] = size(out_array);

%% this is ONLY for CKD 4 and 6, since we multiply MT-CKD 1 by Scott adjustments
if ((divide == +1) & (min(outwave) <= 2830) & (max(outwave) >= 605))
  figure(2); clf; 
  if ((selfmult >= 0.999999) & (formult <= 0.00000001) & ~strcmp(FNAMESELF,'NONE'))
    disp('MULTIPLY SELF COEFFS!!')
    xmult = load(FNAMESELF);
    len = length(xmult);
    multiplier = interp1(xmult(:,2),xmult(:,5),outwave); 
    oo = find(outwave <= xmult(1,2) | outwave >= xmult(len,2));
     multiplier(oo) = 1.0;
    plot(outwave,multiplier,'o-')
    figure(3); plot(xmult(:,2),xmult(:,5))
    multiplier = ones(junkx,1)*multiplier;
    out_array = out_array .* multiplier;
  end
  if ((formult >= 0.999999) & (selfmult <= 0.00000001) & ~strcmp(FNAMEFORN,'NONE'))
    disp('MULTIPLY FORN COEFFS!!')
    xmult = load(FNAMEFORN);
    len = length(xmult);
    multiplier = interp1(xmult(:,2),xmult(:,5),outwave);
    oo = find(outwave <= xmult(1,2) | outwave >= xmult(len,2));
    multiplier(oo) = 1.0;
    plot(outwave,multiplier,'o-')
    multiplier = ones(junkx,1)*multiplier;
    out_array = out_array .* multiplier;
  end
end

cpptotal = cputime-cpptotal;
fprintf(1,'cpptotal = %8.6f \n',cpptotal);

%cd /strowdata1/shared/sergio/HITRAN2UMBCLBL/MAKE_CKD
%error('kljkaf')
