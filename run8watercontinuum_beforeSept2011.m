function [outwave,out_array] = ...
          run8watercontinuum(gasID,fmin,fmax,profname,topts);

%%% CKD 0,21,23,24 are Tony Clough's continuums from LBLRTM
%%% CKD 1          is  Tobin/Clough latest water continuums from LBLRTM
%%% CKD 2          is  Tobin/Clough latest water continuums from LBLRTM
%%%                    with Scott's fixes (Apr 2003)
%%%               /asl/packages/sartaV104/Src_tuning/tunmlt_con1_14jul03.txt
%%% CKD 4          is  Tobin/Clough latest water continuums from LBLRTM
%%%                    with Scott's fixes (Apr 2003) applied to SELF + FOREIGN
%%%                    as opposed to CKD2, which only applies to SELF
%%%               /asl/packages/sartaV104/Src_tuning/tunmlt_con1_14jul03.txt
%%% CKD 3          is  Tobin/Clough latest water continuums from LBLRTM
%%%                    with Scott's fixes (Jan 2004) 
%%%        /asl/data/sarta_database/Data_jan04untun/Coef/tunmlt_jan04deliv.txt
%%% CKD 51 : best from RAL data WATER/CONTINUUM/MATFILES/MAY9_BESTFIT
%%%          has CS(T) and CF(T) at T = 296, 243 K
%%% CKD 52 : best from RAL data WATER/CONTINUUM/MATFILES/MAY13
%%%          has CS(T) but CF(T) is indpt of T
%%% CKD 53 : uses latest estimates of continuum eg
%%%          Dave Tobins's x wavenumber points WATER/CONTINUUM
%%%          has CS(T) and CF(T) at T = 296, 243 K
%%% CKD 55,56 : uses latest estimates of continuum, PLUS Scott Hannon 
%%%             improvements using AIRS data (no bias with column water) eg
%%%             Dave Tobins's x wavenumber points WATER/CONTINUUM
%%%             has CS(T) and CF(T) at T = 296, 243 K

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
%   -1     ?       ?          1               qvtanh(c2v/2T)(296/T)(psCS+pfCF)
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

% integer   CKD            continumm no (-1)                   24
%                          yes water : (0,21,23,24,51,55) 
%                          yes O2,N2 : (+1) 
%
% real     selfmult        multiplier for self part of contiuum  0<x<1
%          formult         multiplier for for  part of contiuum  0<x<1
%
% integer  local           use local lineshape                       0
%                                   +1 local defn with chi
%                                    0 local defn w/o chi
%                                   -1 to use run6 defn
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
ffin          = 0.0025;           % <--- note this ---------->
nbox          = 1;                % <--- note this ---------->

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

allowedparams = [{'ffin'},     {'nbox'},         {'CKD'}, ...
                 {'divide'},   {'selfmult'},     {'formult'},   {'local'}];

%read in {topts}
if nargin == 5
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
if (~ismember(CKD,[0 1 2 3 4 5 6 21 23 24 50 51 52 53 55 56 60]))
  disp('valid CKD versions = [1 2 3 4 5 6 - 0 21 23 24-50 51 52 53 55 56 60]');
  error('invalid CKD version! please retry!')
  end
%%%%%%%%%%%%%%%%%%% IS THIS INTERACTIVE SESSION OR CLUNK THRU PROFILE %%%%%%%%%
useruser=-1;
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

load_gas_profile

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
MaxLen=200010;

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

dff=ffin*nbox;
if ((local == 0)|(local==1))          %local lineshape for water 
  if (ismember(CKD,[0 1 2 3 4 5 6 21 23 24 50]))
    fprintf(1,'\n doing water continuum for CKD %3i lineshape = %2i \n',CKD,local);
    dff=ffin*nbox;

    if (length(outwave) <= MaxLen)          %can do it in one big chunk!
      for jj=MinLayer:Step:MaxLayer %INNER LOOP 1..100 = bottom -> top
        nn=(jj-MinLayer)/Step+1;
        if (CKD_0 == 50)
          scum = mst50(outwave,temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 51)
          %%%scum = mst51(outwave,temperature,press,partpress,...
          scum = mst51_tobin(outwave,temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 52)
          scum = mst52(outwave,temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 53)
          scum = mst53(outwave,temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 55)
          scum = mst55(outwave,temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 56)
          scum = mst56(outwave,temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 60)
          scum = mst60(outwave,temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        else
          fprintf(1,'using calconwater_loc in one gulp ... ii [p T pp] CKD = %3i %8.6e %8.6e %8.6f %3i\n',jj,press(jj),partpress(jj),temperature(jj),CKD)
          scum=calconwater_loc(1,length(outwave),outwave,dff,length(press),...
          temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
          end
        if (divide  == -1)
          tempfreq=ones(size(scum));
        else
          tempfreq = AVOG*GasAmt(jj)*outwave .* ...
                     tanh(c2*outwave/2/temperature(jj))* ...
                     (296/temperature(jj));
          if ((selfmult >= 0.999999) & (formult <= 0.00000001))
            tempfreq=tempfreq * partpress(jj);
          elseif ((formult >= 0.999999) & (selfmult <= 0.00000001))
            tempfreq=tempfreq * (press(jj)-partpress(jj));
            end
          end
        if max(tempfreq) <= 1e-20
          scum=scum./(tempfreq+1e-20);   %% else e get divide by 0 = NAN
        else
          scum=scum./(tempfreq);
          end
        out_array(nn,:)=out_array(nn,:)+scum;
        end

    else %have to break it into frequency intervals
      disp('using calconwater_loc in mini gulps...')
      for jj=MinLayer:Step:MaxLayer %LOOP 1..100 = bottom -> top
        nn=(jj-MinLayer)/Step+1;
        index0=1:(MaxLen-10);
        number=floor(length(outwave)/(MaxLen-10));

        for kk=1:number             %LOOP OVER FREQUENCY INTERVALS
          index=index0+(kk-1)*(MaxLen-10);
          if (CKD_0 == 50)
            scum = mst50(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 51)
            %%%scum = mst51(outwave(index),temperature,press,partpress,...
            scum = mst51_tobin(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 52)
            scum = mst52(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 53)
            scum = mst53(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 55)
            scum = mst55(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 56)
            scum = mst56(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 60)
            scum = mst60(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 5)
            scum = tobin5(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          else
            scum=calconwater_loc(1,length(index),outwave(index),dff,...
              length(press),...
              temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
            end
          if (divide  == -1)
            tempfreq=ones(size(scum));
          else
            tempfreq = AVOG*GasAmt(jj)*outwave(index) .* ...
                       tanh(c2*outwave(index)/2/temperature(jj))* ...
                       (296/temperature(jj));
            if ((selfmult >= 0.999999) & (formult <= 0.00000001))
              tempfreq=tempfreq * partpress(jj);
            elseif ((formult >= 0.999999) & (selfmult <= 0.00000001))
              tempfreq=tempfreq * (press(jj)-partpress(jj));
              end
            end
          if max(tempfreq) <= 1e-20
            scum=scum./(tempfreq+1e-20);   %% else e get divide by 0 = NAN
          else
            scum=scum./(tempfreq);
            end
          out_array(nn,index)=out_array(nn,index)+scum;
          end

        %%see if anything left over
        if (index(length(index))+1 <= length(outwave)) 
          index=index(length(index))+1:length(outwave);
          if (CKD_0 == 50)
            scum = mst50(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 51)
            %%%scum = mst51(outwave(index),temperature,press,partpress,...
            scum = mst51_tobin(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 52)
            scum = mst52(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 53)
            scum = mst53(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 55)
            scum = mst55(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 56)
            scum = mst56(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 60)
            scum = mst60(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          elseif (CKD_0 == 5)
            scum = tobin5(outwave(index),temperature,press,partpress,...
                       GasAmt,CKD_0,selfmult,formult,jj,profname);
          else
            scum=calconwater_loc(1,length(index),outwave(index),dff,...
                length(press),...
                temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
            end
          if (divide  == -1)
            tempfreq=ones(size(scum));
          else
            tempfreq = AVOG*GasAmt(jj)*outwave(index) .* ...
                       tanh(c2*outwave(index)/2/temperature(jj))* ...
                       (296/temperature(jj));
            if ((selfmult >= 0.999999) & (formult <= 0.00000001))
              tempfreq=tempfreq * partpress(jj);
            elseif ((formult >= 0.999999) & (selfmult <= 0.00000001))
              tempfreq=tempfreq * (press(jj)-partpress(jj));
              end
            end
          if max(tempfreq) <= 1e-20
            scum=scum./(tempfreq+1e-20);   %% else e get divide by 0 = NAN
          else
            scum=scum./(tempfreq);
            end
          out_array(nn,index)=out_array(nn,index)+scum;
          end
        end
      end            %%%%if then else loop
    end
  end      

if (local == -1)          %lorentz lineshape for water
  if (ismember(CKD,[1 2 3 4 5 24 50 51 52 53 55 56 60]))
    error('v 1,2,3,4,5,24,50-56,60 cannot be used with lorentz lineshape!!');
    end
  if (ismember(CKD,[0 21 23]))
    fprintf(1,'\n doing water continuum for CKD %3i Lorentz lineshape = %2i \n',CKD,local);
    if (length(outwave) <= MaxLen)          %can do it in one big chunk!
      for jj=MinLayer:Step:MaxLayer %INNER LOOP 1..100 = bottom -> top
        nn=(jj-MinLayer)/Step+1;
        scum=calconwater(1,length(outwave),outwave,dff,length(press),...
          temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
        if (divide  == -1)
          tempfreq=ones(size(scum));
        else
          tempfreq = AVOG*GasAmt(jj)*outwave .* ...
                     tanh(c2*outwave/2/temperature(jj))* ...
                     (296/temperature(jj));
          if ((selfmult >= 0.999999) & (formult <= 0.00000001))
            tempfreq=tempfreq * partpress(jj);
          elseif ((formult >= 0.999999) & (selfmult <= 0.00000001))
            tempfreq=tempfreq * (press(jj)-partpress(jj));
            end
          end
        if max(tempfreq) <= 1e-20
          scum=scum./(tempfreq+1e-20);   %% else e get divide by 0 = NAN
        else
          scum=scum./(tempfreq);
          end
        out_array(nn,:)=out_array(nn,:)+scum;
        end
    else %have to break it into frequency intervals
      fprintf(1,'breaking into smaller intervals, length (MaxLen-10) \n');
      for jj=MinLayer:Step:MaxLayer %LOOP 1..100 = bottom -> top
        nn=(jj-MinLayer)/Step+1;
        index0=1:(MaxLen-10);
        number=floor(length(outwave)/(MaxLen-10));
        for kk=1:number             %LOOP OVER FREQUENCY INTERVALS
          index=index0+(kk-1)*(MaxLen-10);
          scum=calconwater(1,length(index),outwave(index),dff,length(press),...
              temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
          if (divide  == -1)
            tempfreq=ones(size(scum));
          else
            tempfreq = AVOG*GasAmt(jj)*outwave(index) .* ...
                       tanh(c2*outwave(index)/2/temperature(jj))* ...
                       (296/temperature(jj));
            if ((selfmult >= 0.999999) & (formult <= 0.00000001))
              tempfreq=tempfreq * partpress(jj);
            elseif ((formult >= 0.999999) & (selfmult <= 0.00000001))
              tempfreq=tempfreq * (press(jj)-partpress(jj));
              end
            end
          if max(tempfreq) <= 1e-20
            scum=scum./(tempfreq+1e-20);   %% else e get divide by 0 = NAN
          else
            scum=scum./(tempfreq);
            end
          out_array(nn,index)=out_array(nn,index)+scum;
          end

        %%see if anything left over
        if (index(length(index))+1 <= length(outwave)) 
          index=index(length(index))+1:length(outwave);
          scum=calconwater(1,length(index),outwave(index),dff,length(press),...
                temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
          if (divide  == -1)
            tempfreq=ones(size(scum));
          else
            tempfreq = AVOG*GasAmt(jj)*outwave(index) .* ...
                       tanh(c2*outwave(index)/2/temperature(jj))* ...
                       (296/temperature(jj));
            if ((selfmult >= 0.999999) & (formult <= 0.00000001))
              tempfreq=tempfreq * partpress(jj);
            elseif ((formult >= 0.999999) & (selfmult <= 0.00000001))
              tempfreq=tempfreq * (press(jj)-partpress(jj));
              end
            end
          if max(tempfreq) <= 1e-20
            scum=scum./(tempfreq+1e-20);   %% else e get divide by 0 = NAN
          else
            scum=scum./(tempfreq);
            end
          out_array(nn,index)=out_array(nn,index)+scum;
          end

        end
      end            %%%%if then else loop
    end
  end      

cpptotal=cputime-cpptotal;
fprintf(1,'cpptotal = %8.6f \n',cpptotal);
