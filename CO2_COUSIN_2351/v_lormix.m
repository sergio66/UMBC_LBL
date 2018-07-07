function [voimix]=v_lormix(iii,freqq,f,ymix,temperature,w_forq,...
                                w_selfq,strenqt,NIF,birn,ratio); 

%this is really the same as VoigtMix2.m used in run6co2 LBL code
%which is first order lorentz mixing
%big difference is that "stuff" DNE

%NIF  = L for lorentz
%       V  for voigt actually vanHuber
%birn = y for include birnbaum
%       n otherwise

global duration frequency_shift
global pressure_self 
global pressure_for 
global temperature_ref
global density Boltzmann speed_light path_length pressure_ref
global c0 cspan

amt=path_length(iii)*101325*pressure_self(iii)/1e9/8.314674269981136/...
    temperature(iii);

frequency_shift=0;

pself=pressure_self(iii); 
pfor=pressure_for(iii); 
%pressure_ref=stuff.pressure_ref; 
%temperature_ref=stuff.temperature_ref; 
plen=path_length(iii); 
bsm=1.0; 
mass=44;

%duration=stuff.duration;

%c do the doppler widths first 
      k=1.380658e-23;
      c_light=2.99792458e8;        %ms-1
      amu=1.6605402e-27;            %nucleon mass/kg
      mass=44;                       %change to kg

%c alpha_doppler=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light)
%c r2=2*log(2)*k/amu
      r2=11526.218;

K_voigt_mixing=density*pself/pressure_ref*temperature_ref*... 
        plen/temperature(iii);  %%%%%%%NOTE no factor of pi 
K_lorentz_mixing=density*pself/pressure_ref*temperature_ref*... 
        plen/temperature(iii)/pi;  %%%%%%%NOTE factor of pi 
%%%%%thus there is a factor of pi difference in defn of voigt vs lorentz
%%%%% also see Dave Tobin's voigter.m and GENLN2 manual eqn 4.8, bott pg 32
no_lines=length(freqq); no_pts=length(f); 
k=zeros(no_pts,1); 
 
%%%%%%%%%%%%%%%%%% Compute Kvoi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
voimix=zeros(length(f),1);

c2 = 1.4387863;        %K/ cm-1  from Genln2 manual
freqq_shift=freqq+frequency_shift/100;

w_tot=(pself*w_selfq+pfor*w_forq)/pressure_ref;

if ((NIF == 'V') | (NIF == 'v'))
  NIFNIF=1;
elseif ((NIF == 'L') | (NIF == 'l'))
  NIFNIF=-1;
  end
NIFNIF = -1; %always do Lorentz
NIFNIF = +1; %try VHH

birnbirn=0;
if ((birn == 'b') | (birn == 'B') ) 
   birnbirn=1;
elseif ((birn == 'c') | (birn == 'C') )  
  birnbirn=-1;
elseif ( ((birn == 'c') | (birn == 'C')) & ((NIF == 'F')  | (NIF == 'f')) )
  error('cannot have FULL and cousin!!');
  end

summm=sum(abs(ymix));  %see if the ymix coeffs are nonzero 
if ((summm > 1.0e-6) & ((birn == 'c') | (birn == 'C')) ) 
  error('cannot have mixing AND cousin!!!') 
  end 

%%orig code
%voimix=doVmix(f,freqq_shift,w_tot,temperature,duration,...
%              pfor,pself,mass,strenqt,ymix,...
%              K_voigt_mixing,K_lorentz_mixing,bsm,birnbirn,NIFNIF);

dodo = [c0 cspan temperature duration pfor pself];
voimix=doVmix2(f,freqq_shift,w_tot,temperature,duration,...
              pfor,pself,mass,strenqt,ymix,...
              K_voigt_mixing,K_lorentz_mixing,bsm,birnbirn,NIFNIF,c0,cspan);
