function [voimix]=voigt_mix_simple(freqq,f,ymix,jq,temperature,...
                                   w_forq,w_selfq,... 
                                  strenqt,stuff,amt,NIF,IO,birn); 
%IO   = '0' or '1' for no/yes to linemixing

%NIF  = L for lorentz
%       V  for voigt actually vanHuber
%birn = y for include birnbaum
%       n otherwise
frequency_shift=stuff.frequency_shift;
density=stuff.density; 
pressure_self=stuff.pressure_self; 
pressure_for=stuff.pressure_for; 
pressure_ref=stuff.pressure_ref; 
temperature_ref=stuff.temperature_ref; 
path_length=stuff.path_length; 
bsm=stuff.bsm; 
mass=stuff.mass_CO2;
mass=44;
duration=stuff.duration;

%c do the doppler widths first 
      k=1.380658e-23;
      c_light=2.99792458e8;        %ms-1
      amu=1.6605402e-27;            %nucleon mass/kg
      mass=44;                       %change to kg

%c alpha_doppler=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light)
%c r2=2*log(2)*k/amu
      r2=11526.218;

K_voigt_mixing=density*pressure_self/pressure_ref*temperature_ref*... 
        path_length/temperature;  %%%%%%%NOTE no factor of pi 
K_lorentz_mixing=density*pressure_self/pressure_ref*temperature_ref*... 
        path_length/temperature/pi;  %%%%%%%NOTE factor of pi 
%%%%%thus there is a factor of pi difference in defn of voigt vs lorentz
%%%%% also see Dave Tobin's voigter.m and GENLN2 manual eqn 4.8, bott pg 32
no_lines=length(freqq); no_pts=length(f); 
k=zeros(no_pts,1); 
 
frequency_shift=0; 

%%%%%%%%%%%%%%%%%% Compute Kvoi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('computing voigt absorption coefficient')
%%%

%%voi=zeros(length(f),1);
%%lor=zeros(length(f),1);
%%voi2=zeros(length(f),1);
voimix=zeros(length(f),1);
%%lormix=zeros(length(f),1);
%%voimix2=zeros(length(f),1);

c2 = 1.4387863;        %K/ cm-1  from Genln2 manual
freqq_shift=freqq+frequency_shift/100;

w_tot=(pressure_self*w_selfq+pressure_for*w_forq)/pressure_ref; 

for i=1:length(freqq)
  %this is the broadening

  %   alpha_doppler=freqq(i)/c_light*sqrt(r2*temperature/mass);
  %   [freqq(i) strenqt(i)*K_voigt_mixing w_tot(i) alpha_doppler]

  if ((NIF == 'V') | (NIF == 'v'))
    [wr,wi]=vhh1RI(f',freqq(i),temperature,mass,w_tot(i));
    wr=strenqt(i)*wr;
    wi=strenqt(i)*wi;
    temp=(wr'+ymix(i)*wi');  %%ymix=zeros if IO set to 0
    %no need to do the f/freqq(i)*tanh(f/T)/tanh(freqq(i)/T) because this
    %is already done by the vhh mex file

%   this is to use the new VOIGTER that GENLN2 uses
%   [wr2,wi2]=vhh2RI(f',freqq(i),temperature,mass,w_tot(i));
%   wr2=strenqt(i)*wr2;
%   wi2=strenqt(i)*wi2;
%   voi2=voi2+wr2';
%   voimix2=voimix2+(wr2'+ymix(i)*wi2'); 
%   %no need to do the f/freqq(i)*tanh(f/T)/tanh(freqq(i)/T) because this
%   %is already done by the vhh mex file

  elseif ((NIF == 'L') | (NIF == 'l'))
    temp=strenqt(i)*(f/freqq_shift(i)).*(w_tot(i)+...
       (f-freqq_shift(i))*ymix(i))./((f-freqq_shift(i)).^2+(w_tot(i)).^2);
    %%% no need to do 
    %%%temp=temp.*f/freqq(i).*tanh(f*c2/2/temperature)./ ...
    %%%                    tanh(freqq(i)*c2/2/temperature);
    %%% as tanh ~ 1, and we already do f/freqq(i)
  else
    error('need to do either voigt or lorentz')
    end
  
  % compute chi function 
  if ((birn == 'b') | (birn == 'B') )        %birnbaum
    %this is the broadening
    chi=birnbaum(f',freqq_shift(i),w_tot(i),temperature,duration);
    chi=chi';
  elseif ((birn == 'c') | (birn == 'C') )    %cousin
    %this is the broadening
    chi=cousin1(f',freqq_shift(i),w_tot(i),temperature,...
               pressure_for,pressure_self);
    chi=chi';
  else
    chi=ones(size(temp));
    end

  %%%%%% this is straight from DAVE TOBIN's code klormix3.m in
  %%%%%% /salsify/scratch4/Strow/Tobin_home/tobin/Co2q/B_sigpie/Run/15um
  if (IO == '1')
    ind=find(abs(f-freqq_shift(i))>15);
    temp(ind)=0.5*temp(ind);
    end

  voimix=voimix+temp.*chi;

  end

if ((NIF == 'V') | (NIF == 'v'))
  voimix=voimix*K_voigt_mixing*bsm;
else
  voimix=voimix*K_lorentz_mixing*bsm;
  end

