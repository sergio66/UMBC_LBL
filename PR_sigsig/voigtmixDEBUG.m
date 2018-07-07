function [voimix]=voigtmix(freqq,f,ymix,jq,temperature,w_forq,w_selfq,... 
                     strenqt,stuff,amt,NIF,birn); 

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

%c alpha_doppler=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light)
%c r2=2*log(2)*k/amu
      r2=11526.218;

K_voigt_mixing=density*pressure_self/pressure_ref*temperature_ref*... 
        path_length/temperature  %%%%%%%NOTE no factor of pi 
K_lorentz_mixing=density*pressure_self/pressure_ref*temperature_ref*... 
        path_length/temperature/pi  %%%%%%%NOTE factor of pi 
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

  stuff.beta 
  stuff.duration 
  plot(freqq); title('freqr'); pause 
  plot(freqq,ymix); title('ymix'); pause 
  plot(freqq,jq); title('jr'); pause 
  plot(freqq,w_forq); title('forr'); pause 
  plot(freqq,w_selfq); title('selfr'); pause 
  plot(freqq,strenqt); title('strenrt'); pause 
   
%w_tot=(pressure_self*w_selfq+pressure_for*w_forq)/pressure_ref; 
%plot(freqq,w_tot); title('braodening'); pause;

for i=1:length(freqq)
  %this is the broadening
  w_tot=(pressure_self*w_selfq(i)+pressure_for*	w_forq(i))/pressure_ref; 

  %   alpha_doppler=freqq(i)/c_light*sqrt(r2*temperature/mass);
  %   [freqq(i) strenqt(i)*K_voigt_mixing w_tot alpha_doppler]

  if ((NIF == 'V') | (NIF == 'v'))
    [wr,wi]=vhh1RI(f',freqq(i),temperature,mass,w_tot);
    wr=strenqt(i)*wr;
    wi=strenqt(i)*wi;
    temp=(wr'+ymix(i)*wi');  %%ymix=zeros if IO set to 0
    %no need to do the f/freqq(i)*tanh(f/T)/tanh(freqq(i)/T) because this
    %is already done by the vhh mex file

%   this is to use the new VOIGTER that GENLN2 uses
%   [wr2,wi2]=vhh2RI(f',freqq(i),temperature,mass,w_tot);
%   wr2=strenqt(i)*wr2;
%   wi2=strenqt(i)*wi2;
%   voi2=voi2+wr2';
%   voimix2=voimix2+(wr2'+ymix(i)*wi2'); 
%   %no need to do the f/freqq(i)*tanh(f/T)/tanh(freqq(i)/T) because this
%   %is already done by the vhh mex file

  elseif ((NIF == 'L') | (NIF == 'l'))
    temp=strenqt(i)*(f/freqq_shift(i)).*(w_tot+...
       (f-freqq_shift(i))*ymix(i))./((f-freqq_shift(i)).^2+(w_tot).^2);
    %%% no need to do 
    %%%temp=temp.*f/freqq(i).*tanh(f*c2/2/temperature)./ ...
    %%%                    tanh(freqq(i)*c2/2/temperature);
    %%% as tanh ~ 1, and we already do f/freqq(i)
  else
    error('need to do either voigt or lorentz')
    end
  
  % compute chi function 
  if ((birn == 'Y') | (birn == 'y'))
    Am=1; 
    ii=sqrt(-1); 
    tau0=0.72/temperature; 
    tau2=duration; 
    dnu=f-freqq_shift(i); 
    zz=sqrt((w_tot^2 + dnu.^2).*(tau0^2 + tau2^2)); 
    ex=exp(tau2*w_tot + tau0*dnu); 
    chi=Am.*zz.*(-pi/2).*real(besselh(1,ii*zz)).*ex; 
  else
    chi=ones(size(temp));
    end

  kkk=(temp.*chi >= 0.0);
  [freqq(i) dnu(1) dnu(length(dnu)) chi(1) chi(length(chi)) ...
     zz(1) zz(length(zz)) ex(1) ex(length(ex))  strenqt(i) ymix(i) w_tot ...
     temp(1) temp(length(temp))]
  figure(1)
  subplot(211): plot(f,temp.*chi.*kkk*K_lorentz_mixing);
  subplot(212); plot(f,chi); 
  voimix=voimix+temp.*chi;
  figure(2)
  plot(f,voimix*K_lorentz_mixing); pause(0.1)
  end

tau0
tau2
birn
stuff.beta
duration
%%voi=voi*K_voigt_mixing*bsm;

if ((NIF == 'V') | (NIF == 'v'))
  voimix=voimix*K_voigt_mixing*bsm;
else
  voimix=voimix*K_lorentz_mixing*bsm;
  end

