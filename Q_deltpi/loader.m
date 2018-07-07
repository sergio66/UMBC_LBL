function [jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff]=...
     loader(temperature,band,path_length,ptotal,pself)

stuff.band=band;

stuff.p1=2.0649e-02; 
stuff.p2=1.5840e-01; 

stuff.path_length=path_length;       %path length in cm
stuff.pressure_self=pself;           %already in atm 
stuff.pressure_for=ptotal-pself;     %already in atm 

%/sals/scratch4/Strow/Tobin_home/tobin/Co2q/B_deltpie/Mfiles/Fit/nfive400.mat
%this corresponds to file 3 in loader-mod in  
%/salsify/data/Strow/Tobin/Co2q_B_deltpie
%beta_pi_self   = 0.5407505;  % fitted from pure Pi-Sigma data
%beta_pi_air    = 0.5995593;  % fitted from air-broadened Pi-Sigma data
%beta_delt_self = 0.4880067;  % fitted from pure Pi-Delta data
%beta_delt_air  = 0.72172644;  % fitted from air-broadened Pi-Delta data

[duration_pure,duration_for,beta_pi_self,beta_pi_air,...
    beta_delt_self,beta_delt_air] = co2_param(band,ptotal,pself);

stuff.beta_pi_self   = beta_pi_self;
stuff.beta_pi_air    = beta_pi_air;
stuff.beta_delt_self = beta_delt_self;
stuff.beta_delt_air  = beta_delt_air;

stuff.bsm = 1;
stuff.btz = 1.4387863;
stuff.B0  = 0.4;
% temperature_ref must go along with the value of density_ref (at 273.15 k)
stuff.temperature_ref = 273.15;
stuff.speed_light     = 2.99792458e8;
stuff.pressure_ref    = 1;
stuff.density         = 2.6867e19;
stuff.Boltzmann       = 1.380662e-23;
stuff.mass_proton     = 1.6726485e-27;
stuff.mass_CO2        = 44*stuff.mass_proton;

stuff.duration = 0.000000000000000; %just a dummy for this Q branch
duration = (pself*duration_pure+(ptotal-pself)*duration_for)/ptotal;
stuff.duration = duration;

stuff.frequency_shift=0.00000000; %just a dummy for this Q branch

if (band == 740) 
  load ../CO2_MATFILES/hit740 
elseif (band == 668) 
  load ../CO2_MATFILES/hit668 
elseif (band == 2093) 
  load ../CO2_MATFILES/hit2093 
else 
  error('wrong band in loading the hitBLAH file!!!') 
  end 
 
isotope=1;
if band==668
  v_l=2;v_u=4;
  end
if band==2093
  v_l=2;v_u=14;
  end
if band==740
  v_l=4;v_u=8;
  end

stuff.isotope=isotope;

%DAVES
index=find(j_lower(:,5)=='Q' & v_lower==v_l & v_upper==v_u & iso==isotope);

j_lowerq(:,1:9) = j_lower(index,1:9);		% J of lower state (strings)
junk_str        = j_lowerq(:,6:8);
for n = 1:length(junk_str);
  eval(['jq(n)=' junk_str(n,:) ';']);
  end

jq         = jq';
freqq      = freq(index);		% frequencies (cm-1)
strenq     = stren(index);		% line strengths
elowerq    = elower(index);		% lower state total energy
w_forq     = w(index);			% Air broadened widths(296K)
w_fortempq = w_temp(index);		% Temperature correction
					%      coefficients to air-
					%         broadened widths
w_selfq    = w_s(index);		% Self-broadened widths (296K)
% w_selftempq=0.685			% coefficients for self-widths

if jq(1)>jq(2) 
  jq         = flipud(jq);
  freqq      = flipud(freqq);
  strenq     = flipud(strenq);
  elowerq    = flipud(elowerq);
  w_selfq    = flipud(w_selfq);
  w_forq     = flipud(w_forq);
  w_fortempq = flipud(w_fortempq);
  end

tpr   = temperature;
[a1,b1,c1,d1]=qqttiippss(isotope);

Qt    = a1 + b1*tpr + c1*tpr^2 + d1*tpr^3;
Q296  = a1 + b1*296 + c1*296^2 + d1*296^3;

numer   = Q296*exp(stuff.btz*elowerq/296).*(1-exp(-stuff.btz*freqq/tpr));
denom   = Qt*exp(stuff.btz*elowerq/tpr).*(1-exp(-stuff.btz*freqq/296));
strenqt = strenq.*numer./denom;

w_forq  = w_forq.*(296/tpr).^w_fortempq;
w_selfq = w_selfq.*(296/tpr)^0.685;

[jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq]=... 
    orderer(jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq); 

stuff.population   = zeros(size(jq));
stuff.population_t = zeros(size(jq));

clear accuracy dipole elower freq gas_id iso j_lower j_upper line_status
clear p_shift reference stren v_lower v_upper w w_s w_temp j_upper index
clear hitfile v_l v_u mass_proton
