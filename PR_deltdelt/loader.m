function [jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr,stuff]=... 
    loader(temperature,band,path_length,ptotal,pself,prb); 

global p2311_21_jmax pr2351_jmax p2350_jmax r2350_jmax

stuff.band=band;

stuff.p1=2.0649e-02; 
stuff.p2=1.5840e-01; 

stuff.path_length=path_length;       %path length in cm
stuff.pressure_self=pself;           %already in atm 
stuff.pressure_for=ptotal-pself;     %already in atm 

if (prb == 'R')
  if (band == 2310)
    v_l = 4;
    v_u = 24;
    isotope = 1;
  elseif (band == 2311)
    v_l = 4;
    v_u = 24;
    isotope = 2;
    end
  stuff.prb = 'R';
elseif (prb == 'P')
  if (band == 2310)
    v_l = 4;
    v_u = 24;
    isotope = 1;  
  elseif (band == 2311)
    v_l = 4;
    v_u = 24; 
    isotope = 2;
    end
  stuff.prb = 'P';
else
  error('need prb ==== P or R')
  end

stuff.isotope=isotope;

%%%%%% use papaya,mango for dummy beta params
[duration_pure,duration_for,beta_pure,beta_for,papaya,mango]=...
        co2_param(band,ptotal,pself);

frequency_shift = -5.189868419218914e-01;
frequency_shift = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bsm = 1;
stuff.frequency_shift    = frequency_shift;
band_strength_multiplier = bsm; 

stuff.beta_for  = beta_for;
stuff.beta_pure = beta_pure;

beta = (pself*beta_pure+(ptotal-pself)*beta_for)/ptotal; 
if (beta > 1)
  beta = 1.0;
  end
duration = (pself*duration_pure+(ptotal-pself)*duration_for)/ptotal; 

stuff.duration = duration;
stuff.beta     = beta;

stuff.bsm  = 1;
stuff.btz  = 1.4387863;
stuff.B0   = 0.4;
% temperature_ref must go along with the value of density_ref (at 273.15 k)
stuff.temperature_ref = 273.15;
stuff.speed_light     = 2.99792458e8;
stuff.pressure_ref    = 1;
stuff.density         = 2.6867e19;
stuff.Boltzmann       = 1.380662e-23;
stuff.mass_proton     = 1.6726485e-27;
stuff.mass_CO2        = 44*stuff.mass_proton;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load /carrot/users/tobin/Co2pr/Hit_co2pr/hit_43_92.mat
%load /salsify/scratch4/Strow/Tobin_home/tobin/Co2pr/Hit_co2pr/hit_43_92.mat

if (band == 2310) 
  load ../CO2_MATFILES/hit2310
elseif (band == 2311) 
  load ../CO2_MATFILES/hit2311
else 
  error('wrong band in loading the hitBLAH file!!!') 
  end 

%%%%%%%%%%%%%%%%%%%%% Sort out lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
junk_str=j_lower(:,6:8);
for n=1:length(junk_str);
  eval(['tmp(n)=' junk_str(n,:) ';']);
        end
  tmp=tmp';

ttt=2000;       %max j
%ttt=input('enter ttt ');
if ((stuff.prb == 'P') & (length(intersect(band,[2311,2321])) > 0))
  ttt=p2311_21_jmax;
elseif ((stuff.prb == 'P') & (length(intersect(band,[2351])) > 0))
  ttt=pr2351_jmax;
elseif ((stuff.prb == 'R') & (length(intersect(band,[2351])) > 0))
  ttt=pr2351_jmax;
elseif ((stuff.prb == 'P') & (length(intersect(band,[2350,2320,2321])) > 0))
  ttt=p2350_jmax;
elseif ((stuff.prb == 'R') & (length(intersect(band,[2350])) > 0))
  ttt=r2350_jmax;
  end

%index=find(j_lower(:,5)==prb & v_lower==v_l & v_upper==v_u & iso==isotope);
index=find(j_lower(:,5)==prb & v_lower==v_l & v_upper==v_u & iso==isotope ...
            & tmp <= ttt);

% Now select out these lines from the rest.
j_lowerr(:,1:9)=j_lower(index,1:9);		% J of lower state (strings)
junk_str=j_lowerr(:,6:8);
for n=1:length(junk_str);
	eval(['jr(n)=' junk_str(n,:) ';']);
end
jr=jr';
freqr=freq(index);				% frequencies (cm-1)
strenr=stren(index);				% line strengths
elowerr=elower(index);				% lower state total energy
w_forr=w(index);				% Air broadened widths (296K)
w_fortempr=w_temp(index);			% Temperature correction
						%      coefficients to air-
						%         broadened widths
w_selfr=w_s(index);				% Self-broadened widths (296K)
% w_selftempq=0.685				% coefficients for self-widths


%%%%%%%%%%%%%%%%%% Correct widths to desired temperature %%%%%%%%%%%%%%
w_forr=w_forr.*(296/temperature).^w_fortempr;
w_selfr=w_selfr.*(296/temperature)^0.685;

%%%%%%%% Flip arrays (if needed) so that lowest J's appear first  %%%%%%%%%% 
if jr(1)>jr(2) 
   j_lowerr=flipud(j_lowerr);freqr=flipud(freqr);strenr=flipud(strenr);
   elowerr=flipud(elowerr);w_selfr=flipud(w_selfr);jr=flipud(jr);
   w_forr=flipud(w_forr);
end

%%%%%%%%%%%%%%% Adjust strengths to the correct temperature %%%%%%%%%%%%
% NOTE that the stimulated emiision terms ARE included.
[a1,b1,c1,d1]=qqttiippss(isotope);
Qt=a1 + b1*temperature + c1*temperature^2 + d1*temperature^3;
Q296=a1 + b1*296 + c1*296^2 + d1*296^3;

numer=Q296*exp(stuff.btz*elowerr/296).*(1-exp(-stuff.btz*freqr/temperature));
denom=Qt*exp(stuff.btz*elowerr/temperature).*(1-exp(-stuff.btz*freqr/296));
strenrt=strenr.*numer./denom;

[jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr]=... 
    orderer(jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr); 
stuff.population=zeros(size(jr));
stuff.population_t=zeros(size(jr));
