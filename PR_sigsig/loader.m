function [jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr,stuff]=... 
    loader(temperature,band,path_length,ptotal,pself,prb); 

format long e

global p2311_21_jmax p2350_jmax r2350_jmax pr2351_jmax

stuff.band = band;

stuff.p1 = 2.0649e-02; 
stuff.p2 = 1.5840e-01; 

stuff.path_length   = path_length;       %path length in cm
stuff.pressure_self = pself;           %already in atm 
stuff.pressure_for  = ptotal-pself;     %already in atm 

if (prb == 'R')
  if (band == 2350)
    v_l = 1;
    v_u = 9;
    isotope = 1;
  elseif (band == 2351)
    v_l = 1;
    v_u = 9;
    isotope = 2;
  elseif (band == 2352)
    v_l = 1;
    v_u = 9;
    isotope = 3;
  elseif (band == 2353)
    v_l = 3;
    v_u = 23;
    isotope = 1;
  elseif (band == 2354)
    v_l = 5; 
    v_u = 25;
    isotope=1;
    end
  stuff.prb='R';
elseif (prb == 'P')
  if (band == 2350)
    v_l = 1;
    v_u = 9;
    isotope = 1;
  elseif (band == 2351)
    v_l = 1;
    v_u = 9;
    isotope = 2;
  elseif (band == 2352)
    v_l = 1;
    v_u = 9;
    isotope = 3;
  elseif (band == 2353)
    v_l = 3;
    v_u = 23;
    isotope = 1;
  elseif (band == 2354)
    v_l = 5;
    v_u = 25;
    isotope = 1;
    end
  stuff.prb='P';
else
  error('need prb ==== P or R')
  end

stuff.isotope = isotope;

%%%%%%%%%%%%%%%%%%%%%%%%%% do line mixing parameters %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% use papaya,mango for dummy beta params
[duration_pure,duration_for,beta_pure,beta_for,papaya,mango]=...
        co2_param(band,ptotal,pself);

frequency_shift = 0;

stuff.beta_pure = beta_pure;
stuff.beta_for  = beta_for;
beta            = (pself*beta_pure+(ptotal-pself)*beta_for)/ptotal;
if (beta > 1)
  beta = 1.0;
  end
duration        = (pself*duration_pure+(ptotal-pself)*duration_for)/ptotal;

bsm = 1;   band_strength_multiplier = bsm; 

stuff.frequency_shift = frequency_shift;
stuff.beta            = beta;
stuff.duration        = duration;

stuff.bsm = 1;
stuff.btz = 1.4387863;
stuff.B0  = 0.4;
% temperature_ref must go along with the value of density_ref (at 273.15 k)
stuff.temperature_ref = 273.15;
stuff.speed_light     = 2.99792458e8;
stuff.pressure_ref    = 1;

stuff.density = 2.6867e19;     %%%%loschmidt number

%%%p/kT == 101300/1.38e-23/273 = 2.6889e25 m-3 = 2.6889e19 cm-3
%%%%%%stuff.density = 2.6867775e19;

stuff.Boltzmann   = 1.380662e-23;
stuff.mass_proton = 1.6726485e-27;
stuff.mass_CO2    = 44*stuff.mass_proton;

if (band == 2350) 
  load ../CO2_MATFILES/hit2350
elseif (band == 2351) 
  load ../CO2_MATFILES/hit2351
elseif (band == 2352) 
  load ../CO2_MATFILES/hit2352
elseif (band == 2353) 
  load ../CO2_MATFILES/hit2353
elseif (band == 2354) 
  load ../CO2_MATFILES/hit2354
else 
  error('wrong band in loading the hitBLAH file!!!') 
  end 

%%%%%%%%%%%%%%%%%%%%% Sort out lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
junk_str = j_lower(:,6:8);
for n=1:length(junk_str);
  eval(['tmp(n)=' junk_str(n,:) ';']);
        end
  tmp=tmp';

ttt = 2000;       %max j
%ttt=input('enter ttt ');
if ((stuff.prb == 'P') & (length(intersect(band,[2311,2321])) > 0))
  ttt=p2311_21_jmax;
elseif ((stuff.prb == 'P') & (length(intersect(band,[2351])) > 0))
  ttt=pr2351_jmax;
elseif ((stuff.prb == 'R') & (length(intersect(band,[2351])) > 0))
  ttt=pr2351_jmax;
elseif ((stuff.prb == 'P') & (length(intersect(band,[2350,2320,2310])) > 0))
  ttt=p2350_jmax;
elseif ((stuff.prb == 'R') & (length(intersect(band,[2350])) > 0))
  ttt=r2350_jmax;
  end

%index=find(j_lower(:,5)==prb & v_lower==v_l & v_upper==v_u & iso==isotope);
index=find(j_lower(:,5)==prb & v_lower==v_l & v_upper==v_u & iso==isotope ...
            & tmp <= ttt);

% Now select out these lines from the rest.
j_lowerr(:,1:9) = j_lower(index,1:9);		% J of lower state (strings)
junk_str        = j_lowerr(:,6:8);
for n = 1:length(junk_str);
	eval(['jr(n)=' junk_str(n,:) ';']);
  end
jr = jr';
freqr   = freq(index);				% frequencies (cm-1)
%%% for some reson, this was commented out when I looked on March 3, 2004
%%% leave it COMMENTED OUT, as all other bands have df = 0.0
%%% freqr   = freq(index) + ptotal*p_shift(index);  %freq shift
strenr  = stren(index);				% line strengths
elowerr = elower(index);			% lower state total energy
w_forr  = w(index);				% Air broadened widths (296K)
w_fortempr = w_temp(index);			% Temperature correction
						%      coefficients to air-
						%         broadened widths
w_selfr = w_s(index);				% Self-broadened widths (296K)
% w_selftempq=0.685				% coefficients for self-widths

%%%%%%%%%%%%%%%%%% Correct widths to desired temperature %%%%%%%%%%%%%%
w_forr  = w_forr.*(296/temperature).^w_fortempr;
w_selfr = w_selfr.*(296/temperature)^0.685;

%%%%%%%% Flip arrays (if needed) so that lowest J's appear first  %%%%%%%%%% 
if jr(1)>jr(2) 
   j_lowerr = flipud(j_lowerr);
   freqr    = flipud(freqr);
   strenr   = flipud(strenr);
   elowerr  = flipud(elowerr);
   w_selfr  = flipud(w_selfr);
   jr       = flipud(jr);
   w_forr   = flipud(w_forr);
end

%%%%%%%%%%%%%%% Adjust strengths to the correct temperature %%%%%%%%%%%%
% NOTE that the stimulated emiision terms ARE included.
[a1,b1,c1,d1] = qqttiippss(isotope);       	
Qt   = a1 + b1*temperature + c1*temperature^2 + d1*temperature^3;
Q296 = a1 + b1*296         + c1*296^2         + d1*296^3;


%%%%%%%this is just to test against qtips
bbb = -1;
if bbb > 0
  [A,B,C,D,G] = qtips(2,[1 2 3 4 5 6 7 8]);
  A11 = A(isotope); B11 = B(isotope); C11 = C(isotope); D11 = D(isotope);
  [A11 B11 C11 D11];
  qfcn = Q296./Qt;
  E_li = elowerr; 
  T    = temperature; 
  v0   = freqr; 
  s00  = strenr; %(don't worry about amt)
  %s00=s0*6.022045e26;                        %or could do amt=amt*6.022e26
  c2   = 1.4387863;                           %K/ cm-1  from Genln2 manual
  sb   =exp(-c2*E_li/T)./exp(-c2*E_li/296.0);  %boltzman factor at T
  se   = (1-exp(-c2*v0/T))./(1-exp(-c2*v0/296.0)); %adjust for detailed balance
  %strength = amt*(qfcn'.*s00.*sb.*se);
  strength  = (qfcn'.*s00.*sb.*se);
  end

%%%%%%!!!!!!!!!! the orig code!!!!!!!!!!!!!
numer = Q296*exp(stuff.btz*elowerr/296).*(1-exp(-stuff.btz*freqr/temperature));
denom = Qt*exp(stuff.btz*elowerr/temperature).*(1-exp(-stuff.btz*freqr/296));
strenrt = strenr.*numer./denom;
%%%%%%%%%%

[jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr]=... 
    orderer(jr,elowerr,w_forr,w_selfr,freqr,strenrt,strenr); 
stuff.population   = zeros(size(jr));
stuff.population_t = zeros(size(jr));

