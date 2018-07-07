function [jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff]=... 
    loader662(temperature,band,path_length,ptotal,pself,ii); 

%%%%%%%%%%%%%%%%%%%%  loader.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file will:
%	(1) Load in constants and other parameters
%	(2) Load the HITRAN *.mat data file 
%	(3) Select out the lines with the following qualifications:
%		a) Is a Q-branch line.
%		b) Has the correct user-specified upper and lower
%		   vibrational states.
%		c) Has the main isotope. (iso=1)
%		d) Has J<= 50
%	(4) Correct widths to desired temprature.
%
%    The main idea behind loading the data in from the original hittomat
%  file every time the program is run is to make sure that the results are
%  reproducible and that anyone else who would like to run the programs
%  will know exactly how this original data is obtained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stuff.evenodd=ii;
stuff.p1=2.0649e-02; 
stuff.p2=1.5840e-01; 

stuff.path_length=path_length;       %path length in cm 
stuff.pressure_self=pself;           %already in atm  
stuff.pressure_for=ptotal-pself;     %already in atm  

%beta=input('Enter beta : ');  %0.63; 
%  Beta is proportional to the probability of a collision causing a transition
%  to  go from (f to f) as compared to (f to e).  i.e. to switch symmetry.
% from /salsify/scratch4/Strow/Tobin_home/tobin/Co2q/B_sigpie/Run/15um/run2.m
% beta_self=0.5358;            %self broadened        %%%%%used to be 0.5
% beta_for=0.6002;             %foreign broadened     %%%%%used to be 0.62

%%%%%% use papaya,mango for dummy beta params
[duration_pure,duration_for,beta_self,beta_for,papaya,mango]=...
       co2_param(band,ptotal,pself);
beta_pi=(pself*beta_self+(ptotal-pself)*beta_for)/ptotal;
stuff.beta=beta_pi;

stuff.bsm=1; 
stuff.btz=1.4387863; 
stuff.B0=0.4; 
% temperature_ref must go along with the value of density_ref (at 273.15 k) 
stuff.temperature_ref=273.15; 
stuff.speed_light=2.99792458e8; 
stuff.pressure_ref=1; 
stuff.Boltzmann=1.380662e-23; 
stuff.density=2.6867e19;
stuff.density=101325/stuff.Boltzmann/stuff.temperature_ref/1e6;
stuff.mass_proton=1.6726485e-27; 
if ((band ~= 648) & (band ~= 662))
  stuff.mass_CO2=44*stuff.mass_proton; 
elseif (band == 648)
  stuff.mass_CO2=45*stuff.mass_proton; 
elseif (band == 662)
  stuff.mass_CO2=46*stuff.mass_proton; 
  end
stuff.band=band;

stuff.duration=0.0000000000000; %just a dummy for these Q branches
duration=(pself*duration_pure+(ptotal-pself)*duration_for)/ptotal;
stuff.duration=duration;

stuff.frequency_shift=0.000000; %just a dummy for these Q branches

%%%%%%%%%%%%%%%%%%%%% Load in hittomat-produced file %%%%%%%%%%%%%%%%%%%%%%%%
if (band == 618)
  load ../CO2_MATFILES/hit618
elseif (band == 648)
  load ../CO2_MATFILES/hit648
elseif (band == 662)
  load ../CO2_MATFILES/hit662
elseif (band == 667)
  load ../CO2_MATFILES/hit667
elseif (band == 720)
  load ../CO2_MATFILES/hit720
elseif (band == 791)
  load ../CO2_MATFILES/hit791
elseif (band == 2080)
  load ../CO2_MATFILES/hit2080
else
  error('wrong band in loading the hitBLAH file!!!')
  end

%%%%%%%%%%%%%%%%%%%%% Sort out lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter in lower and upper vibrational states:
% (Note that the user should enter the INDEX of the vib. state, not the actual
% state.  (for example, enter "2", not "01101")  Details on obtaining these
% indices from the known vibrational level and vice-versa can be found in the
% file Hit_co2q/README_viblevels.)

isotope=1;
if band==618
        v_l=2;v_u=3;
end
if band==648
        v_l=1;v_u=2; isotope=2;
end
if band==662
        v_l=1;v_u=2; isotope=3;
end
if band==667
        v_l=1;v_u=2;
end
if band==720
        v_l=2;v_u=5;
end
if band==791
        v_l=3;v_u=8;
end        
if band==2080
	v_l=1;v_u=8;
end

junk_str=j_lower(:,6:8);
for n=1:length(junk_str);
        eval(['tmp(n)=' junk_str(n,:) ';']);
        end
tmp=tmp';

stuff.isotope=isotope;

% Find indices of lines which meet the above requirements:
index=find(j_lower(:,5)=='Q' & v_lower==v_l & v_upper==v_u & iso==isotope);

% Now select out these lines from the rest.
j_lowerq(:,1:9)=j_lower(index,1:9);		% J of lower state (strings)
%jq=num(j_lowerq(:,6:8)); 			% J of lower state (integer)

junk_str=j_lowerq(:,6:8);
for n=1:length(junk_str);
        eval(['jq(n)=' junk_str(n,:) ';']);
        end
jq=jq';

freqq=freq(index);				% frequencies (cm-1)
strenq=stren(index);				% line strengths (296K)
elowerq=elower(index);				% lower state TOTAL energy
w_forq=w(index);				% Air broadened widths (296K)
w_fortempq=w_temp(index);			% Temperature correction 
						% coefficients for w_forq.
w_selfq=w_s(index);				% Self-broadened widths (296K)
w_selftempq=0.685;				% coefficient for self-widths

%%%%%%%%%%now do ODD vs EVEN
if (mod(ii,2) == 0) %even j's
  index=find(mod(jq,2) == 0);
else
  index=find(mod(jq,2) == 1);
  end
jq=jq(index);
freqq=freqq(index);
strenq=strenq(index);
elowerq=elowerq(index);
w_forq=w_forq(index);
w_fortempq=w_fortempq(index);
w_selfq=w_selfq(index);

%%%%%%%% Flip arrays (if needed) so that lowest J's appear first  %%%%%%%%%% 
if jq(1)>jq(2) 
   j_lowerq=flipud(j_lowerq);freqq=flipud(freqq);strenq=flipud(strenq);
   elowerq=flipud(elowerq);w_selfq=flipud(w_selfq);jq=flipud(jq);
   w_forq=flipud(w_forq);
end

%%%%%%%%%%%%%%%%%% Correct widths to desired temperature %%%%%%%%%%%%%%
% This next line corrects the air broadened widths at 296K, wq, to widths
% at the correct temperature.  It is done by using the equation:
%
%	wq= wq*(296/T)^w_tempq         from AGLG paper
w_forq=w_forq.*(296/temperature).^w_fortempq;

w_selfq=w_selfq.*(296/temperature)^w_selftempq;

%%%%%%%%%%%%%%% Adjust strengths to the correct temperature %%%%%%%%%%%%
% NOTE that the stimulated emiision terms ARE included.

% These are the partition function constants for CO2 obtained from qtips.pas
% for the third isotope
[a1,b1,c1,d1]=qqttiippss(isotope);

% Partition functions evaluated at desired temperature and 296K.
Qt=a1 + b1*temperature + c1*temperature^2 + d1*temperature^3;
Q296=a1 + b1*296 + c1*296^2 + d1*296^3;

% See AFGL HITRAN database paper for details
numer=Q296*exp(stuff.btz*elowerq/296).*(1-exp(-stuff.btz*freqq/temperature));
denom=Qt*exp(stuff.btz*elowerq/temperature).*(1-exp(-stuff.btz*freqq/296));
strenqt=strenq.*numer./denom;

[jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq]=...  
    orderer(jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq);  

stuff.population=zeros(size(jq)); 
stuff.population_t=zeros(size(jq)); 

clear accuracy dipole elower freq gas_id iso j_lower j_upper line_status
clear p_shift reference stren v_lower v_upper w w_s w_temp j_upper index
clear hitfile v_l v_u mass_proton









