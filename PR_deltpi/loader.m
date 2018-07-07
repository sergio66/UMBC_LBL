function [jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq,stuff]=... 
    loader(temperature,band,path_length,ptotal,pself,prb); 

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

stuff.p1=2.0649e-02; 
stuff.p2=1.5840e-01; 
 
stuff.path_length=path_length;       %path length in cm 
stuff.pressure_self=pself;           %already in atm  
stuff.pressure_for=ptotal-pself;     %already in atm  

stuff.prb=prb;

if (band == 668)
  if (prb == 'R')
    v_l=2;v_u=4;isotope=1;
    stuff.prb='R';
  elseif (prb == 'P')
    v_l=2;v_u=4;isotope=1;
    stuff.prb='P';
  elseif (prb == 'Q')
    v_l=2;v_u=4;isotope=1;
    stuff.prb='Q';
%    error('need prb ==== P or R')
    end
  end

if (band == 740)
  if (prb == 'R')
    v_l=4;v_u=8;isotope=1;
    stuff.prb='R';
  elseif (prb == 'P')
    v_l=4;v_u=8;isotope=1;
    stuff.prb='P';
  elseif (prb == 'Q')
    v_l=4;v_u=8;isotope=1;
    stuff.prb='Q';
%    error('need prb ==== P or R')
    end
  end

stuff.isotope=isotope;



%%%%%%%%%%%%%%%%%%%%%%%%% THIS WAS ORIGINAL CODE %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% just use one dummy beta here!!!!!!!!!
%%%%%%%%%%%% this beta is NOT used in the overall PR_deltpi coding!!!!!!!
%%%%%%%%%%%% only important parameter is duration_self,duration_for
%use this d.of.c parameter everywhere!!!!!!!
%duration_self = 9.244381131571670e-03;
%duration_for  = 3.872568524276932e-03;

%[duration_self,duration_for,beta_self,beta_for,papaya,mango]=co2_param(band);
%beta_pi=(pself*beta_self+(ptotal-pself)*beta_for)/ptotal;  
%duration=(pself*duration_self+(ptotal-pself)*duration_for)/ptotal;  
%stuff.beta=beta_pi;  
%%%%%%%%%%%%%%%%%%%%%%%%% THIS WAS ORIGINAL CODE %%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% THIS IS NEW CODE ( FROM Q_DELTPI) %%%%%%%%%%%%%%%%%%%%%%%%%%
%beta_self=0.5488; %input('Enter beta_self : ');
%beta_for=0.6762;  %input('Enter beta_for : ');
%stuff.beta_pi_self   = 0.5358;  % fitted from pure Pi-Sigma data
%stuff.beta_pi_air    = 0.6002;  % fitted from air-broadened Pi-Sigma data
%stuff.beta_delt_self = 0.5709;  % fitted from pure Pi-Delta data
%stuff.beta_delt_air  = 0.7776;  % fitted from air-broadened Pi-Delta data
%stuff.beta_pi_self   = 0.4;  % fitted from pure Pi-Sigma data
%stuff.beta_pi_air    = 0.4;  % fitted from air-broadened Pi-Sigma data
%stuff.beta_delt_self = 0.4;  % fitted from pure Pi-Delta data
%stuff.beta_delt_air  = 0.4;  % fitted from air-broadened Pi-Delta data

%/sals/scratch4/Strow/Tobin_home/tobin/Co2q/B_deltpie/Mfiles/Fit/nfive400.mat
%this corresponds to file 3 in loader-mod in  
%/salsify/data/Strow/Tobin/Co2q_B_deltpie
beta_pi_self   = 0.5407505;  % fitted from pure Pi-Sigma data
beta_pi_air    = 0.5995593;  % fitted from air-broadened Pi-Sigma data
beta_delt_self = 0.4880067;  % fitted from pure Pi-Delta data
beta_delt_air  = 0.72172644;  % fitted from air-broadened Pi-Delta data

%%%%%% duration of collision VERY important here, so instead of the following
%%%%%%which is in Q_deltpi 
%%%% [papaya,mango,beta_pi_self,beta_pi_air,beta_delt_self,beta_delt_air] = ...
%%%%               co2_param(band);
%%%% ie ignore use papaya,mango and use real variables instead

[duration_self,duration_for,beta_pi_self,beta_pi_air,...
       beta_delt_self,beta_delt_air] = co2_param(band,ptotal,pself);

duration=(pself*duration_self+(ptotal-pself)*duration_for)/ptotal;  

stuff.beta_pi_self   = beta_pi_self;
stuff.beta_pi_air    = beta_pi_air;
stuff.beta_delt_self = beta_delt_self;
stuff.beta_delt_air  = beta_delt_air;
%%%%%%%%%% THIS IS NEW CODE ( FROM Q_DELTPI) %%%%%%%%%%%%%%%%%%%%%%%%%%

stuff.duration        = duration;
stuff.frequency_shift = 0.0;

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

%%%%%%%%this is NEW and needed for the PR branches
stuff.band=band;

%%%%%%%%%%%%%%%%%%%%% Load in hittomat-produced file %%%%%%%%%%%%%%%%%%%%%%%%
% Details on how these files were created are found in subdirectory /Hit_co2q
% These files contain the HITRAN line parameters (1990 version).
%disp('bands:  618, 667, 720, 791 2080')
%band=input('What band do you want? ');
%val(['load /salsify/users/tobin/Co2q/Hit_co2q/co2_' num2str(band) '_new.mat'
%string='/carrot/users/tobin/Hittomat/Q_params_92/';
%eval(['load ' string 'co2_' num2str(band) '_92.mat'])

%estring='load /carrot/users/tobin/Hittomat/Q_params_92/co2_'; 

%estring='load /beet/users/tobin/Hittomat/Q_params_92/co2_'; 
%eval([estring num2str(band) '_92']) 
%load /beet/users/tobin/Hittomat/Src/hit_co2_15um
if (band == 740)
  load ../CO2_MATFILES/hit740
elseif (band == 668)
  load ../CO2_MATFILES/hit668
else
  error('wrong band in loading the hitBLAH file!!!')
  end
 
%%%%%%%%%%%%%%%%%%%%% Sort out lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter in lower and upper vibrational states:
% (Note that user should enter the INDEX of the vib. state, not the actual
% state.  (for example, enter "2", not "01101")  Details on obtaining these
% indices from the known vibrational level and vice-versa can be found in the
% file Hit_co2q/README_viblevels.)

if band==668
        v_l=2;v_u=4;
        end
if band==2093
        v_l=2;v_u=14;
        end
if band==740
        v_l=4;v_u=8;
        end


%%%%%%%%%%%%%%%%%%%%% Sort out lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk_str=j_lower(:,6:8);
for n=1:length(junk_str);
  eval(['tmp(n)=' junk_str(n,:) ';']);
  end
tmp=tmp';

% Find indices of lines which meet the above requirements:
%%index=find(j_lower(:,5)==prb & v_lower==v_l & v_upper==v_u & iso==1 & ...
%%           tmp<=100);
index=find(j_lower(:,5)==prb & v_lower==v_l & v_upper==v_u & iso==1);

%%%Q branch code has tmp <= 50; 
%%%his code on /beet/users/tobin/Co2q/B_sigpie/Run/15um/PR720 has tmp <= 100

% Now select out these lines from the rest.
j_lowerq(:,1:9)=j_lower(index,1:9);		% J of lower state (strings)
junk_str=j_lowerq(:,6:8);
for n=1:length(junk_str);
  eval(['jq(n)=' junk_str(n,:) ';']);
  end
jq=jq';

freqq=freq(index);			      % frequencies (cm-1)
strenq=stren(index);			      % line strengths (296K)
elowerq=elower(index);			      % lower state TOTAL energy
w_forq=w(index);			      % Air broadened widths (296K)
w_fortempq=w_temp(index);		      % Temperature correction 
					      % coefficients for w_forq.
w_selfq=w_s(index);			      % Self-broadened widths (296K)
w_selftempq=0.685;			      % coefficient for self-widths

%%%%%%%% Flip arrays (if needed) so that lowest J's appear first  %%%%%%%%%% 
if jq(1)>jq(2) 
   j_lowerq=flipud(j_lowerq);freqq=flipud(freqq);strenq=flipud(strenq);
   elowerq=flipud(elowerq);w_selfq=flipud(w_selfq);jq=flipud(jq);
   w_forq=flipud(w_forq);
end

%%%%%%%%%%%%%%%%%% Correct widths to desired temperature %%%%%%%%%%%%%%
%
% This next line corrects the air broadened widths at 296K, wq, to widths
% at the correct temperature.  It is done by using the equation:
%
%	wq= wq*(296/T)^w_tempq         from AGLG paper
w_forq=w_forq.*(296/temperature).^w_fortempq;

% Now I will adjust the self-broadened widths to the correct temperature
% using the exponent n=0.685 which was obtained by Liu in her thesis
% as the best value to use for self-broadened CO2 widths.
w_selfq=w_selfq.*(296/temperature)^w_selftempq;

%%%%%%%%%%%%%%% Adjust strengths to the correct temperature %%%%%%%%%%%%
% NOTE that the stimulated emiision terms ARE included.
[a1,b1,c1,d1]=qqttiippss(isotope);

% Partition functions evaluated at desired temperature and 296K.
Qt=a1 + b1*temperature + c1*temperature^2 + d1*temperature^3;
Q296=a1 + b1*296 + c1*296^2 + d1*296^3;

% See AFGL HITRAN database paper for details
numer=Q296*exp(stuff.btz*elowerq/296).*(1-exp(-stuff.btz*freqq/temperature));
denom=Qt*exp(stuff.btz*elowerq/temperature).*(1-exp(-stuff.btz*freqq/296));
strenqt=strenq.*numer./denom;

%[jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq]=...  
%    orderer(jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq);  

stuff.population=zeros(size(jq)); 
stuff.population_t=zeros(size(jq)); 

clear accuracy dipole elower freq gas_id iso j_lower j_upper line_status
clear p_shift reference stren v_lower v_upper w w_s w_temp j_upper index
clear hitfile v_l v_u mass_proton










