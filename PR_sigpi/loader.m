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

if (band == 618)
  if (prb == 'R')
    v_l=2;v_u=3;isotope=1;
    stuff.prb='R';
  elseif (prb == 'P')
    v_l=2;v_u=3;isotope=1;
    stuff.prb='P';
  elseif (prb == 'Q')
    v_l=2;v_u=3;isotope=1;
    stuff.prb='Q';
%    error('need prb ==== P or R')
    end
  end

if (band == 648)
  if (prb == 'R')
    v_l=1;v_u=2;isotope=2;
    stuff.prb='R';
  elseif (prb == 'P')
    v_l=1;v_u=2;isotope=2;
    stuff.prb='P';
  elseif (prb == 'Q')
    v_l=1;v_u=2;isotope=2;
    stuff.prb='Q';
%    error('need prb ==== P or R')
    end
  end

if (band == 662)
  if (prb == 'R')
    v_l=1;v_u=2;isotope=3;
    stuff.prb='R';
  elseif (prb == 'P')
    v_l=1;v_u=2;isotope=3;
    stuff.prb='P';
  elseif (prb == 'Q')
    v_l=1;v_u=2;isotope=3;
    stuff.prb='Q';
%    error('need prb ==== P or R')
    end
  end

if (band == 667)
  if (prb == 'R')
    v_l=1;v_u=2;isotope=1;
    stuff.prb='R';
  elseif (prb == 'P')
    v_l=1;v_u=2;isotope=1;
    stuff.prb='P';
  elseif (prb == 'Q')
    v_l=1;v_u=2;isotope=1;
    stuff.prb='Q';
%    error('need prb ==== P or R')
    end
  end

if (band == 720)
  if (prb == 'R')
    v_l=2;v_u=5;isotope=1;
    stuff.prb='R';
  elseif (prb == 'P')
    v_l=2;v_u=5;isotope=1;
    stuff.prb='P';
  elseif (prb == 'Q')
    v_l=2;v_u=5;isotope=1;
    stuff.prb='Q';
%    error('need prb ==== P or R')
    end
  end

if (band == 791)
  if (prb == 'R')
    v_l=3;v_u=8;isotope=1;
    stuff.prb='R';
  elseif (prb == 'P')
    v_l=3;v_u=8;isotope=1;
    stuff.prb='P';
  elseif (prb == 'Q')
    v_l=3;v_u=8;isotope=1;
    stuff.prb='Q';
%    error('need prb ==== P or R')
    end
  end

if (band == 1932)
  if (prb == 'R')
    v_l=1;v_u=6;isotope=1;
    stuff.prb='R';
  elseif (prb == 'P')
    v_l=1;v_u=6;isotope=1;
    stuff.prb='P';
  elseif (prb == 'Q')
    v_l=1;v_u=6;isotope=1;
    stuff.prb='Q';
%    error('need prb ==== P or R')
    end
  end

if (band == 2080)
  if (prb == 'R')
    v_l=1;v_u=8;isotope=1;
    stuff.prb='R';
  elseif (prb == 'P')
    v_l=1;v_u=8;isotope=1;
    stuff.prb='P';
  elseif (prb == 'Q')
    v_l=1;v_u=8;isotope=1;
    stuff.prb='Q';
%    error('need prb ==== P or R')
    end
  end

if (band == 2129)
  if (prb == 'R')
    v_l=2;v_u=15;isotope=1;
    stuff.prb='R';
  elseif (prb == 'P')
    v_l=2;v_u=15;isotope=1;
    stuff.prb='P';
  elseif (prb == 'Q')
    v_l=2;v_u=15;isotope=1;
    stuff.prb='Q';
%    error('need prb ==== P or R')
    end
  end


stuff.isotope=isotope;

%%%%%% use papaya,mango for dummy beta params
[duration_self,duration_for,beta_self,beta_for,papaya,mango]=...
     co2_param(band,ptotal,pself);

beta_pi  = (pself*beta_self+(ptotal-pself)*beta_for)/ptotal; 
duration = (pself*duration_self+(ptotal-pself)*duration_for)/ptotal; 

stuff.beta            = beta_pi; 
stuff.duration        = duration;
stuff.frequency_shift = 0.0;

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
if ((band ~= 648) & (band ~= 662)) 
  stuff.mass_CO2=44*stuff.mass_proton;  
elseif (band == 648) 
  stuff.mass_CO2=45*stuff.mass_proton;  
elseif (band == 662) 
  stuff.mass_CO2=46*stuff.mass_proton;  
  end 
stuff.band = band; 

%%%%%%%%this is NEW and needed for the PR branches
stuff.band = band;

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
if (band == 618)
  %load hit618
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
elseif (band == 1932)
  load ../CO2_MATFILES/hit1932
elseif (band == 2080)
  load ../CO2_MATFILES/hit2080
elseif (band == 2129)
  load ../CO2_MATFILES/hit2129
else
  error('wrong band in loading the hitBLAH file!!!')
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
index=find(j_lower(:,5)==prb & v_lower==v_l & v_upper==v_u & iso==isotope);

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
% These are the partition function constants for CO2 obtained from qtips.pas
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










