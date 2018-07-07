function [stuff] = ...
      loader_jjohns(temperature,band,path_length,ptotal,pself,prb,region); 

%%%%%%%%%%%%%%% same as loader.m in PR_sigsig directory, except that it
%%%%%%%%%%%%%%% loads in the co2_params_JJOHNS files, depending on REGION

format long e

global p2311_21_51_jmax

stuff.band = band;

stuff.p1=2.0649e-02; 
stuff.p2=1.5840e-01; 

stuff.path_length=path_length;       %path length in cm
stuff.pressure_self=pself;           %already in atm 
stuff.pressure_for=ptotal-pself;     %already in atm 

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
    isotope = 1;
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

%band
%prb
stuff.isotope = isotope;

%%%%%%%%%%%%%%%%%%%%%%%%%% do line mixing parameters %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% use papaya,mango for dummy beta params
%%%%%% call co2_param_JJOHNS2002(band,ptotal,pself);
[duration_pure,duration_for,beta_pure,beta_for,papaya,mango]=...
        co2_param_JJOHNS2002(band,ptotal,pself,region);

frequency_shift = 0;

stuff.beta_pure = beta_pure;
stuff.beta_for  = beta_for;
beta = (pself*beta_pure+(ptotal-pself)*beta_for)/ptotal;
if (beta > 1)
  beta=1.0;
  end
duration = (pself*duration_pure+(ptotal-pself)*duration_for)/ptotal;

bsm = 1;  band_strength_multiplier = bsm; 

stuff.frequency_shift = frequency_shift;
stuff.beta            = beta;
stuff.duration        = duration;

