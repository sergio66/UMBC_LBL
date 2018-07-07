%%this is for the RAL data, which is good from 2230- cm-1

%%%%%%%%%%%%%%%%%%%  loader.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('*******************************************************************');
disp('i hope that you did the following (for datafits run7co2_JJOHNS) ...');
disp('    (a) set dofudge = +1 in PR_sigsig/driver4um');
disp('    (b) set dofudge = +1 in PR_sigsig/yrun_sigsig');
disp('    (c) set co2_param.m -> co2_param_AIRS2002.m');
disp('    (d) set iWeight = -1 in PR_sigsig/driver4um');
disp('    (e) set dofudge = +1 in CO2_COMMON/co2_param ');
disp('*******************************************************************');

format long e

disp('------------------------------------------------------------------')
disp(' #  File Name  T(K)     P(CO2)  P(tot) P.L.(m)  REGION   resolution')
disp(' co2a :                  HPa     HPa                       cm-1');
disp('------------------------------------------------------------------')
disp(' 1  airs 01    296.0    26.80    26.80  512.746    4um     0.0038')
disp(' 2  airs 41    296.0    26.80   161.40  512.746    4um     0.0038')
disp(' 3  airs 71    296.0    26.80   560.30  512.746    4um     0.0038')
disp(' 4  airs 91    296.0    26.80   961.70  512.746    4um     0.0038')
disp('-------------------------------------------------------------------')

disp('------------------------------------------------------------------')
disp(' #  File Name  T(K)     P(CO2)  P(tot) P.L.(m)  REGION   resolution')
disp(' co2a :                  HPa     HPa                       cm-1');
disp('------------------------------------------------------------------')
disp(' 1  airs 01    291.0    26.80   1.07    512.746    4um     0.0038')
disp(' 2  airs 11    276.0    26.80   0.8040  512.746    4um     0.0038')
disp(' 3  airs 21    259.0    26.80   0.5780  512.746    4um     0.0038')
disp(' 4  airs 31    241.0    26.80   0.3940  512.746    4um     0.0038')
disp(' 5  airs 41    221.0    26.80   0.2500  512.746    4um     0.0038')
disp(' 6  airs 51    216.0    26.80   0.1449  512.746    4um     0.0038')
disp(' 7  airs 61    216.0    26.80   0.0734  512.746    4um     0.0038')
disp('-------------------------------------------------------------------')

disp('Pick which data files to include in the fit.  Enter 0 to quit')
stop=0;
i=0;
while stop==0
  junk=input('Enter raw data file #: ');
  if junk ~= 0
    i=i+1;
    pick_no(i)=junk;
  else
    stop=1;
    end
  end

load /asl/data/ral/Co2a/COUSIN2351/cousin2351.mat

K_back_all = zeros(size(kcousin2351));
nitrogen_all = zeros(size(kcousin2351));

fname = '/home/sergio/SPECTRA/CO2_COUSIN_2351/IPFILES/some_co2'; 
fudge = ones(1,7);

ipfile = load(fname);
ps = ipfile(:,3);
pf = ipfile(:,2);   %%actually, this is ptotal
pf = pf - ps;
temp = ipfile(:,4);
fitted_t = temp;

%%% need path length in cm
pl = ipfile(:,5);
MGC=8.314674269981136;
pl = pl*1e9*MGC/101325 .* temp ./ps;

for i=1:length(pick_no);
  ni = num2str(pick_no(i));
  NI = ni;
  ipfileIN = ['/home/sergio/SPECTRA/CO2_COUSIN_2351/IPFILES/some_co2_set' ni]
  temperature(i)  = fitted_t(pick_no(i));
  path_length(i)  = pl(pick_no(i));
  pressure_self(i)= ps(pick_no(i));
  pressure_for(i) = pf(pick_no(i));
  end

clear fil filb filn2 filw fitted_t pl ps pf

% Load in raw data and background calcs
maxnpts=0;
for i=1:length(pick_no)
  fprintf(1,'loaded in synthetic cousin data ... haha fibber ... \n');
  ff = fr';

  ni=num2str(pick_no(i))
 
  pah = ['rawdata = kcousin2351(' ni ',:);'];
  eval(pah)
  rawdata = exp(-rawdata);

  ind=find((ff >= 2255) & (ff <= 2305));
  ind=find((ff >= 2270) & (ff <= 2290));
  ff=ff(ind);
  rawdata=rawdata(ind)';

  chosen = 1 : 10 : length(ff);

  fdg=fudge(pick_no(i));

  ni=num2str(i)

  %%%this is the data, at the lab spacing; 
  eval(['f' ni '= ff(chosen);']);
  eval(['t_rawdata' ni '= rawdata(chosen)*fdg;']);
  eval(['maxnpts=max([ maxnpts length(ff(chosen))]);']);
  disp(['f, t_rawdata loaded in.']);

  %%%this is the KCARTA runs at 0.0025 spacing
  eval(['nitrogen0_' ni ' = nitrogen_all(' num2str(pick_no(i)) ',:);']);
  eval(['weak_run7co2' ni '  = K_back_all (' num2str(pick_no(i)) ',:);']);
  eval(['clear lor lor_mod '])
  eval(['length(f' ni ')'])
  eval(['gah0 = weak_run7co2' ni ';']);
  disp(['f, t_rawdata, weak_run7co2: ' ni ' loaded in.'])
  
  clf
  plot(ff(chosen),rawdata(chosen),fr,exp(-gah0),'r');
  title('co2 background lines (r), co2 data (b)')
  pause(0.1)

  %%%now we need to interpolate the weak background lines, plus N2
  %%%onto the lab data spacing
  eval(['nitrogen' ni  ' = interp1(fr,nitrogen0_' ni ',f' ni ');']);
  eval(['k_co2weak' ni ' = interp1(fr,gah0,f' ni ');']);

  red_num=input('Reduce # data points by factor of : ','s');
  pah1 = ['[f' ni ',t_rawdata' ni ',k_co2weak' ni ',nitrogen' ni '] ='];
  pah2 = ['reduce4(f' ni ',t_rawdata' ni ',k_co2weak' ni ];
  pah3 = [', nitrogen' ni ',' red_num ');'];
  pah = [pah1 pah2 pah3];
  eval(pah);

  eval(['maxnpts=min([ maxnpts length(f' ni ')]);']);
  eval(['gah = k_co2weak' ni ';']);
  eval(['plot(f' ni ',t_rawdata' ni ',f' ni ', exp(-gah),''r'')'])
  title('after reduction co2 background lines (r), co2 data (b)')
  pause(0.1);

  end

clear ff

f          = zeros(maxnpts,length(pick_no));
t_rawdata  = zeros(maxnpts,length(pick_no));
k_co2weak  = zeros(maxnpts,length(pick_no));
K_back     = zeros(maxnpts,length(pick_no));

for i=1:length(pick_no)
  ni=num2str(pick_no(i));
  ni=num2str(i);
  eval(['f(1:length(f' ni '),' ni ')=f' ni ';'])
  eval(['t_rawdata(1:length(f' ni '),' ni ')=t_rawdata' ni ';'])
  eval(['k_co2weak(1:length(f' ni '),' ni ') = k_co2weak' ni ';'])
  eval(['nitrogen(1:length(f' ni '),' ni ') = nitrogen' ni ';'])
  eval(['clear f' ni ' t_rawdata' ni ' k_co2weak' ni])
  end

%whos
%pause

clear estring fn i junk ni stop

global betafudge docfudge
global beta beta_pure beta_for bsm duration frequency_shift fudge
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor pick_no
global trans_amplr populationr population_tr t_rawdata;
global trans_amplp populationp population_tp strenpt w_forp w_selfp
global voi_back K_back strenrt w_forr w_selfr

%% GLOBAL VARIABLES
btz=1.4387863;B0=0.4;temperature_ref=273.15;pressure_ref=1;density=2.6867e19;
Boltzmann=1.380662e-23;mass_proton=1.6726485e-27;mass_CO2=44*mass_proton;
speed_light=2.99792458e8;

%%%%%%%%%%%%%%%%%%%%% do the pipi and deltdelt !!!!!!!!! %%%%%%%%%%%%
pipi_deltdelt_R = zeros(size(fr));  %%%%%!!NEW
temp_pipi_deltdelt_R = zeros(size(fr));  %%%%%!!NEW

for i=1:length(pick_no)
  ni=num2str(i);
  eval(['blah = k_co2weak(:,' ni ') +  nitrogen(:,' ni ');']);
  blah=blah + pipi_deltdelt_R(:,i);
  eval(['K_back(:,' ni ')=blah;']);
  end

%%%%%%%%%%%%%%%%% fit the strong SigSig stuff for 2351 band %%%%%%%%%%%%
clear topts
topts.band = [2351];
topts.PQRallowed = [-11 -12 -13 +11 +12 +13];

clear sopts;
sopts = optimset('LargeScale', 'off', ...
                 'MaxIter', 400, ...
                 'MaxFunEvals', 400, ...
                 'TolFun', 1e-6, ...
                 'TolX', 1e-6, ...
                 'DiffMinChange', 1e-9, ...
                 'Display', 'iter',...
                 'LevenbergMarquardt','on');

X0 = [1.0 1.0 0.0];
X0 = [0.2 1.0 0.0];
betafudge = X0(1);
docfudge  = X0(2);
freqshift = X0(3);
labdata.f     = f(:,1);
labdata.data  = t_rawdata;
sim.backgnd   = blah;
         %%%kcarta simulation at 0.0025 spacing
sim.pp_dd_KC  = temp_pipi_deltdelt_R;    
sim.weakco2   = zeros(size(sim.pp_dd_KC)); 
sim.n2_KC     = zeros(size(sim.weakco2));
         %%%kcarta simulation at 0.0025 spacing
sim.ipfile    = ipfileIN;

plot(labdata.f,labdata.data); title('labdata');

hahaha = sim.pp_dd_KC + sim.weakco2 + sim.n2_KC;

pooh0 = fit_run7(X0,topts,labdata,sim);
X = lsqnonlin(@fit_run7,X0,[],[],sopts,topts,labdata,sim);
% X = [0.358935 1.175938 -0.000021 ]   %%%file 1
% X = [0.200163 0.998898 -0.000039 ]   %%%file 2
poohF = fit_run7(X,topts,labdata,sim);

saver

%%%% to load the saved data and compare do run display_fits.m

disp('*******************************************************************');
disp('now you have to do the opposite (for general run7co2) ...');
disp('    (a) set dofudge = -1 in PR_sigsig/driver4um');
disp('    (b) set dofudge = -1 in PR_sigsig/yrun_sigsig');
disp('    (c) set co2_param.m -> co2_param_AIRS2002.m');
disp('    (d) set iWeight = +1 in PR_sigsig/driver4um');
disp('    (e) set dofudge = -1 in CO2_COMMON/co2_param ');
disp('*******************************************************************');
disp('*******************************************************************');
disp('*******************************************************************');
disp('optimal setting seem to be ...');
disp('    (a) set dofudge = -1 in PR_sigsig/driver4um');
disp('    (b) set dofudge = -1 in PR_sigsig/yrun_sigsig');
disp('    (c) set co2_param.m -> co2_param_AIRS2002.m');
disp('    (d) set iWeight = *** -1 ***  in PR_sigsig/driver4um ');
disp('    (e) set dofudge = -1 in CO2_COMMON/co2_param ');
disp('*******************************************************************');


