%%%%%%%%%%%%%%%%%%%  loader.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chosen1; chosen1NEW;
chosen2; chosen2NEW; %chosen2NEWNEW;
chosen3; chosen3NEW; 
chosen4; chosen4NEW; 
chosen5;
chosen6;
chosen7;

format long e

disp('------------------------------------------------------------------')
disp(' #  File Name  T(K)     P(CO2)  P(tot) P.L.(m)  REGION   resolution')
disp(' co2a :                  HPa     HPa                       cm-1');
disp('------------------------------------------------------------------')
disp(' 1  CDAB7S     296.0    26.80    26.80  512.746    4um     0.0038')
disp(' 2  CDAB8S     296.0    26.80   161.40  512.746    4um     0.0038')
disp(' 3  CDAB9S     296.0    26.80   560.30  512.746    4um     0.0038')
disp(' 4  CDAB10S    296.0    26.80   961.70  512.746    4um     0.0038')
disp(' 5  CDBC1S     296.0     0.24    0.24   32.746     4um    0.0038')
disp(' 6  CDBC2S     296.0     0.20    0.20   32.746     4um    0.0038')
disp(' 7  CDBC3S     296.0     0.20    0.20   32.746     4um    0.0038')
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

% note that the weak bkgnds files and the N2 files can be computed on a 
% 0.0025 spacing, using run6 and appropriate run6co2, but then the 
% spectra have to be spline interpolated onto the appropriate grid!
% this grid can easily be found by loading in eg
% /salsify/scratch3/Sergio/RAL_DATA/CO2A/file1n2.dat etc!

% TO MAKE THE N2 BACKGND, remember first file has no N2, so k=0
% [fr,k]=run6(22,2370,2490,0.0005,0.1,0.5,1,1,2,25,5,1e-28,1e-28,'V',1,...
%                '../SPECTRA/IPFILES/RAL_CO2/ral4n2');
% kk=zeros(1,48000); kk(2:4,:)=k; clear k
% load /salsify/scratch3/Sergio/RAL_DATA/CO2A/file1n2.dat
% kki=interp1(fr,kk',file1n2(:,1));
% d=[file1n2(:,1) kki(:,1)]; whos; plot(d(:,1),exp(-d(:,2)))
% d=[file2n2(:,1) kki(:,2)]; whos; plot(d(:,1),exp(-d(:,2)))
% d=[file3n2(:,1) kki(:,3)]; whos; plot(d(:,1),exp(-d(:,2)))
% d=[file4n2(:,1) kki(:,4)]; whos; plot(d(:,1),exp(-d(:,2)))
% save /salsify/scratch3/Sergio/RAL_DATA/CO2A/file1n2_avg.dat d -ascii -double
% or to use NEW N2, straight from LBLRTM v5.10 with 79/21 N2/O2 continuum 
%save blah/file1n2_avg_Apr00.dat d -ascii -double     
% or to use NEW N2, straight from LBLRTM v5.10 but more general
%save blah/file1n2_avg_Apr00_newN2.dat d -ascii -double     
% and repeat for i=1,2,3,4 so you save all 4 files!

% TO MAKE THE CO2 BACKGND, LOOK AT  removeCO2*.m AND FOLLOW DIRECTIONS
%%%%%   removeCO2*.m      main loop     linemix loop    ||    RESULT
%%%% ---------------------------------------------------------------
%% k1     ORIG              ON             OFF         ||      T
%% k2     ORIG              OFF            ON          ||      S + W
%% k3     RAL               OFF            ON          ||      S
%for n=1,2,3 as defined above, do
% [fr,kn]=run6co2(2,2380,2480,0.0005,0.1,0.5,1,1,2,150,5,1e-28,1e-28,'F',...
%                 '1','b','../SPECTRA/IPFILES/RAL_CO2/ral4');
%% from which the background lines to be used for the 4umn fit files (file1w,
%% file2w,file3w,file4w) is given by B = T+ W = k1 + (k2-k3)
%%
%% eg
%% load /salsify/scratch3/Sergio/RAL_DATA/CO2A/file3w_avg.dat
%% bb=spline(fr',b(3,:),file3w_avg(:,1)');            where b=k1+(k2-k3);
%% d=[file3w_avg(:,1) bb']; whos
%% save /salsify/scratch3/Sergio/RAL_DATA/CO2A/file3w_avg_Apr00.dat d -ascii -double

avg=1;

if (avg < 0)
  % Data parameters for CO2A   ..... these were the unaveraged ones used in 
  % the fits done prior to Oct 20, 1999
  fil(1,:)='CDAB7S';    filb(1,:)='CDFA1B';     filw(1,:)='file1w';
    filn2(1,:)='file1n2';        ps(1)=26.80; pf(1)=26.80; 
    fudge(1)=0.9999;
  fil(2,:)='CDAB8S';    filb(2,:)='CDFA1B';     filw(2,:)='file2w';
    filn2(2,:)='file2n2';      ps(2)=26.80; pf(2)=161.40;  
    fudge(2)=0.9875;
  fil(3,:)='CDAB9S';    filb(3,:)='CDFA1B';     filw(3,:)='file3w';
    filn2(3,:)='file3n2';      ps(3)=26.80; pf(3)=560.03; 
    fudge(3)=0.975;
  fil(4,:)='CDABXS';    filb(4,:)='CDFA1B';     filw(4,:)='file4w';
    filn2(4,:)='file4n2';      ps(4)=26.80; pf(4)=961.70; 
    fudge(4)=0.96;
  for i=1:4
    fitted_t(i)=296.0;
    pl(i)=512.746*100;           %change to cm as this is what run6co2 needs
    pf(i)=pf(i)-ps(i);           %foreign = total - self
    pf(i)=pf(i)/1013.25;         %change to atm as this is what run6co2 needs
    ps(i)=ps(i)/1013.25;         %change to atm as this is what run6co2 needs
    end
elseif (avg == 0)
  % Data parameters for CO2A   ..... these are the averaged ones used in 
  % the fits done Oct 20, 1999
  % the averaging is Mike Wolk going thru the parameter files to find the
  % avg temperature, pressure etc instead of nominal reported temps
  fil(1,:)='CDAB7S';    filb(1,:)='CDFA1B';     filw(1,:)='file1w_avg';
    filn2(1,:)='file1n2_avg';    
    ps(1)=26.80; pf(1)=26.80; temp(1)=18.5143;
    fudge(1)=0.9999;
  fil(2,:)='CDAB8S';    filb(2,:)='CDFA1B';     filw(2,:)='file2w_avg';
    filn2(2,:)='file2n2_avg';    
    ps(2)=26.80; pf(2)=161.40;  temp(2)=18.5143;
    fudge(2)=0.9875;
  fil(3,:)='CDAB9S';    filb(3,:)='CDFA1B';     filw(3,:)='file3w_avg';
    filn2(3,:)='file3n2_avg';    
    ps(3)=26.80; pf(3)=560.30; temp(3)=18.7129;
    fudge(3)=0.975;
  fil(4,:)='CDABXS';    filb(4,:)='CDFA1B';     filw(4,:)='file4w_avg';
    filn2(4,:)='file4n2_avg';    
    ps(4)=26.80; pf(4)=961.75;  temp(4)=18.8914;
    fudge(4)=0.96;
  for i=1:4
    fitted_t(i)=temp(i)+273.15;
    pl(i)=512.746*100;           %change to cm as this is what run6co2 needs
    pf(i)=pf(i)-ps(i);           %foreign = total - self
    pf(i)=pf(i)/1013.25;         %change to atm as this is what run6co2 needs
    ps(i)=ps(i)/1013.25;         %change to atm as this is what run6co2 needs
    end
elseif (avg > 0)
  % Data parameters for CO2A   ..... these are the averaged ones used in 
  % the fits done April 2000
  % the averaging is Mike Wolk going thru the parameter files to find the
  % avg temperature, pressure etc instead of nominal reported temps
  fil(1,:)='CDAB7S';    filb(1,:)='CDFA1B';     filw(1,:)='file1w_avg_Apr00';
    filn2(1,:)='file1n2_avg_Apr00_newN2';    
    %filn2(1,:)='file1n2_avg_Apr00';    
    ps(1)=26.80; pf(1)=26.80; temp(1)=18.5143;
    fudge(1)=0.9999;
  fil(2,:)='CDAB8S';    filb(2,:)='CDFA1B';     filw(2,:)='file2w_avg_Apr00';
    filn2(2,:)='file2n2_avg_Apr00_newN2';    
    %filn2(2,:)='file2n2_avg_Apr00';    
    ps(2)=26.80; pf(2)=161.40;  temp(2)=18.5143;
    fudge(2)=0.9875;
  fil(3,:)='CDAB9S';    filb(3,:)='CDFA1B';     filw(3,:)='file3w_avg_Apr00';
    filn2(3,:)='file3n2_avg_Apr00_newN2';    
    %filn2(3,:)='file3n2_avg_Apr00';    
    ps(3)=26.80; pf(3)=560.30; temp(3)=18.7129;
    fudge(3)=0.975;
  fil(4,:)='CDABXS';    filb(4,:)='CDFA1B';     filw(4,:)='file4w_avg_Apr00';
    filn2(4,:)='file4n2_avg_Apr00_newN2';    
    %filn2(4,:)='file4n2_avg_Apr00';    
    ps(4)=26.80; pf(4)=961.75;  temp(4)=18.8914;
    fudge(4)=0.96;
  for i=1:4
    fitted_t(i)=temp(i)+273.15;
    pl(i)=512.746*100;           %change to cm as this is what run6co2 needs
    pf(i)=pf(i)-ps(i);           %foreign = total - self
    pf(i)=pf(i)/1013.25;         %change to atm as this is what run6co2 needs
    ps(i)=ps(i)/1013.25;         %change to atm as this is what run6co2 needs
    end
  end

% Data parameters for CO2A  15 um data  ..... these are the unaveraged ones 
% used in the fits done prior to Oct 20, 1999. but since the 15 um data is
% too optically thick to do any worthwhile fits, I haven't bothered putting
% in the averaged values.
if avg < 0
  fil(5,:)='CDBC1S';    filb(5,:)='CDEB2B';     filw(5,:)='file5w';
    filn2(5,:)='file1n2';      ps(10)=00.20; pf(10)=00.20; 
    fudge(5)=1.0;
  fil(6,:)='CDBC2S';    filb(6,:)='CDEB2B';     filw(6,:)='file6w';
    filn2(6,:)='file1n2';      ps(11)=00.20; pf(11)=00.20;  
    fudge(6)=1.0;
  fil(7,:)='CDBC3S';    filb(7,:)='CDEB2B';     filw(7,:)='file7w';
    filn2(7,:)='file1n2';      ps(12)=13.42; pf(12)=13.42; 
    fudge(7)=1.0;
  for i=5:7
    fitted_t(i)=296.0;
    pl(i)=32.746*100;            %change to cm as this is what run6co2 needs
    pf(i)=pf(i)-ps(i);           %foreign = total - self
    pf(i)=pf(i)/1013.25;         %change to atm as this is what run6co2 needs
    ps(i)=ps(i)/1013.25;         %change to atm as this is what run6co2 needs
    end
  end

for i=1:length(pick_no);
  temperature(i)=296.0;
  temperature(i)=fitted_t(pick_no(i));    %new on Oct 20, 1999
  fil_nam(i,:)=fil(pick_no(i),:);         %transmittance data with gas
  fil_nam_back(i,:)=filb(pick_no(i),:);   %transmittance data w/o gas
  fil_nam_weak(i,:)=filw(pick_no(i),:);   %weak background lines from run6co2
  fil_nam_n2(i,:)=filn2(pick_no(i),:);    %N2 background lines from run6
  path_length(i)=pl(pick_no(i));
  pressure_self(i)=ps(pick_no(i));
  pressure_for(i)=pf(pick_no(i));
  end

clear fil filb filn2 filw fitted_t pl ps pf

% Load in raw data and background calcs
maxnpts=0;
for i=1:length(pick_no)
  fprintf(1,'loading in RAL data ... yawn! ... \n');
  fn=fil_nam(i,:);
  fnb=fil_nam_back(i,:);
  fnw=fil_nam_weak(i,:);
  fn2=fil_nam_n2(i,:);

  ni=num2str(i);

  % the stuff used to be on /salsify/scratch3/Sergio/RAL_DATA/CO2A+-CO2B
  % and was all in lowercase. now it is in uppercase eg 
  % cdab7s.dat  --> CDAB7S.DAT 
  if (pick_no(i) <= 4)
    estring='load /asl/data/ral/Co2a/asciidat/'; 
          eval([estring fn '.DAT']);
    estring='load /asl/data/ral/Co2a/asciidat/'; 
          eval([estring fnb '.DAT']);
    estring='load /asl/data/ral/Co2a/asciidat/'; 
          eval([estring fnw '.DAT']);
    estring='load /asl/data/ral/Co2a/asciidat/'; 
          eval([estring fn2 '.DAT']);
  elseif (pick_no(i) >4)
    estring='load /asl/data/ral/Co2b/asciidat/'; 
          eval([estring fn '.DAT']);
    estring='load /asl/data/ral/Co2b/asciidat/'; 
          eval([estring fnb '.DAT']);
    estring='load /asl/data/ral/Co2b/asciidat/'; 
          eval([estring fnw '.DAT']);
    estring='load /asl/data/ral/Co2a/asciidat/'; 
          eval([estring fn2 '.DAT']);
    end

  eval(['ff =' fn '(:,1);']);
  eval(['rawdata = ' fn '(:,2)./' fnb '(:,2);']);
  ind=find((ff >= 2380) & (ff <= 2480));
  ff=ff(ind);
  rawdata=rawdata(ind);

  switch pick_no(i)
    case{1}
      chosen=chosen11; 
    case{2}
      chosen=chosen22; 
    case{3}
      chosen=chosen33; 
    case{4}
      chosen=chosen44; 
    case{5}
      chosen=chosen55; 
    case{6}
      chosen=chosen66; 
    case{7}
      chosen=chosen77; 
    otherwise
      err('unknown case!')
    end

  fdg=fudge(pick_no(i))
  %%%fdg=1.0

  eval(['f' ni '= ff(chosen);']);
  eval(['t_rawdata' ni '= rawdata(chosen)*fdg;']);
  eval(['K_backk' ni '=' fnw '(chosen,2);']);
  eval(['nitrogen' ni '=' fn2 '(chosen,2);']);

  eval(['maxnpts=max([ maxnpts length(ff(chosen))]);']);
  disp(['f, t_rawdata  ' ni ' loaded in.']);
  eval(['clear lor lor_mod ' fn])
  disp(['f, t_rawdata, K_backk: ' ni ' loaded in.'])
  eval(['length(f' ni ')'])

  eval(['gah =' fnw '(chosen,2);']);
  plot(ff(chosen),rawdata(chosen),ff(chosen),exp(-gah),'r');
  title('co2 background lines (r), co2 data (b)')

  red_num=input('Reduce # data points by factor of : ','s');
  eval(['[f' ni ',t_rawdata' ni ',K_backk' ni ',nitrogen' ni ']=reduce4(f' ni ',t_rawdata' ni ',K_backk' ni ', nitrogen' ni ',' red_num ');'])

  eval(['maxnpts=min([ maxnpts length(f' ni ')]);']);

  end

clear ff

f         = zeros(maxnpts,length(pick_no));
t_rawdata = zeros(maxnpts,length(pick_no));
K_backk   = zeros(maxnpts,length(pick_no));
K_back    = zeros(maxnpts,length(pick_no));

for i=1:length(pick_no)
  ni=num2str(i);
  eval(['f(1:length(f' ni '),' ni ')=f' ni ';'])
  eval(['t_rawdata(1:length(f' ni '),' ni ')=t_rawdata' ni ';'])
  eval(['K_backk(1:length(f' ni '),' ni ') = K_backk' ni ';'])
  eval(['nitrogen(1:length(f' ni '),' ni ') = nitrogen' ni ';'])
  eval(['clear f' ni ' t_rawdata' ni ' K_backk' ni])
  end

%whos
%pause

clear estring fn i junk ni stop

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
cd PIPI
disp('doing the R PIPI branch ... ');
run1_R;

cd ../DELTDELT
disp('doing the R DELTDELT branch ... ');
run1_R;

cd ..

for i=1:length(pick_no)
  ni=num2str(i);
  eval(['blah = K_backk(:,' ni ') +  nitrogen(:,' ni ');']);
  blah=blah + pipiR_lor_cousin(:,i) + dedeR_lor_cousin(:,i);
  eval(['K_back(:,' ni ')=blah;']);
  end


%%%%%%%%%%%%%%%%%%%%% Load in hittomat-produced file %%%%%%%%%%%%%%
%%%load /beet/users/tobin/Co2pr/Hit_co2pr/hit_43_92.mat
load ../Hit_co2pr/hit_43_92.mat

%%%%%%%%%%%%%%%%%%%%% Sort out lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('sorting out foreground lines')
v_l=1;v_u=9;isotope=1;prb='R';
index=find(j_lower(:,5)==prb & v_lower==v_l & v_upper==v_u & iso==isotope);

% Now select out these lines from the rest.
j_lowerr(:,1:9)=j_lower(index,1:9);		% J of lower state (strings)
junk_str=j_lowerr(:,6:8);
for n=1:length(junk_str);
  eval(['jr(n)=' junk_str(n,:) ';']);
  end
jr=jr';
freqr=freq(index);
strenr=stren(index);
elowerr=elower(index);
w_forr=w(index);
w_fortempr=w_temp(index);
w_selfr=w_s(index);

 %%%%%%%%%%%%%%%%%% Correct widths to desired temperature %%%%%%%%%%%%%%
disp('correct widths to desired temperature')
for i=1:length(pick_no)
  w_fr(:,i)=w_forr.*(296/temperature(i)).^w_fortempr;
  w_sr(:,i)=w_selfr.*(296/temperature(i))^0.685;
  end
w_forr=w_fr;w_selfr=w_sr;clear w_fr w_sr

%%%%%%%% Flip arrays (if needed) so that lowest J's appear first  %%%%%%%%%% 
if jr(1)>jr(2)
  j_lowerr=flipud(j_lowerr);freqr=flipud(freqr);strenr=flipud(strenr);
  elowerr=flipud(elowerr);w_selfr=flipud(w_selfr);jr=flipud(jr);
  w_forr=flipud(w_forr);
  end
 
%%%%%%%%%%%%%%% Adjust strengths to the correct temperature %%%%%%%%%%%%
% NOTE that the stimulated emiision terms ARE included.
disp('adjusting strengths to desired temperature')
        a1=-0.2199475485E+01;   b1=0.9675055715E+00;
        c1=-0.8082711378E-03;   d1=0.2803987451E-05;
for i=1:length(pick_no)
  Qt=a1 + b1*temperature(i) + c1*temperature(i)^2 + d1*temperature(i)^3;
  Q296=a1 + b1*296 + c1*296^2 + d1*296^3;
  numer=Q296*exp(btz*elowerr/296).*(1-exp(-btz*freqr/temperature(i)));
  denom=Qt*exp(btz*elowerr/temperature(i)).*(1-exp(-btz*freqr/296));
  strenrt(:,i)=strenr.*numer./denom;
  end

%%%%%%%%%%%%%%%%%%%%% Sort out lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('sorting out foreground lines')
v_l=1;v_u=9;isotope=1;prb='P';
index=find(j_lower(:,5)==prb & v_lower==v_l & v_upper==v_u & iso==isotope);
% Now select out these lines from the rest.
j_lowerp(:,1:9)=j_lower(index,1:9);		% J of lower state (strings)
junk_str=j_lowerp(:,6:8);
for n=1:length(junk_str);
  eval(['jp(n)=' junk_str(n,:) ';']);
  end
jp=jp';freqp=freq(index);strenp=stren(index);elowerp=elower(index);
w_forp=w(index);w_fortempp=w_temp(index);w_selfp=w_s(index);

%%%%%%%%%%%%%%%%%% Correct widths to desired temperature %%%%%%%%%%%%%%
disp('correct widths to desired temperature')
for i=1:length(pick_no)
  w_fp(:,i)=w_forp.*(296/temperature(i)).^w_fortempp;
  w_sp(:,i)=w_selfp.*(296/temperature(i))^0.685;
  end

w_forp=w_fp;w_selfp=w_sp;clear w_fp w_sp

%%%%%%%% Flip arrays (if needed) so that lowest J's appear first  %%%%%%%%%% 
if jp(1)>jp(2) 
  j_lowerp=flipud(j_lowerp);freqp=flipud(freqp);strenp=flipud(strenp);
  elowerp=flipud(elowerp);w_selfp=flipud(w_selfp);jp=flipud(jp);
  w_forp=flipud(w_forp);
  end

for i=1:length(pick_no)
  Qt=a1 + b1*temperature(i) + c1*temperature(i)^2 + d1*temperature(i)^3;
  Q296=a1 + b1*296 + c1*296^2 + d1*296^3;
  numer=Q296*exp(btz*elowerp/296).*(1-exp(-btz*freqp/temperature(i)));
  denom=Qt*exp(btz*elowerp/temperature(i)).*(1-exp(-btz*freqp/296));
  strenpt(:,i)=strenp.*numer./denom;
  end

%%%%%%%%%%%%%%%%%%%%% Clear unwanted variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear accuracy dipole elower freq gas_id iso j_lower j_upper line_status
clear p_shift reference stren v_lower v_upper w w_s w_temp j_upper index
clear hitfile v_l v_u mass_proton

clear w_fortempp w_fortempr voi_back prb numer n junk_str
clear j_lowerp j_lowerr isotope indexp ind estring a1 b1 c1 d1 Qt Q296



