%% same as find_qnew_isotopes_slow except it does only ONE call to
%% F77 code for all isotopes, as opposed to find_qnew_isotopes_slow which 
%% makes one f77 code call per isotope

[x,y] = size(mass_info);

thisdir = pwd;
thisdir(length(thisdir)+1) = '/';

%% replace  num2str(round(rand(1,1)*1e9)) w/ num2str(floor(100000000*rand))

iUseDatabase = +1;
if iUseDatabase < 0
  %% go ahead and compute Q at 296 K for all isotopes
  rstrIN = [num2str(floor(100000000*rand)) '.' num2str(floor(100000000*rand))];
  rstrOUT= [num2str(floor(100000000*rand)) '.' num2str(floor(100000000*rand))];
  infile  = [thisdir 'new_bt_f77_in.' rstrIN];
  outfile = [thisdir 'new_bt_f77_out.' rstrOUT];

  %%      print *,'Enter iMol iNumISO  rT : ' 
  theT    = 296.0;
  junkstr = [num2str(gasID) ' ' num2str(x) ' ' num2str(theT)];
  fid = fopen(infile,'w');
  fprintf(fid,'%s \n',junkstr);
  fclose(fid);

  new_bt_f77 = ['!/home/sergio/SPECTRA/Global_Data_HITRAN2004/'];
  new_bt_f77 = [new_bt_f77 'BD_TIPS_2003_allisotopes.x '];
  new_bt_f77 = [new_bt_f77 ' < ' infile ' > ' outfile];
  eval(new_bt_f77)

  qt_xtemp = textread(outfile,'%s');
  for ii = 1 : x
    qt_xx = qt_xtemp{5+ii};
    qt_xx = str2num(qt_xx);
    Q296(ii) = qt_xx;
    end
  rmer = ['!/bin/rm ' infile ' ' outfile]; eval(rmer);
else
  Q296 = mass_info(:,3)';
  end

%% go ahead and compute Q at "tempr" K for all isotopes
rstrIN  = [num2str(floor(100000000*rand)) '.' num2str(floor(100000000*rand))];
rstrOUT = [num2str(floor(100000000*rand)) '.' num2str(floor(100000000*rand))];
infile  = [thisdir 'new_bt_f77_in.' rstrIN];
outfile = [thisdir 'new_bt_f77_out.' rstrOUT];

%%      print *,'Enter iMol iNumISO  rT : ' 
theT    = tempr;
junkstr = [num2str(gasID) ' ' num2str(x) ' ' num2str(theT)];
fid = fopen(infile,'w');
fprintf(fid,'%s \n',junkstr);
fclose(fid);

new_bt_f77 = ['!/home/sergio/SPECTRA/Global_Data_HITRAN2004/'];
new_bt_f77 = [new_bt_f77 'BD_TIPS_2003_allisotopes.x '];
new_bt_f77 = [new_bt_f77 ' < ' infile ' > ' outfile];
eval(new_bt_f77)

qt_xtemp = textread(outfile,'%s');
for ii = 1 : x
  qt_xx = qt_xtemp{5+ii};
  qt_xx = str2num(qt_xx);
  Qt(ii) = qt_xx;
  end
rmer = ['!/bin/rm ' infile ' ' outfile]; eval(rmer);

% check to see that Q296 == mass_info(:,3)
% ratioCalc2Database = (Q296'./mass_info(:,3))'

% Evaluate partition functions at desired temperature and 296K 
%Qt   = a1 + b1*T   + c1*T^2   + d1*T^3; 
%Q296 = a1 + b1*296.0 + c1*(296.0^2) + d1*(296.0^3); 
%qfcn = Q296./Qt;

qfcnALL = Q296./Qt;
