%% same as find_qnew_isotopes except it makes one f77 code call per isotope
%% as opposed to find_qnew_isotopes which makes does only ONE call to
%% F77 code for all isotopes

[x,y] = size(mass_info);

iUseDatabase = +1;
if iUseDatabase < 0
  %% go ahead and compute Q at 296 K
  for ii = 1 : x
    %rstrIN =[num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
    %rstrOUT=[num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
    %infile  = ['new_bt_f77_in.' rstrIN];
    %outfile = ['new_bt_f77_out.' rstrOUT];
    infile  = mktemp('new_bt_f77_in.');
    outfile = mktemp('new_bt_f77_out.');

    %%      print *,'Enter iMol iISO  iG  rT : ' 
    theT    = 296.0;
    junkstr = [num2str(gasID) ' ' num2str(ii) ' ' num2str(mass_info(ii,2))];
    junkstr = [junkstr ' ' num2str(theT)];
    fid = fopen(infile,'w');
    fprintf(fid,'%s \n',junkstr);
    fclose(fid);

    new_bt_f77 = ['!/home/sergio/SPECTRA/Global_Data_HITRAN2004/BD_TIPS_2003_try.x '];
    new_bt_f77 = [new_bt_f77 ' < ' infile ' > ' outfile];
    eval(new_bt_f77)

    qt_x = textread(outfile,'%s');
    qt_x = qt_x{7};
    qt_x = str2num(qt_x);

    rmer = ['!/bin/rm ' infile ' ' outfile]; eval(rmer);
    Q296(ii) = qt_x;
    end
else
  Q296 = mass_info(:,3)';
  end

for ii = 1 : x
  rstrIN  = [num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
  rstrOUT = [num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
  infile  = ['new_bt_f77_in.' rstrIN];
  outfile = ['new_bt_f77_out.' rstrOUT];

  %%      print *,'Enter iMol iISO  iG  rT : ' 
  theT    = tempr;
  junkstr = [num2str(gasID) ' ' num2str(ii) ' ' num2str(mass_info(ii,2))];
  junkstr = [junkstr ' ' num2str(theT)];
  fid = fopen(infile,'w');
  fprintf(fid,'%s \n',junkstr);
  fclose(fid);

  new_bt_f77 = ['!/home/sergio/SPECTRA/Global_Data_HITRAN2004/BD_TIPS_2003_try.x '];
  new_bt_f77 = [new_bt_f77 ' < ' infile ' > ' outfile];
  eval(new_bt_f77)

  qt_x = textread(outfile,'%s');
  qt_x = qt_x{7};
  qt_x = str2num(qt_x);

  rmer = ['!/bin/rm ' infile ' ' outfile]; eval(rmer);
  Qt(ii) = qt_x;
  end

% check to see that Q296 == mass_info(:,3)
% ratioCalc2Database = (Q296'./mass_info(:,3))'

% Evaluate partition functions at desired temperature and 296K 
%Qt   = a1 + b1*T   + c1*T^2   + d1*T^3; 
%Q296 = a1 + b1*296.0 + c1*(296.0^2) + d1*(296.0^3); 
%qfcn = Q296./Qt;

qfcnALL = Q296./Qt;
