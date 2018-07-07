%% same as find_qnew_isotopes_slow except it does only ONE call to
%% F77 code for all isotopes, as opposed to find_qnew_isotopes_slow which 
%% makes one f77 code call per isotope

addpath /asl/matlib/aslutil

otherinfoX = load('Global_Data_HITRAN2012/molparam1.txt.simple');
otherinfo  = otherinfoX(:,1);
otherinfo4 = otherinfoX(:,4);

ngas = 47;

mass00 = load('mass12.dat');
%% mass00a = mass00(01:32,1:2);
%% mass00b = mass00(33:length(mass00),:);
mass00a = mass00(01:ngas,1:2);
mass00b = mass00(ngas+1:length(mass00),:);

cumsumms = cumsum(mass00a(:,2));
cumsumms = [0; cumsumms];

iaT = 156:10:376;
iaT = 70:05:405;

fidX = fopen('qtips12.m','w');

str = ['function [a,b,c,d,g] = qtips(gasid,liso);'];
fprintf(fidX,'%s\n',str);
fprintf(fidX,' \n');

str = '% this was made by find_qnewABCD_H12.m';
fprintf(fidX,'%s\n',str);
str = '% returns nuclear degeneracy factors g and Gamache''s internal partition';
fprintf(fidX,'%s\n',str);
str = '% sum coefficients a, c, c, and d.  ';
fprintf(fidX,'%s\n',str);
str = '%C  PURPOSE        TOTAL INTERNAL PARTITION SUMS';
fprintf(fidX,'%s\n',str);
fprintf(fidX,' \n');

for ii = 1 : ngas
  strss = (cumsumms(ii)+1:cumsumms(ii+1));
  %strss = mass00b(strss,4); strss = strss';
  strss = otherinfo4(strss); strss = strss';
  strss = num2str(strss);
  str = ['if gasid == ' num2str(ii) '; g = [' strss ']; end'];
  fprintf(fidX,'%s\n',str);
end

fprintf(fidX,' \n');
str = ['g = g'';'];
fprintf(fidX,'%s\n',str);

str = ['if (length(g) ~= liso)'];
fprintf(fidX,'%s\n',str);
str = ['  fprintf(1,''there is inconsistency in number of isotopes \n''); '];
fprintf(fidX,'%s\n',str);
str = ...
  ['  fprintf(1,''in qtips.m %3i and that in mass.dat %3i \n'',length(g),liso);'];

fprintf(fidX,'%s\n',str);
str = ['  error(''please check!!!!!!'')'];
fprintf(fidX,'%s\n',str);
str = 'end';
fprintf(fidX,'%s\n',str);

fprintf(fidX,' \n');

for iGasID = 1 : ngas
  strss = (cumsumms(iGasID)+1:cumsumms(iGasID+1));
  strss = otherinfo(strss)';
  strss = num2str(strss);
  clear qfcnALL Y A B C D q296  matr
  for iT = 1 : length(iaT)

    iNumIso = mass00a(iGasID,2);

    %% go ahead and compute Q at "tempr" K for all isotopes
    %rstrIN =[num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
    %rstrOUT=[num2str(round(rand(1,1)*1e9)) '.' num2str(round(rand(1,1)*1e9))];
    %infile  = ['new_bt_f77_in.' rstrIN];
    %outfile = ['new_bt_f77_out.' rstrOUT];
    infile  = mktemp('new_bt_f77_in.');
    outfile = mktemp('new_bt_f77_out.');

    %%      print *,'Enter iMol iNumISO  rT : ' 
    theT    = iaT(iT);
    junkstr = [num2str(iGasID) ' ' num2str(iNumIso) ' ' num2str(theT)];
    fid = fopen(infile,'w');
    fprintf(fid,'%s \n',junkstr);
    fclose(fid);

    new_bt_f77 = ['!/home/sergio/SPECTRA/Global_Data_HITRAN2012/'];
    %new_bt_f77 = [new_bt_f77 'BD_TIPS_2012_allisotopes.x '];
    new_bt_f77 = [new_bt_f77 'TIPS_2012_allisotopes.x'];
    new_bt_f77 = [new_bt_f77 ' < ' infile ' > ' outfile];
    eval(new_bt_f77)

    qt_xtemp = textread(outfile,'%s');
    for ii = 1 : iNumIso
      qt_xx = qt_xtemp{6+ii};
      qt_xx = str2num(qt_xx);
      Qt(ii) = qt_xx;
    end
    rmer = ['!/bin/rm ' infile ' ' outfile]; eval(rmer);

    qfcnALL(iT,:) = Qt;
  end 

  q296 = find(mass00b(:,1) == iGasID); q296 = mass00b(q296,5);

  figure(1); clf;
  figure(2); clf;
  for iI = 1 : iNumIso
    [P,S] = polyfit(iaT,qfcnALL(:,iI)',3);
    YiaT = polyval(P,iaT);
    figure(1); plot(iaT,qfcnALL(:,iI),iaT,YiaT,'--'); hold on
    figure(2); plot(iaT,100*(qfcnALL(:,iI)-YiaT')./qfcnALL(:,iI)); hold on
    Y(iI) = polyval(P,296);
    A(iI) = P(4);
    B(iI) = P(3);
    C(iI) = P(2);
    D(iI) = P(1);
  end
  figure(1); title(num2str(iGasID)); 
  figure(2); title(num2str(iGasID)); %%pause(1);

  format short e
  popo = Y./q296';
  matr = [A;B;C;D]; matr = matr';

  str = ['% fitted/actual at 296K = '];
  fprintf(fidX,'%s \n',str);
  str = ['%' num2str(popo)];
  fprintf(fidX,'%s \n',str);

  str = ['if gasid == ' num2str(iGasID)];
  fprintf(fidX,'%s \n',str);
  str = ['    abcd = [ ...'];
  fprintf(fidX,'%s \n',str);
  for ii = 1 : iNumIso-1
    fprintf(fidX,'  %9.5e , %9.5e, %9.5e , %9.5e ; \n',matr(ii,:));
  end
  ii = iNumIso;
  fprintf(fidX,'  %9.5e , %9.5e, %9.5e , %9.5e ]; \n',matr(ii,:));
  format

  str = ['  disp(''The isotopes are ' strss ''');'];
  fprintf(fidX,'%s \n',str);
  str = ['end'];
  fprintf(fidX,'%s \n \n',str);

end

str = 'a = abcd(:,1);';   fprintf(fidX,'%s \n',str);
str = 'b = abcd(:,2);';   fprintf(fidX,'%s \n',str);
str = 'c = abcd(:,3);';   fprintf(fidX,'%s \n',str);
str = 'd = abcd(:,4);';   fprintf(fidX,'%s \n',str);

fclose(fidX);
