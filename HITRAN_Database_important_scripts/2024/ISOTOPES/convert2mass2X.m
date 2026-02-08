%data = load('mass20_0.dat');
%data = load('mass24_0.dat');
%ngas = 55;    %% should be able to find this

dataIN = load('molparam.txt');
ngas = 61;    %% should be able to find this

boo = find(dataIN(:,3) == 0 & dataIN(:,4) == 0 & dataIN(:,5) == 0);
if length(boo) ~= ngas
  fprintf(1,'you said there are %3i gases but only found %3i gases in mass20_0.dat \n',ngas,length(boo))
  error('bah try again')
else
  fprintf(1,'found %3i gases \n',ngas);
end

%% now check the isotopes are correct
[m,n] = size(dataIN);
numiso = sum(dataIN(boo,2));

ok = m-numiso;
if ok == ngas
  fprintf(1,'you claim there are %3i gases, %3i isotopes and there are %4i gasID/iso and %4i info lines GOOD \n',ngas,numiso,ok,m-length(boo));
else  
  dataIN(boo,[ 1 2])
  fprintf(1,'you claim there are %3i gases, %3i isotopes and there are %4i gasID/iso and %4i info lines BADD \n',ngas,numiso,ok,m-length(boo));
  error('try again');
end  

disp('done cheking molparam.txt ... now proceeding to make massXY.dat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout_fid = fopen('mass24.dat','w');
fprintf(fout_fid,'%% see /asl/data/hitran/H202X/ISOTOPES/convert2mass2Xm or /umbc/xfs3/strow/asl/rta/hitran/H202X/ISOTOPES/convert2mass2X.m\n');
fprintf(fout_fid,'%% gasID   numiso   junk1       junk2        junk3     junk4 \n');

%% write out [gasID numISO -1.000000e+00 0.000000e+00 0.000000e+00 0000]
for ii = 1 : length(boo)
  if ii == dataIN(boo(ii),1)
    fprintf(fout_fid,'%6i  %6i -1.000000e+00 0.000000e+00 0.000000e+00 0000\n',dataIN(boo(ii),1),dataIN(boo(ii),2));
  else
    [ii dataIN(boo(ii),1)]
    error('gasID');
  end
end

%% write out [gasID    mol.mass    abundanace   gj    Q(296) iso]
fprintf(fout_fid,'%% gasID   mol.mass    abundanace     gj    Q(296)       iso \n');
iCnt = 0;
for ii = 1 : length(boo)
  if ii < ngas
    booII = boo(ii):boo(ii+1);
    booII = booII(2:end-1);    
  else
    booII = boo(ii):m;
    booII = booII(2:end);    
  end
  booII;
  iCnt = iCnt + length(booII);
  junkdataOUT = dataIN(booII,:);
  junkdataOUT = junkdataOUT(:,[5 2 4 3 1]);
  junkdataOUT = [ones(length(booII),1)*ii junkdataOUT];
  fprintf(fout_fid,'%6i  %10.5f   %8.5e   %3i   %8.5e   %4i \n',junkdataOUT');
end
fclose(fout_fid);

if iCnt ~= numiso
  error('iCnt ~= numiso')
else
  fprintf(1,'done dumping out %4i isotopes for %3i gases \n',numiso,ngas);
end  
disp('copy to /home/sergio/git/UMBC_LBL/MASS_ISOTOPES')
