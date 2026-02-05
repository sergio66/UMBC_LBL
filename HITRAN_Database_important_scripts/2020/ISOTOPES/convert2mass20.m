data = load('mass20_0.dat');

ngas = 55;    %% should be able to find this
boo = find(data(:,3) == 0 & data(:,4) == 0 & data(:,5) == 0);
if length(boo) ~= ngas
  fprintf(1,'you said there are %3i gases but only found %3i gases in mass20_0.dat \n',ngas,length(boo))
  error('bah try again')
else
  fprintf(1,'found %3i gases \n',ngas);
end

%% now check the isotopes are correct
[m,n] = size(data);
numiso = sum(data(boo,2));

ok = m-numiso;
if ok == ngas
  fprintf(1,'you claim there are %3i gases, %3i isotopes and there are %4i gasID/iso and %4i info lines GOOD \n',ngas,numiso,ok,m-length(boo));
else  
  data(boo,[ 1 2])
  fprintf(1,'you claim there are %3i gases, %3i isotopes and there are %4i gasID/iso and %4i info lines BADD \n',ngas,numiso,ok,m-length(boo));
  error('try again');
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('mass20.dat','w');
fprintf(fid,'%% see /asl/data/hitran/H2020/ISOTOPES/convert2mass20.m \n');
fprintf(fid,'%% gasID   numiso   junk1       junk2        junk3     junk4 \n');

%% write out [gasID numISO -1.000000e+00 0.000000e+00 0.000000e+00 0000]
for ii = 1 : length(boo)
  if ii == data(boo(ii),1)
    fprintf(fid,'%6i  %6i -1.000000e+00 0.000000e+00 0.000000e+00 0000\n',data(boo(ii),1),data(boo(ii),2));
  else
    [ii data(boo(ii),1)]
    error('gasID');
  end
end

%% write out [gasID    mol.mass    abundanace   gj    Q(296) iso]
fprintf(fid,'%% gasID   mol.mass    abundanace     gj    Q(296)       iso \n');
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
  junkdata = data(booII,:);
  junkdata = junkdata(:,[5 2 4 3 1]);
  junkdata = [ones(length(booII),1)*ii junkdata];
  fprintf(fid,'%6i  %10.5f   %8.5e   %3i   %8.5e   %4i \n',junkdata');
end
fclose(fid);
if iCnt ~= numiso
  error('iCnt ~= numiso')
else
  fprintf(1,'done dumping out %4i isotopes for %3i gases \n',numiso,ngas);
end  
