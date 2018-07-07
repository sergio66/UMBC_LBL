function high_res = do_hartmann(fmin,fmax,ffin,nbox,profile,ii,NIFx,iVersHartmann);

iUse = (fmax-fmin)/ffin + 1;
iUse = round(iUse);
iUse = 1 : iUse;

%{
thedir = pwd;
xstr1 = num2str(floor(100000000*rand));
xstr2 = num2str(floor(100000000*rand));
xstr3 = num2str(floor(100000000*rand));
xstr = [xstr1 '.' xstr2 '.' xstr3];
fnameIN   = [thedir '/hartmann_input' xstr '.in'];
fnameOUT  = [' ''' thedir '/hartmann_input' xstr '.out'' '];
fnameOUT1 = [thedir '/hartmann_input' xstr '.out'];
fnameUGH  = [thedir '/hartmann_input' xstr '.ugh'];
%}

fnameINX = mktemp('hartmann_input');
fnameIN   = [fnameINX '.in'];
fnameOUT  = [' ''' fnameINX '.out'' '];
fnameOUT1 = [fnameINX '.out'];
fnameUGH  = [fnameINX '.ugh'];

fnameOUT  = ['''' fnameOUT1 '''']

fprintf(1,'fnameINX = %s \n',fnameINX);
fprintf(1,'fnameIN  = %s \n',fnameIN);
fprintf(1,'fnameOUT = %s \n',fnameOUT);
fprintf(1,'fnameOUT1 = %s \n',fnameOUT1);
fprintf(1,'fnameUGH = %s \n',fnameUGH);

%data  = [fmin-ffin*dn fmax-(nbox-dn)*ffin ffin profile(ii,:)]
dn = round((nbox-1)/2);
%c      print *,'units : cm-1  cm-1  cm-1  atm   atm   K    kilomoles/cm2 ' 
%c      print *,'Enter : sgmin sgmax dsg   pTot  pCO2  Temp qamt          : ' 
fid = fopen(fnameIN,'w');
if iVersHartmann == 2007
  %c  print *,'units : cm-1  cm-1  cm-1  atm   atm   K    kilomoles/cm2 ' 
  %c  print *,'Enter : sgmin sgmax dsg   pTot  pCO2  Temp qamt          : ' 
  data  = [fmin fmax ffin profile(ii,:)];
  fprintf(fid,'%10.6e %10.6e %10.6e %10.6e %10.6e %10.6f %10.6e \n',data);
  fprintf(fid,'%s \n',fnameOUT);
  fprintf(fid,'%3i %8i\n',NIFx,length(iUse));
elseif iVersHartmann == 2010
  %c  print *,'units : cm-1  cm-1  cm-1  atm   atm   atm K    kilomoles/cm2 ' 
  %c  print *,'Enter : sgmin sgmax dsg   pTot  pCO2  pWV Temp qamt          : ' 
  profx = profile(ii,:);
  prof(1:2) = profx(1:2);
  prof(3)   = 0;            %% use VMR of water vapor = 0 for 4,15 um?????
  prof(4:5) = profx(4:5);
  data  = [fmin fmax ffin profx];
  fprintf(fid,'%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6f %10.6e\n',data);
  fprintf(fid,'%s \n',fnameOUT);
  fprintf(fid,'%3i %8i\n',NIFx,length(iUse));
else
  iVersHartmann
  error('only have 2007 and 2010 versions of Hartmanns codes')
end 

fclose(fid);
morer = ['!more ' fnameIN]; eval(morer)

if iVersHartmann == 2007
  cder = ['!cd /home/sergio/SPECTRA/JMHARTMANN/LM_PQR_CO2_2.0/Source/'];
  eval(cder);
  %runcode = ['!/home/sergio/SPECTRA/JMHARTMANN/LM_PQR_CO2_2.0/Source/'];
  runcode = ['./loop_code_umbclbl_fast.x < ' fnameIN ' >& ' fnameUGH];
  eval(runcode);
elseif iVersHartmann == 2010
  cder = ['!cd /home/sergio/SPECTRA/JMHARTMANN_2010/Source/'];
  eval(cder);
  %runcode = ['!/home/sergio/SPECTRA/JMHARTMANN_2010/Source/'];
  %runcode = [runcode 'loop_code_umbclbl_fast.x < ' fnameIN ' >& ' fnameUGH];
  runcode = ['!cd /home/sergio/SPECTRA/JMHARTMANN_2010/Source/; ./loop_code_umbclbl_fast.x < ' fnameIN ' >& ' fnameUGH];
  eval(runcode);
else
  iVersHartmann
  error('only have 2007 and 2010 versions of Hartmanns codes')
end

cder = ['!cd ' pwd];
eval(cder);

ee = exist(fnameOUT);
ee1 = exist(fnameOUT1);
fprintf(1,'ee ee1 = %2i %2i \n',ee,ee1)
if ee == 0 & ee1 == 0
  error('WHAT???? neither fnameOUT nor fnameOUT1 exist!!!')
elseif ee1 > 0
  disp('reading in fnameOUT1, which we expect!!!')
  high_res = load(fnameOUT1);
elseif ee > 0
  disp('reading in fnameOUT, which surprises us!!!')
  high_res = load(fnameOUT);
end

[m,n]=size(high_res);
if m ~=  length(iUse)
  [m length(muse)]
  error('m ~=  length(iUse)')
end

high_res = high_res(iUse,2);

rmer = ['!/bin/rm ' fnameIN ' ' fnameOUT1 ' ' fnameUGH ' ' fnameINX];
eval(rmer);
