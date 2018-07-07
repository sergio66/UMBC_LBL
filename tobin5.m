function xcont = tobin5(freq,temperature,press,partpress,...
                     GasAmt,selfmult,formult,jj)
%%%% function xcont = tobin5(freq,temperature,press,partpress,...
%%%%                     GasAmt,selfmult,formult,jj)

q   = GasAmt(jj);
tmp = temperature(jj);
p   = press(jj);
pp  = partpress(jj);

array = [1 p pp tmp q];
theprofname = '/home/sergio/SPECTRA/IPFILES/water_mst';
theprofname = './water_mst';
fid=fopen(theprofname,'w');
fprintf(fid,'%4i    %12.8e  %12.8e  %12.8f  %12.8e \n',array');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path,'/home/sergio/SPECTRA/CKDLINUX');

fname = '/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDSelf5.bin';
fname = '/asl/data/kcarta/ckd/CKDSelf5.bin';
[kself,ff,tt] = contread(fname);

fname = '/asl/data/kcarta/KCARTADATA/General/CKDieee_le/CKDFor5.bin';
fname = '/asl/data/kcarta/ckd/CKDFor5.bin';
[kforn,ff,tt] = contread(fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now compute the self and foreign continuum optical depths
c2 = 1.4387863;
AVOG = 6.022045E+26;

  temper = temperature(jj);
  p      = press(jj);
  pp     = partpress(jj);
  q      = GasAmt(jj); 

  %%% find where temp lies between
  kkM = find(tt <= temper); kkM = kkM(length(kkM)); tM = tt(kkM);
  kkP = find(temper < tt);  kkP = kkP(1);           tP = tt(kkP);
  
  %%self coeffs
  ksM = interp1(ff,kself(kkM,:),freq);
  ksP = interp1(ff,kself(kkP,:),freq);
  slope = (ksP-ksM)/(tP-tM);
  ks = slope*(temper-tM) + ksM;

  %%self coeffs
  kfM = interp1(ff,kforn(kkM,:),freq);
  kfP = interp1(ff,kforn(kkP,:),freq);
  slope = (kfP-kfM)/(tP-tM);
  kf = slope*(temper-tM) + kfM;
 
  tempself   = AVOG*q*freq.*tanh(c2*freq/2/tmp)*(296/tmp);
  mstself    = ks.*tempself*(pp)*selfmult; 

  tempforn   = AVOG*q*freq.*tanh(c2*freq/2/tmp)*(296/tmp);
  mstforn    = kf.*tempforn*(p-pp)*formult; 

  xcont = mstself + mstforn;

%plot(freq,xcont)
%whos freq xcont press
%xcont = xcont';

%disp('removing junk')
eval(['!rm ./water_mst']);