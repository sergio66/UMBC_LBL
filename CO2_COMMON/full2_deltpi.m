function [fullmix4]=full2_deltpi(freqq,f,W_plus,jq,w_selfq,w_forq,...
                                 temperature,trans_ampl,population_t,stuff,...
                                 birn,ratio,strenqt,layeramt,ymix1,band,prb) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  full_mix4.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file computes the full mixing absorption coefficient, 
%       using FORTRAN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global v12p1

density         = stuff.density; 
pressure_self   = stuff.pressure_self; 
pressure_for    = stuff.pressure_for; 
pressure_ref    = stuff.pressure_ref; 
temperature_ref = stuff.temperature_ref; 
path_length     = stuff.path_length; 
bsm             = stuff.bsm; 
duration        = stuff.duration;
frequency_shift = stuff.frequency_shift;

freqq_shift0 = frequency_shift/100+freqq;
freqq0        = freqq;
jq0           = jq;
w_selfq0      = w_selfq;
w_forq0       = w_forq;
trans_ampl0   = trans_ampl;
population_t0 = population_t;
strenqt0      = strenqt;
ymix10        = ymix1;

K_scale_mixing = density*pressure_self/pressure_ref*temperature_ref*...
	         path_length/temperature/pi;
no_lines       = length(freqq); 
no_pts         = length(f);
k              = zeros(no_pts,1);

freqq_shift    = freqq+frequency_shift/100;

birnbirn = 0;
if ((birn == 'b') | (birn == 'B') ) 
  birnbirn = 1;
elseif ((birn == 'c') | (birn == 'C') )  
  error('cannot have mixing AND cousin!!');
  birnbirn = -1;
  end

fprintf(1,'%4i Qdeltpi bands (or 2322 pipi): reorder stuff ... \n',band);
%plot(freqq,log10(strenqt)); title('1');
%ratio
%pause;

%%if (band == 2322 & (prb == 'p' | prb == 'P') & v12p1 < 0)
if (band == 2322 & (prb == 'p' | prb == 'P'))
  bwah = -1;      %%%%%%% HAVE THIS SET TO -1 to order things correctly,
else
  bwah = +1;          %%%%%%% HAVE THIS SET TO +1 to order things correctly
  end
plotbwah = -1;

if (ratio >= 1e-2) 
  disp('ratio >= 1e-2')
  H     = diag(freqq_shift)+sqrt(-1)*W_plus;
  [A,L] = eig(H);

  %to test no mixing ....
  %w_tot=(pressure_self*w_selfq+pressure_for*w_forq)/pressure_ref; 
  %H     = diag(freqq_shift)+sqrt(-1)*diag(w_tot);
  %[A,L] = eig(H);

  if bwah > 0
    [blah,ii]    = sort(freqq_shift);
    if ((band == 740) & (prb ~= 'R'))
      ii=1:length(ii);
      end
    freqq_shift  = freqq_shift(ii);
    w_selfq      = w_selfq(ii);
    w_forq       = w_forq(ii);
    trans_ampl   = trans_ampl(ii);
    population_t = population_t(ii);
    strenqt      = strenqt(ii);
    jq           = jq(ii);

    diagLorig = diag(L);
    diagL     = diag(L); 
    [blah,ii] = sort(diagL); 
    diffy     = freqq_shift(16)-freqq_shift(10);
    [szm,szn] = size(ii);
    if diffy < 0
      if (szm > szn)
        ii=flipud(ii);
       else
        ii=fliplr(ii);
        end
       end
    if ((band == 740) & (prb ~= 'R'))
      ii=1:length(ii);
      end
    diagL     = diagL(ii); 
    L         = diag(diagL);
    B         = A;
    A         = zeros(size(A));

    for jj = 1:length(ii)
      A(:,jj) = B(:,ii(jj));
      end
    %clf; mesh(abs(inv(A)*H*A)); title('after reodering A'); pause(1)

    if plotbwah > 0
      baba = [jq freqq_shift abs(diag(L)) ... 
          freqq_shift./abs(diagLorig) freqq_shift./abs(diag(L))];
      fprintf(1,'%3i %12.6f %12.6f %8.6f %8.6f \n',baba')
      fprintf(1,'       \n');
      subplot(311); plot(1:length(freqq_shift),freqq_shift,'+')
      title(num2str(band));
      subplot(312); plot(1:length(freqq_shift)-1,diff(freqq_shift))
      subplot(313); plot(1:length(freqq_shift),freqq_shift./abs(diag(L)),...
                       1:length(freqq_shift),freqq_shift./abs(diagLorig))
      pause
      end

    arg1 = trans_ampl'*A;

    %divide populations by freqq(i) and then multiply abs.coef. by f to get
    %rid of stimulated emission.
    population_tp = population_t./freqq_shift;
    arg2 = inv(A)*diag(population_tp)*trans_ampl;

    w_tot = (pressure_self*w_selfq+pressure_for*w_forq)/pressure_ref; 

    diagL = diag(L);

    if plotbwah > 0
      clf; 
      l1 = 1:length(freqq_shift);
      %orig        semilogy(freqq_shift,abs(arg1.*arg2'),freqq_shift,strenqt);
      subplot(211);semilogy(abs(diagL),abs(arg1.*arg2'),'+',...
                          freqq_shift,strenqt,'+');
      subplot(212);semilogy(l1,abs(arg1.*arg2'),l1,strenqt);
      title(num2str(band)); pause
      end

    %this is a call to a Mex file
    %here, birn = sum(full)* sum(lor*birn)/sum(lor)
    %fullmix4 =dofull(f,freqq_shift,w_tot,temperature,duration,pressure_for,...
    %            pressure_self,arg1,arg2,diagL,K_scale_mixing,bsm,birnbirn);

    %this is a call to a Mex file
    %here, birn  =  sum(full*birn) which is what Dave Tobin did in his thesis
    fullmix4 = dofullNEWBIRN(f,freqq_shift,w_tot,temperature,duration,...
      pressure_for,pressure_self,arg1,arg2,diagL,K_scale_mixing,bsm,birnbirn);
  end

else
  disp('ratio < 1e-2')
  %be more careful .. use only 2/3 strongest lines for the mixing calcs
  [gg,hh]  =  sort(strenqt); 
  maxx = floor(0.33*length(strenqt));
  indx = hh(maxx:length(strenqt));
  indx1 = hh(1:maxx-1);

  %be more careful .. use only 2/3 lowest j lines for the mixing calcs
  [gg,hh] = sort(jq); maxx = floor(0.99*length(jq));

  indx = hh(1:maxx-1);
  indx1 = hh(maxx:length(jq));

  H = diag(freqq_shift(indx)) + sqrt(-1)*W_plus(indx,indx);
  [A,L] = eig(H);

  if bwah > 0
    [blah,ii]    = sort(freqq_shift);
    if ((band == 740) & (prb ~= 'R'))
      ii = 1:length(ii);
      end
    freqq_shift  = freqq_shift(ii);
    w_selfq      = w_selfq(ii);
    w_forq       = w_forq(ii);
    trans_ampl   = trans_ampl(ii);
    population_t = population_t(ii);
    strenqt      = strenqt(ii);
    jq           = jq(ii);

    diagLorig = diag(L);
    diagL = diag(L); [blah,ii] = sort(diagL); 
    diffy     = freqq_shift(16)-freqq_shift(10);
    [szm,szn] = size(ii);
    if diffy < 0
      if (szm > szn)
        ii = flipud(ii);
       else
        ii = fliplr(ii);
        end
       end
    if ((band == 740) & (prb ~= 'R'))
      ii = 1:length(ii);
      end
    diagL = diagL(ii); 
    L = diag(diagL);
    B = A;
    for jj = 1:length(ii)
      A(:,jj) = B(:,ii(jj));
      end

    if plotbwah > 0
      baba = [jq(indx) freqq_shift(indx) abs(diag(L)) ... 
          freqq_shift(indx)./abs(diagLorig) freqq_shift(indx)./abs(diag(L))];
      fprintf(1,'%3i %12.6f %12.6f %8.6f %8.6f \n',baba')
      fprintf(1,'       \n');
      subplot(311); plot(1:length(freqq_shift(indx)),freqq_shift(indx),'+')
      title(num2str(band));
      subplot(312); 
      plot(1:length(freqq_shift(indx))-1,diff(freqq_shift(indx)),'+')
      subplot(313); 
      plot(1:length(freqq_shift(indx)),freqq_shift(indx)./abs(diag(L)),...
            1:length(freqq_shift(indx)),freqq_shift(indx)./abs(diagLorig))
      pause
      end
    end

  arg1 = trans_ampl(indx)'*A;

  %divide populations by freqq(i) and then multiply abs.coef. by f to get
  %rid of stimulated emission.
  population_tp(indx) = population_t(indx)./freqq(indx);
  %% save /home/sergio/KCARTA/SRC/NONLTE2/sergio/NONLTE_MATLAB/a.mat A
  arg2 = inv(A)*diag(population_tp(indx))*trans_ampl(indx);

  w_tot = ...
     (pressure_self*w_selfq(indx)+pressure_for*w_forq(indx))/pressure_ref; 

  diagL = diag(L);

  if (plotbwah > 0)
    clf; semilogy(freqq_shift(indx),abs(arg1.*arg2'),...
                 freqq_shift(indx),strenqt(indx)); 
    title(num2str(band)); pause
    end

  %this is a call to a Mex file
  %here, birn = sum(full)* sum(lor*birn)/sum(lor)
  %fullmix4=dofull(f,freqq_shift(indx),w_tot(indx),...
  %             temperature,duration,pressure_for,pressure_self,...
  %             arg1,arg2,diagL,K_scale_mixing,bsm,birnbirn);

  %this is a call to a Mex file
  %here, birn = sum(full*birn) which is what Dave Tobin did in his thesis
  fullmix4a = dofullNEWBIRN(f,freqq_shift(indx),w_tot(indx),...
                          temperature,duration,pressure_for,pressure_self,...
                          arg1,arg2,diagL,K_scale_mixing,bsm,birnbirn);

  %this is first order line mixing
  %voi1a=voigtmix2(freqq_shift(indx),f,ymix1(indx),jq(indx),...  
  %             temperature,w_forq,w_selfq,strenqt,stuff,layeramt,'V',birn);   

  %this is pure lorentz
  %ymix=zeros(size(freqq_shift(indx)));
  %voia=voigtmixRatio(freqq_shift(indx),f,ymix,jq(indx),temperature,...
  %   w_forq(indx),w_selfq(indx),strenqt(indx),stuff,layeramt,'V',birn,1.0);


  %now do the weaker lines
  ymix = zeros(size(freqq_shift(indx1)));
  nomix = voigtmixRatio(freqq_shift(indx1),f,ymix,jq(indx1),temperature,...
     w_forq(indx1),w_selfq(indx1),strenqt(indx1),stuff,layeramt,'V',birn,1.0);

%  subplot(211); plot(freqq,jq,'+'); title('jq'); grid; 
%  subplot(212); semilogy(freqq,strenqt,'+'); title('jq'); grid; pause
%  clf
%  plot(f,exp(-(fullmix4a)),f,exp(-(nomix))); title('T1'); grid; pause
%  plot(f,exp(-(voi1a)),f,exp(-(nomix))); title('T2'); grid; pause
%  plot(f,exp(-(voia)),f,exp(-(nomix))); title('T3'); grid; pause
%  plot(f,fullmix4a,f,voi1a,f,voia,f,nomix); title('k'); grid; pause
%  plot(f,fullmix4a./voia,f,voi1a./voia); title('k/klor'); grid; pause

  b = find(fullmix4a<0);
  fullmix4a(b) = 0.0;

  fullmix4 = fullmix4a+nomix;
  end

freqq        = freqq0;
jq           = jq0;
w_selfq      = w_selfq0;
w_forq       = w_forq0;
trans_ampl   = trans_ampl0;
population_t = population_t0;
strenqt      = strenqt0;
ymix1        = ymix10;

%plot(f,exp(-fullmix4)); pause(1)
%plot(f,fullmix4); pause

