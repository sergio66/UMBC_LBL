%%%%%%%%%%%%%%%%%%%%%%%%%%% klormix_birn.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file computes:
%		(1) Lorentz with 1st order mixing absorption coefficient
%			with and without Birnbaum Chi factor
%		(2) Lorentz absorption coefficient
%			with and without Birnbaum Chi factor
%	It also checks SUM{Si*Yi} to see if detailed balance worked.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing absorption coeffs.')
freqr_shift=freqr+frequency_shift/100;
freqp_shift=freqp+frequency_shift/100;

for iii=1:length(pick_no)
  ind=find(f(:,iii)~=0);
  fnow=f(ind,iii);trd=t_rawdata(ind,iii);kback=K_back(ind,iii);

  K_scale_lor=density*pressure_self(iii)/pressure_ref*temperature_ref*...
	path_length(iii)/temperature(iii)/pi;

  lor=zeros(length(fnow),1);		  % lorentz
  lor_birn=zeros(length(fnow),1);	  % lorentz & birnbaum
  lormix=zeros(length(fnow),1);  	  % lorentz & 1st order mixing
  lormix_birn=zeros(length(fnow),1);   % lorentz & 1st order mixing & birnbaum

%  for i=1:length(freqr_shift)
  for i=1:-1
    w_tot=(pressure_self*w_selfr(i)+pressure_for*w_forr(i))/pressure_ref;
   temp=strenrt(i)*(w_tot)*(fnow/freqr_shift(i))./((fnow-freqr_shift(i)).^2+...
        (w_tot).^2);
    temp2=strenrt(i)*(fnow/freqr_shift(i)).*(w_tot+...
       (fnow-freqr_shift(i))*ymixr(i))./((fnow-freqr_shift(i)).^2+(w_tot).^2);
    lor=lor+temp;
    lormix=lormix+temp2;

    Am=1;ii=sqrt(-1);tau0=0.72/temperature;tau2=duration;
    dnu=fnow-freqr_shift(i);zz=sqrt((w_tot^2 + dnu.^2).*(tau0^2 + tau2^2));
    ex=exp(tau2*w_tot + tau0*dnu);
%   chi=Am.*zz.*(-pi/2).*besselh_tobin(1,ii*zz).*ex;
    chi=Am.*zz.*(-pi/2).*real(besselh(1,ii*zz)).*ex;
    lor_birn=lor_birn+temp.*chi;lormix_birn=lormix_birn+temp2.*chi;

    [freqr(i) dnu(1) dnu(length(dnu)) chi(1) chi(length(chi)) ... 
        zz(1) zz(length(zz)) ex(1) ex(length(1)) strenrt(i) ymixr(i) w_tot ...
        temp2(1) temp2(length(temp2))]
    figure(1)
    subplot(211); plot(fnow,temp2.*chi.*K_scale_lor);
    subplot(212): plot(fnow,chi); 
    figure(2); plot(fnow,lormix_birn*K_scale_lor); pause(0.1);
    end

  for i=1:length(freqp_shift)
    w_tot=(pressure_self*w_selfp(i)+pressure_for*w_forp(i))/pressure_ref;
  temp=strenpt(i)*(w_tot)*(fnow/freqp_shift(i))./((fnow-freqp_shift(i)).^2 +...
	(w_tot).^2);
    temp2=strenpt(i)*(fnow/freqp_shift(i)).*(w_tot+...
      (fnow-freqp_shift(i))*ymixp(i))./((fnow-freqp_shift(i)).^2+(w_tot).^2);
    lor=lor+temp;lormix=lormix+temp2;

    Am=1;ii=sqrt(-1);tau0=0.72/temperature;tau2=duration;
    dnu=fnow-freqp_shift(i);zz=sqrt((w_tot^2 + dnu.^2).*(tau0^2 + tau2^2));
    ex=exp(tau2*w_tot + tau0*dnu);
%	chi=Am.*zz.*(-pi/2).*besselh_tobin(1,ii*zz).*ex;
   chi=Am.*zz.*(-pi/2).*real(besselh(1,ii*zz)).*ex;
   lor_birn=lor_birn+temp.*chi;lormix_birn=lormix_birn+temp2.*chi;

    figure(1)
    [freqp(i) dnu(1) dnu(length(dnu)) chi(1) chi(length(chi)) ... 
        zz(1) zz(length(zz)) ex(1) ex(length(ex)) strenpt(i) ymixp(i) w_tot ...
        temp2(1) temp2(length(temp2))]
    subplot(211); plot(fnow,temp2.*chi.*K_scale_lor);
    subplot(212): plot(fnow,chi); 
    figure(2); plot(fnow,lormix_birn*K_scale_lor);pause(0.1);
   end

  lor=lor*K_scale_lor*bsm+kback*bsm;
  lor_birn=lor_birn*K_scale_lor*bsm+kback*bsm;
  lormix=lormix*K_scale_lor*bsm+kback*bsm;
  lormix_birn=lormix_birn*K_scale_lor*bsm+kback*bsm;

  pflag=0;
  if pflag==1
    clg;subplot(2,1,1),plot(f,t_rawdata,'r',f,exp(-lormix_birn),'b')
    ylabel('T');title('Minimizing Obs-Calc(1st order mixing & Birnbaum)');
    subplot(2,1,2),plot(f,t_rawdata-exp(-lormix_birn),'b')
    ylabel('Diff.');xlabel('cm-1');drawnow
    end

  end

tau0
tau2
