%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  tdif_birn.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file computes the difference between 1st order mixing transmission
%       with birnbaum and actual transmission.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function final_diff=tdif_birn(bstart,freqr,freqp,jp,f,flag4)

global beta beta_pure beta_for bsm duration frequency_shift
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor
global trans_amplr populationr population_tr t_rawdata frequency_shift
global trans_amplp populationp population_tp strenpt w_forp w_selfp
global voi_back pick_no
global K_back strenrt w_forr w_selfr abscoef 
global sigsig_comp sigsigR sigsigP ymixR wwwR

load strenths
[beta duration]

for i=1:length(temperature)
  eval(['global W_plus_r' num2str(i)])
  eval(['global W_plus_p' num2str(i)])
  end

i=0;
if flag4(1)==1
  i=i+1;bsm=bstart(i);
  end

if flag4(2)==1
  i=i+1;beta=bstart(i);
  end
if flag4(3)==1
  i=i+1;duration=bstart(i);
  end
if flag4(4)==1
  i=i+1;frequency_shift=bstart(i);
  end

band_strength_multiplier=bsm;

if (beta > 0)
  beta=min([ beta  1]);
else
  beta=0.95;
  end

if (duration < 0)
  duration=5.0e-3;
  end
if (duration > 1e-2)
  duration=5e-3;
  end

diff=[];

for iii=1:length(temperature)
  ind=find(f(:,iii)~=0);
  fnow=f(ind,iii);
  trd=t_rawdata(ind,iii);
  kback=K_back(ind,iii);

  K_scale_lor=density*pressure_self(iii)/pressure_ref*...
		temperature_ref*path_length(iii)/temperature(iii)/pi;
  abscoefP=zeros(length(fnow),1);  
  abscoefR=zeros(length(fnow),1);  
  abscoef=zeros(length(fnow),1);  

  %%%%%%%%%%%%% R-branch lines %%%%%%%%%%%%%%%%%%%%%%%%%%%
%  fprintf(1,'the pressures are \n');
%  [pressure_self(iii) pressure_for(iii)]
%  fprintf(1,'the strengths are \n');
%  strenrt(:,iii)
%  fprintf(1,'the w forr are \n');
%  w_forr(:,iii)
%  fprintf(1,'the w selfr are \n');
%  w_selfr(:,iii)

  freqr_shift=freqr+frequency_shift/100;
  no_lines=length(freqr_shift);Y_1stmix=zeros(1,no_lines);
  ZZ=(beta+(1-beta)*eye(no_lines,no_lines));
  eval(['W_plus2=W_plus_r' num2str(iii) '.*ZZ;']);clear ZZ
  for n=1:no_lines
    for m=1:no_lines
      if m~=n
        Y_1stmix(n)=Y_1stmix(n)+trans_amplr(m,iii)*(W_plus2(m,n))/...
            (freqr_shift(n)-freqr_shift(m));
        end
      end
    end
  ymix=2*Y_1stmix./trans_amplr(:,iii)';
  ymixR=ymix;
  for i=1:length(freqr_shift)
    w_tot=(pressure_self(iii)*w_selfr(i,iii)+pressure_for(iii)*...
		w_forr(i,iii))/pressure_ref;
    wwwR(i)=w_tot;
    temp=strenrt(i,iii)*(fnow/freqr_shift(i)).*...
	    (w_tot+(fnow-freqr_shift(i))*ymix(i))./...
	    ((fnow-freqr_shift(i)).^2+(w_tot).^2);

% want to call the function as 
    Am=1;ii=sqrt(-1);tau0=0.72/temperature(iii);tau2=duration;
    dnu=fnow-freqr_shift(i);
    zz=sqrt((w_tot^2 + dnu.^2).*(tau0^2 + tau2^2));
    ex=exp(tau2*w_tot + tau0*dnu);
    bes=real(besselh(1,ii*zz));
    chi=Am.*zz.*(-pi/2).*bes.*ex;

    figure(2); plot(fnow,chi)
    temp=temp.*chi;
    abscoefR=abscoefR+temp;
    end

  abscoefR=abscoefR*K_scale_lor;
  sigsigR=abscoefR;

  %%%%%%%%%%%%% P-branch lines %%%%%%%%%%%%%%%%%%%%%%%%%%%
  freqp_shift=freqp+frequency_shift/100;no_lines=length(freqp_shift);  
  Y_1stmix=zeros(1,no_lines);
  ZZ=(beta+(1-beta)*eye(no_lines,no_lines));
  eval(['W_plus2=W_plus_p' num2str(iii) '.*ZZ;']);clear ZZ
  for n=1:no_lines
    for m=1:no_lines
      if m~=n
        Y_1stmix(n)=Y_1stmix(n)+trans_amplp(m,iii)*(W_plus2(m,n))/...
	       (freqp_shift(n)-freqp_shift(m));
        end
      end
    end
  ymix=2*Y_1stmix./trans_amplp(:,iii)';
  for i=1:length(freqp_shift)
    w_tot=(pressure_self(iii)*w_selfp(i,iii)+pressure_for(iii)*...
		w_forp(i,iii))/pressure_ref;
    temp=strenpt(i,iii)*(fnow/freqp_shift(i)).*(w_tot+...
		(fnow-freqp_shift(i))*ymix(i))./...
		((fnow-freqp_shift(i)).^2+(w_tot).^2);

    Am=1;ii=sqrt(-1);tau0=0.72/temperature(iii);tau2=duration;
    dnu=fnow-freqp_shift(i);
    zz=sqrt((w_tot^2 + dnu.^2).*(tau0^2 + tau2^2));
    ex=exp(tau2*w_tot + tau0*dnu);
    bes=real(besselh(1,ii*zz));
    chi=Am.*zz.*(-pi/2).*bes.*ex;

    temp=temp.*chi;
    abscoefP=abscoefP+temp;
    end

  abscoefP=abscoefP*K_scale_lor;

  lorcousin=abscoefP;
  ind=find((fnow > 0) & (fnow > max(freqp_shift)+40)); 
  ymix0=zeros(size(freqp_shift));  
  voivoi=vcousin(iii,freqp,fnow',ymix0,temperature(iii),w_forp(:,iii),... 
            w_selfp(:,iii),strenpt(:,iii),'V','c',1.0); 
  lorcousin(ind)=voivoi(ind);
  figure(2); plot(fnow,exp(-lorcousin),fnow,exp(-abscoefP)); pause(0.1)
  figure(1);
  abscoefP=lorcousin;

  sigsigP=abscoefP;

  %%%%%%%%%%%%% add P and R branch lines %%%%%%%%%%%%%%%%%%%%%%%%%%%

  abscoef=abscoefP+abscoefR;
  sigsig_comp=abscoef;

  abscoef=abscoef*bsm+kback*bsm;
  error=trd-exp(-abscoef);

  %figure(iii);clg;
  %plot(fnow,trd,'r',fnow,exp(-abscoef),'g',fnow,error,'b');grid;drawnow
	
  subplot(length(temperature),1,iii)
  plot(fnow,trd,'r',fnow,exp(-abscoef),'g',fnow,error,'b');grid;drawnow
  axis([min(fnow) max(fnow) -0.2 1.2])
  fss=frequency_shift;
%eval(['save fileGL' num2str(pick_no(iii)) ' fnow trd abscoef tau2 beta frequency_shift'])
%eval(['save fileIND' num2str(pick_no(iii)) ' fnow trd abscoef tau2 beta frequency_shift'])
eval(['save file' num2str(pick_no(iii)) ' fnow trd abscoef tau2 beta fss kback sigsigR sigsigP'])

  diff=[diff;error];
  end

nbeta=num2str(beta);nbsm=num2str(bsm);
ndur=num2str(duration);nfs=num2str(frequency_shift);
disp(['bsm=' nbsm '   beta=' nbeta '   tau2=' ndur '   fs=' nfs])
fprintf(1,'beta= %10.8f dofc = %10.8f \n',beta,duration);

final_diff=diff;

