function [fullmix4]=full_mix4(freqq,f,W_plus,jq,w_selfq,w_forq,...
                              temperature,trans_ampl,population_t,stuff,birn) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  full_mix4.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file computes the full mixing absorption coefficient
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
density=stuff.density; 
pressure_self=stuff.pressure_self; 
pressure_for=stuff.pressure_for; 
pressure_ref=stuff.pressure_ref; 
temperature_ref=stuff.temperature_ref; 
path_length=stuff.path_length; 
bsm=stuff.bsm; 
duration=stuff.duration;
frequency_shift=stuff.frequency_shift;

K_scale_mixing=density*pressure_self/pressure_ref*temperature_ref*...
	path_length/temperature/pi;
no_lines=length(freqq); no_pts=length(f);
k=zeros(1,no_pts);

freqq_shift=freqq+frequency_shift/100;

H=diag(freqq_shift)+i*W_plus;
[A,L]=eig(H);
arg1=trans_ampl'*A;

%divide populations by freqq(i) and then multiply abs.coef. by f to get
%rid of stimulated emission.
population_tp=population_t./freqq;
arg2=inv(A)*diag(population_tp)*trans_ampl;

rd=population_t.*trans_ampl;
len=length(rd);

%sss=W_plus.*(ones(len,1)*rd').*(trans_ampl*ones(1,len));
%sss=triu(sss,1);
%denom=diag(W_plus).*rd.*trans_ampl;
%denom=sum(denom);
%sss=1+2*sum(sum(sss))/denom

w_tot=(pressure_self*w_selfq+pressure_for*w_forq)/pressure_ref; 

if ((birn == 'b') | (birn == 'B') )        %birnbaum
  sumlorbirn=0;
  sumlor=0;
  end

for i=1:no_lines
  temp=arg1(i)*arg2(i)./(f-L(i,i));

  if ((birn == 'b') | (birn == 'B') )        %birnbaum
    %this is the broadening
    chi=birnbaum(f,freqq_shift(i),w_tot(i),temperature,duration);
    thelor=lorentz(f,freqq(i),temperature,1.00,w_tot(i));
    sumlor=sumlor+thelor;
    sumlorbirn=sumlorbirn+thelor.*chi;
  elseif ((birn == 'c') | (birn == 'C') )    %cousin
    error('cannot have FULL mixing and cousin!!');
    %this is the broadening
    chi=cousin1(f,freqq_shift(i),w_tot(i),temperature,...
               pressure_for,pressure_self);    
  else
    chi=ones(size(temp));
    end

  if ((birn == 'c') | (birn == 'C') )    %cousin
    temp=temp.*chi;
    end
  k=k+temp;
  end

if ((birn == 'b') | (birn == 'B') )    
  fullmix4=K_scale_mixing*f.*imag(k)*bsm.*(sumlorbirn./sumlor);
else             %we have already done cousin or no chi
  fullmix4=K_scale_mixing*f.*imag(k)*bsm;
  end
