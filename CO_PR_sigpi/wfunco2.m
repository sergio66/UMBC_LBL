function dif=wfun(xx,elower,wq_tape,jall,temperature,stuff)

global quiet

prb=stuff.prb;
B0=stuff.B0;
btz=stuff.btz;
band=stuff.band;

a1=xx(1);a2=xx(2);a3=xx(3);

no_lines=length(elower);        %no_lines=length(jall); 

energy_diff=ones(no_lines,1)*[elower]'-[elower]*ones(1,no_lines);
energy_diff=abs(energy_diff);

K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature);
J=ones(no_lines,1)*(2*jall+1)';        % This is NOT J, it is 2J+1
J=ones(length(jall),1)*(2*jall+1)';    % This is NOT J, it is 2J+1

no_lines
whos K energy_diff  J jall

K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature).*J'./J;

%%%%%%%%% Calculating the widths: (i.e. the diagonal elements of W) %%%%%%%%
%
%          index    1   2   3   4
%   K=           ___________________
%            1   | k00 k10 k20 k30 | 0
%	     2	 | k01 k11 k21 k31 | 1
%	     3	 | k02 k12 k22 k32 | 2
%	     4	 | k03 k13 k23 k33 | 3
%		 ___________________
%                   0   1   2   3    J


if (band ~= 2352)        %all these are isotope 1 or 2
  i_even=find(rem(jall,2)==0);             % i=1,3,5...   J=0,2,4...
  i_odd=find(rem(jall,2)~=0);              % i=2,4,6...   J=1,3,5...

  width_even_lower=-0.5*sum(K(i_even,i_even));    % k00+k02+k04...
  width_odd_upper=-0.5*sum(K(i_odd,i_odd));       % k11+k13+k15...

  if prb=='P'
    width_odd_upper=width_odd_upper(1:length(width_odd_upper)-1);
    end
  widthq=width_even_lower + width_odd_upper;
  dif=(widthq'-wq_tape)';
  minn=min([length(widthq) length(wq_tape)]);

else     %this is isotope 3 ===> all  R=0,1,2,...  P=1,2,3.. allowed
%  minn=min([length(widthq) length(wq_tape)]);
%  widthq=widthq(1:minn);
%  wq_tape=wq_tape(1:minn);

  i_even=find(rem(jall,2)==0);             % i=1,3,5...   J=0,2,4...
  i_odd=find(rem(jall,2)~=0);              % i=2,4,6...   J=1,3,5...
  i_even2=i_even(2:length(i_even)); %   i=  3,5,7 ...   J=  2,4,6... 

  width_even_lower(i_even)=-0.5*sum(K(i_even,i_even));    % k00+k02+k04...
  width_odd_lower(i_odd)=-0.5*sum(K(i_odd,i_odd));    % k00+k02+k04...
  width_even_upper(i_even)=-0.5*sum(K(i_even,i_even));    % k00+k02+k04...
  width_odd_upper(i_odd)=-0.5*sum(K(i_odd,i_odd));    % k00+k02+k04...

  minn=min([length(width_even_lower) length(width_odd_lower)]);
  if (prb == 'R')
    %have R0,1,2,3,4,...
    widthq=(width_even_lower(1:minn)+width_even_upper(1:minn)+...
          width_odd_lower(1:minn)++width_odd_upper(1:minn)); 
  else
    %only have P1,2,3,4,...
    %had 2:minn
    widthq=(width_even_lower(1:minn-1)+width_even_upper(1:minn-1)+...
          width_odd_lower(1:minn-1)++width_odd_upper(1:minn-1)); 
    end

  minn=min([length(widthq) length(wq_tape)]);

  widthq=widthq(1:minn);
  wq_tape=wq_tape(1:minn);
  dif=(widthq'-wq_tape)';
  end

%%% Calculated the difference between fitted widths and data widths%%%%%%%%%%


if (band ~= 2352)
  pflag=quiet;
  if pflag > 0
    clf;
    jrplot=(0:2:max(jall));
    if prb=='P'
      jrplot=(2:2:max(jall));
      end
    subplot(2,1,1);plot(jrplot,wq_tape(1:minn),'*',jrplot,widthq(1:minn)')
    ylabel('Width');title('Fitting for a1,a2,a3 via widths');
    subplot(2,1,2);plot(jrplot,dif)
    xlabel('f');ylabel('Diff.');
    pause(0.2)
    end  

else
  pflag=quiet;
  if pflag>0
    clf;
    jrplot=(0:max(jall));
    if prb=='P'
      jrplot=(1:max(jall));
      end
    jrplot=jrplot(1:minn);
    subplot(2,1,1);plot(jrplot,wq_tape(1:minn),'*',jrplot,widthq(1:minn)')
    ylabel('Width');title('Fitting for a1,a2,a3 via widths');
    subplot(2,1,2);plot(jrplot,dif)
    xlabel('f');ylabel('Diff.');
    pause(0.2)
    end
  end
