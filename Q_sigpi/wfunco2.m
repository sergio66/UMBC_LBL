function dif=wfun(xx,elower,wq_tape,jall,temperature,stuff)

global quiet 

B0=stuff.B0;
btz=stuff.btz;
beta=stuff.beta;
band=stuff.band;

a1=xx(1);a2=xx(2);a3=xx(3);

no_lines=length(elower);	%no_lines=length(jall);

energy_diff=ones(no_lines,1)*[elower]'-[elower]*ones(1,no_lines);
energy_diff=abs(energy_diff);

K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature);
J=ones(no_lines,1)*(2*jall+1)';    % This is NOT J, it is 2J+1
K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature).*J'./J;

%%%%%%%%% Calculating the widths: (i.e. the diagonal elements of W) %%%%%%%%
%
% This calculation is taken from D.P. Edwards and L.L.Strow  "Spectral Line
% Shape Considerations for Limb Temperature Sounders" JGR. Vol 96 1991.
%

%    index   1   2   3   4   5    
%         _______________________
%  K=   1 |  k00 k10 k20 k30 k40 |  0
%       2 |  k01 k11 k21 k31 k41 |  1
%       3 |  k02 k12 k22 k32 k42 |  2
%       4 |  k03 k13 k23 k33 k43 |  3
%       5 |  k04 k14 k24 k34 k44 |  4
%         ________________________
%             0   1   2   3   4     J

if (band ~= 662)
  i_even=find(rem(jall,2)==0);      %   i=1,3,5,7 ...   J=0,2,4,6...
  i_odd=find(rem(jall,2)~=0);       %   i=2,4,6,8 ...   J=1,2,5,7...
  i_even2=i_even(2:length(i_even)); %   i=  3,5,7 ...   J=  2,4,6...

  width_even_lower=-0.5*sum(K(i_even,i_even2));           % k20+k22+k24+...
  width_even_upper=-0.5*beta*sum(K(i_even2,i_even2));     % k22+k24+k26+...
  width_odd_upper=-0.5*(1-beta)*sum(K(i_odd,i_even2));    % k21+k23+k25+...

  widthq=(width_even_lower+width_even_upper+width_odd_upper)';

  %%%%%%% because of the 662 isotope breaking symmetry, we do this
  %%%%%%% ie instead of only Q2,4,6,... being allowed, we can have Q1,2,3,4..
%  minn=min([length(widthq) length(wq_tape)]);
%  widthq=widthq(1:minn);
%  wq_tape=wq_tape(1:minn);

elseif ((band == 662) & (mod(stuff.evenodd,2) == 0))  %even jq's
  i_even=find(rem(jall,2)==0);      %   i=1,3,5,7 ...   J=0,2,4,6...
  i_odd=find(rem(jall,2)~=0);       %   i=2,4,6,8 ...   J=1,3,5,7...
  i_even2=i_even(2:length(i_even)); %   i=  3,5,7 ...   J=  2,4,6...

  width_even_lower=-0.5*beta*sum(K(i_even2,i_even2));      % k20+k22+k24+...
  width_odd_lower=-0.5*(1-beta)*sum(K(i_odd,i_even2));    % k21+k23+k25+...
  width_even_upper=-0.5*beta*sum(K(i_even,i_even2));     % k22+k24+k26+...
  width_odd_upper=-0.5*(1-beta)*sum(K(i_odd,i_even2));    % k21+k23+k25+...

  widthq=(width_even_lower+width_even_upper+width_odd_upper)';

  %%%%%%% because of the 662 isotope breaking symmetry, we do this
  %%%%%%% ie instead of only Q2,4,6,... being allowed, we can have Q1,2,3,4..
elseif ((band == 662) & (mod(stuff.evenodd,2) == 1))  %odd jq's
  i_even=find(rem(jall,2)==0);      %   i=1,3,5,7 ...   J=0,2,4,6...
  i_odd=find(rem(jall,2)~=0);       %   i=2,4,6,8 ...   J=1,3,5,7...
  i_even2=i_even(2:length(i_even)); %   i=  3,5,7 ...   J=  2,4,6...

  width_even_lower=-0.5*(1-beta)*sum(K(i_even2,i_odd));      % k20+k22+k24+...
  width_odd_lower=-0.5*beta*sum(K(i_odd,i_odd));    % k21+k23+k25+...
  width_even_upper=-0.5*(1-beta)*sum(K(i_even,i_odd));     % k22+k24+k26+...
  width_odd_upper=-0.5*beta*sum(K(i_odd,i_odd));    % k21+k23+k25+...

  widthq=(width_even_lower+width_even_upper+width_odd_upper)';

  %%%%%%% because of the 662 isotope breaking symmetry, we do this
  %%%%%%% ie instead of only Q2,4,6,... being allowed, we can have Q1,2,3,4..
  end

% Calculated the difference between fitted widths and data widths%%%%%%%%%%
% The square of the difference will be minimized via the least squares
% optimization method.

dif=(widthq-wq_tape)';   
if (band == 662)
  jqplot=(2:2:no_lines);
  ind=find(jqplot < 76);
  dif=dif(ind);

  %plot the results
  %clf;
  %jqplot=(2:2:no_lines);
  %jqplot=jqplot(ind);
  %subplot(211);plot(jqplot,wq_tape(ind),'*',jqplot,widthq(ind))
  %ylabel('Width');title('Fitting for a1,a2,a3 via widths');
  %subplot(212);         plot(jqplot,dif)
  %xlabel('f');ylabel('Diff.');pause(1)

  end

%lq=1:length(wq_tape);
%plot(lq,widthq,lq,wq_tape)
%pause(0.1)

% wq_tape are the actual air broadened widths from the
                   % HITRAN tape (which have been temperature corrected)
                   % dif will be minimized to get the best values for a1,a2,a3

pflag=quiet;
if pflag > 0
  %plot the results
  clf;
  jqplot=(2:2:no_lines);
  subplot(211);plot(jqplot,wq_tape,'*',jqplot,widthq)
  ylabel('Width');title('Fitting for a1,a2,a3 via widths');
  subplot(212);         plot(jqplot,dif)
  xlabel('f');ylabel('Diff.');pause(1)
end

