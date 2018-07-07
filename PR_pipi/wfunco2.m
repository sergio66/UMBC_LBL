function dif = wfunco2(xx,elower,wq_tape,jall,temperature,stuff)

global quiet

band = stuff.band;

prb = stuff.prb;
B0 = stuff.B0;
btz = stuff.btz;
beta = stuff.beta;

% Variables are:
%  	 xx:		starting values for a1,a2,a3
%
%    elower:		HITRAN energies
%
%     wq_tape:		wq are the temperature adjusted widths
%			for the R-branch lines. (Only even J's)
%
% temperature:		temperature

no_lines = length(elower);  

%%%%%%%%%%%%%%%%%%%%%% Energy difference calculation %%%%%%%%%%%%%%%%%%%%%%%
% Here, the energy difference between levels is calculated

energy_diff = ones(no_lines,1)*[elower]'-[elower]*ones(1,no_lines);
energy_diff = abs(energy_diff);

%%%%%%%%%%%%%%%%%%%%%% Calculation of collision rates %%%%%%%%%%%%%%%%%%%%
a1 = xx(1); a2 = xx(2); a3 = xx(3);
K = a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature);
  % Detailed balance is achieved by using the equation:
  %   Kjj' *(2j'+1)*exp(-Ej'/kT) = Kj'j *(2j+1)*exp(-Ej/kT)
J = ones(no_lines,1)*(2*jall+1)';    % This is NOT J, it is 2J+1
K = tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature).*J'./J;

%%%%%%%%% Calculating the widths: (i.e. the diagonal elements of W) %%%%%%%%
%
% This calculation is taken from D.P. Edwards and L.L.Strow  "Spectral Line
% Shape Considerations for Limb Temperature Sounders" JGR. Vol 96 1991.
%
%    index    1    2    3    4
%         ______________________
% K=   1   | k11  k21  k31  k41 | 1 
%      2   | k12  k22  k32  k42 | 2
%      3   | k13  k23  k33  k43 | 3
%      4   | k14  k24  k34  k44 | 4
%         ______________________
%             1    2    3    4    J                    

start2vs3 = 2;
if ((band == 2322))
  start2vs3 = 3;
  %now put zeros
  aa = 2:2:length(jall); bb = 1:2:length(jall);
  K(aa,bb) = 0.0; K(bb,aa) = 0.0;
  end

i_even  = find(rem(jall,2) == 0);   % i=2,4,6...   J=2,4,6...
i_odd   = find(rem(jall,2) ~= 0);   % i=1,3,5...   J=1,3,5...
i_even1 = i_even-1;              

width_even_lower(i_even) = -0.5*beta*sum(K(i_even,i_even));    % k22+k24+k26...
width_even_lower(i_odd)  = -0.5*(1-beta)*sum(K(i_even,i_odd)); % k12+k14+k16...
width_odd_lower(i_even)  = -0.5*(1-beta)*sum(K(i_odd,i_even)); % K21+k23+k25...
width_odd_lower(i_odd)   = -0.5*beta*sum(K(i_odd,i_odd));      % k11+k13+k15...

width_even_upper(i_even) = -0.5*beta*sum(K(i_even,i_even));    % k22_k24+k26...
width_even_upper(i_odd)  = -0.5*(1-beta)*sum(K(i_even,i_odd)); % k12+k14+k16...
width_odd_upper(i_even)  = -0.5*(1-beta)*sum(K(i_odd,i_even)); % k21+k23+k25...
width_odd_upper(i_odd)   = -0.5*beta*sum(K(i_odd,i_odd));      % k11+k13+k15...

if prb == 'P'
  width_even_lower = width_even_lower(start2vs3:length(width_even_lower));
  width_even_upper = width_even_upper(start2vs3:length(width_even_upper));
  width_odd_lower  = width_odd_lower(start2vs3:length(width_odd_lower));
  width_odd_upper  = width_odd_upper(start2vs3:length(width_odd_upper));
  end

widthq = width_even_upper+width_even_lower+width_odd_upper+width_odd_lower;

%%% Calculated the difference between fitted widths and data widths%%%%%%%%%%
% The square of the difference will be minimized via the least squares
% optimization method.
%

dif = (widthq'-wq_tape)';  % wq are the actual air broadened widths from the
                   % HITRAN tape (which have been temperature corrected)
                   % dif will be minimized to get the best values for a1,a2,a3

pflag = quiet;
if pflag > 0
   %plot the results
   clf;
  jqplot = (1:max(jall));
  if prb == 'P'
    jqplot = (2:max(jall));
    end
  subplot(211);plot(1:length(wq_tape),wq_tape,'*',1:length(widthq),widthq)
  ylabel('Width');title('Fitting for a1,a2,a3 via widths');
  subplot(212);plot(dif)
  xlabel('f');ylabel('Diff.'); 
  pause(0.1)
  end
