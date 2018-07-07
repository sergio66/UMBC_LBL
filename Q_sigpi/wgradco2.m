function df = wgradco2(xx,elower,wq,jall,temperature,stuff)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wgradco2.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This function simply speeds up the optimization done in wfunco2.m by
%    including the gradient in the calculation.
%

B0=stuff.B0;  
btz=stuff.btz;  
beta=stuff.beta;
band=stuff.band;

a1=xx(1);a2=xx(2);a3=xx(3);
no_lines=length(elower); % jall goes from 1 to 51 & elower is for j=1 to 51

% Here, the energy difference between levels is calculated

energy_diff=ones(no_lines,1)*[elower]'-[elower]*ones(1,no_lines);
energy_diff=abs(energy_diff)+eye(no_lines,no_lines);
  
K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature).*...
  (1-eye(no_lines,no_lines));
J=ones(no_lines,1)*(2*jall+1)';     % This is NOT J, it is 2J+1
K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature).*J'./J;

if (band ~= 662)
  i_even=find(rem(jall,2)==0);
  i_odd=find(rem(jall,2)~=0);
  i_even2=i_even(2:length(i_even));

  width_even_lower=-0.5*sum(K(i_even,i_even2));
  width_even_upper=-0.5*beta*sum(K(i_even2,i_even2));
  width_odd_upper=-0.5*(1-beta)*sum(K(i_odd,i_even2));

  df1=(width_even_lower + width_even_upper + width_odd_upper)/a1;

  K1=-K.*log(energy_diff/B0);
  width_even_lower=-0.5*sum(K1(i_even,i_even2));
  width_even_upper=-0.5*beta*sum(K1(i_even2,i_even2));
  width_odd_upper=-0.5*(1-beta)*sum(K1(i_odd,i_even2));

  df2=(width_even_lower + width_even_upper + width_odd_upper)/2*2;

  K2=-K*btz.*energy_diff/temperature;
  width_even_lower=-0.5*sum(K2(i_even,i_even2));
  width_even_upper=-0.5*beta*sum(K2(i_even2,i_even2));
  width_odd_upper=-0.5*(1-beta)*sum(K2(i_odd,i_even2));

  df3=(width_even_upper + width_even_lower + width_odd_upper)/2*2;

  df=[df1' df2' df3']';

elseif ((band == 662) & (mod(stuff.evenodd,2) == 0))  %even jq's
  i_even=find(rem(jall,2)==0);      %   i=1,3,5,7 ...   J=0,2,4,6...
  i_odd=find(rem(jall,2)~=0);       %   i=2,4,6,8 ...   J=1,3,5,7...
  i_even2=i_even(2:length(i_even)); %   i=  3,5,7 ...   J=  2,4,6...

  width_even_lower=-0.5*beta*sum(K(i_even2,i_even2));      % k20+k22+k24+...
  width_odd_lower=-0.5*(1-beta)*sum(K(i_odd,i_even2));    % k21+k23+k25+...
  width_even_upper=-0.5*beta*sum(K(i_even,i_even2));     % k22+k24+k26+...
  width_odd_upper=-0.5*(1-beta)*sum(K(i_odd,i_even2));    % k21+k23+k25+...
  df1=(width_even_lower+width_even_upper+width_odd_upper+width_odd_lower)/a1;

  K1=-K.*log(energy_diff/B0);
  width_even_lower=-0.5*beta*sum(K1(i_even2,i_even2));      % k20+k22+k24+...
  width_odd_lower=-0.5*(1-beta)*sum(K1(i_odd,i_even2));    % k21+k23+k25+...
  width_even_upper=-0.5*beta*sum(K1(i_even,i_even2));     % k22+k24+k26+...
  width_odd_upper=-0.5*(1-beta)*sum(K1(i_odd,i_even2));    % k21+k23+k25+...
  df2=(width_even_lower+width_even_upper+width_odd_upper+width_odd_lower);

  K2=-K*btz.*energy_diff/temperature;
  width_even_lower=-0.5*beta*sum(K2(i_even2,i_even2));      % k20+k22+k24+...
  width_odd_lower=-0.5*(1-beta)*sum(K2(i_odd,i_even2));    % k21+k23+k25+...
  width_even_upper=-0.5*beta*sum(K2(i_even,i_even2));     % k22+k24+k26+...
  width_odd_upper=-0.5*(1-beta)*sum(K2(i_odd,i_even2));    % k21+k23+k25+...
  df3=(width_even_lower+width_even_upper+width_odd_upper+width_odd_lower);

  df=[df1' df2' df3']';

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
  df1=(width_even_lower+width_even_upper+width_odd_upper+width_odd_lower)/a1;

  K1=-K.*log(energy_diff/B0);
  width_even_lower=-0.5*(1-beta)*sum(K1(i_even2,i_odd));      % k20+k22+k24+...
  width_odd_lower=-0.5*beta*sum(K1(i_odd,i_odd));    % k21+k23+k25+...
  width_even_upper=-0.5*(1-beta)*sum(K1(i_even,i_odd));     % k22+k24+k26+...
  width_odd_upper=-0.5*beta*sum(K1(i_odd,i_odd));    % k21+k23+k25+...
  df2=(width_even_lower+width_even_upper+width_odd_upper+width_odd_lower);

  K2=-K*btz.*energy_diff/temperature;
  width_even_lower=-0.5*(1-beta)*sum(K2(i_even2,i_odd));      % k20+k22+k24+...
  width_odd_lower=-0.5*beta*sum(K2(i_odd,i_odd));    % k21+k23+k25+...
  width_even_upper=-0.5*(1-beta)*sum(K2(i_even,i_odd));     % k22+k24+k26+...
  width_odd_upper=-0.5*beta*sum(K2(i_odd,i_odd));    % k21+k23+k25+...
  df3=(width_even_lower+width_even_upper+width_odd_upper+width_odd_lower);

  df=[df1' df2' df3']';

  %%%%%%% because of the 662 isotope breaking symmetry, we do this
  %%%%%%% ie instead of only Q2,4,6,... being allowed, we can have Q1,2,3,4..
  end


if (band==662)
  jqplot=(2:2:no_lines);
  ind=find(jqplot < 76);
  df=[df1(ind)'   df2(ind)'  df3(ind)']';
  end


