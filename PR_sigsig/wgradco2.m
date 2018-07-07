function df=wfun(xx,elower,wr_tape,jall,temperature,stuff) 
 
prb=stuff.prb; 
B0=stuff.B0; 
btz=stuff.btz; 
band=stuff.band;

a1=xx(1);a2=xx(2);a3=xx(3);
no_lines=length(elower);

% Here, the energy difference between levels is calculated


energy_diff=ones(no_lines,1)*[elower]'-[elower]*ones(1,no_lines);
energy_diff=abs(energy_diff)+eye(no_lines,no_lines);
  
K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature).*...
  (1-eye(no_lines,no_lines));
J=ones(no_lines,1)*(2*jall+1)'; % This is NOT J, it is 2J+1
K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature).*J'./J;

if (band ~= 2352)
  i_even=find(rem(jall,2)==0);             % i=1,3,5...   J=0,2,4...
  i_odd=find(rem(jall,2)~=0);              % i=2,4,6...   J=1,3,5...

  width_even_lower=-0.5*sum(K(i_even,i_even));    % k00+k02+k04... 
  width_odd_upper=-0.5*sum(K(i_odd,i_odd));       % k11+k13+k15... 
  if prb=='P'
    width_odd_upper=width_odd_upper(1:length(width_odd_upper)-1);
    end
  df1=(width_even_lower + width_odd_upper)/a1;

  K1=-K.*log(energy_diff/B0);
  width_even_lower=-0.5*sum(K1(i_even,i_even));    % k00+k02+k04... 
  width_odd_upper=-0.5*sum(K1(i_odd,i_odd));       % k11+k13+k15... 
  if prb=='P'
    width_odd_upper=width_odd_upper(1:length(width_odd_upper)-1);
    end
  df2=(width_even_lower + width_odd_upper)/2*2;

  K2=-K*btz.*energy_diff/temperature;
  width_even_lower=-0.5*sum(K2(i_even,i_even));    % k00+k02+k04... 
  width_odd_upper=-0.5*sum(K2(i_odd,i_odd));       % k11+k13+k15... 
  if prb=='P'
    width_odd_upper=width_odd_upper(1:length(width_odd_upper)-1);
    end
  df3=(width_even_lower + width_odd_upper)/2*2;

  df=[df1' df2' df3']';

else
  i_even=find(rem(jall,2)==0);             % i=1,3,5...   J=0,2,4...
  i_odd=find(rem(jall,2)~=0) ;             % i=2,4,6...   J=1,3,5...
  i_even2=i_even(2:length(i_even)); %   i=  3,5,7 ...   J=  2,4,6... 

  width_even_lower(i_even)=-0.5*sum(K(i_even,i_even));    % k00+k02+k04...
  width_odd_lower(i_odd)=-0.5*sum(K(i_odd,i_odd));    % k00+k02+k04...
  width_even_upper(i_even)=-0.5*sum(K(i_even,i_even));    % k00+k02+k04...
  width_odd_upper(i_odd)=-0.5*sum(K(i_odd,i_odd));    % k00+k02+k04...

  minn=min([length(width_even_lower) length(width_odd_lower)]);
  if (prb == 'R')
    width_even_lower=width_even_lower(1:minn);
    width_even_upper=width_even_upper(1:minn);
    width_odd_lower=width_odd_lower(1:minn);
    width_odd_upper=width_odd_upper(1:minn);
  else
    %had 2:minn
    width_even_lower=width_even_lower(1:minn-1);
    width_even_upper=width_even_upper(1:minn-1);
    width_odd_lower=width_odd_lower(1:minn-1);
    width_odd_upper=width_odd_upper(1:minn-1);
    end
  df1=(width_even_lower+width_odd_upper+width_even_upper+width_odd_lower)/a1;


  K1=-K.*log(energy_diff/B0);
  width_even_lower(i_even)=-0.5*sum(K1(i_even,i_even));    % k00+k02+k04...
  width_odd_lower(i_odd)=-0.5*sum(K1(i_odd,i_odd));    % k00+k02+k04...
  width_even_upper(i_even)=-0.5*sum(K1(i_even,i_even));    % k00+k02+k04...
  width_odd_upper(i_odd)=-0.5*sum(K1(i_odd,i_odd));    % k00+k02+k04...

  minn=min([length(width_even_lower) length(width_odd_lower)]);
  if (prb == 'R')
    width_even_lower=width_even_lower(1:minn);
    width_even_upper=width_even_upper(1:minn);
    width_odd_lower=width_odd_lower(1:minn);
    width_odd_upper=width_odd_upper(1:minn);
  else
    width_even_lower=width_even_lower(1:minn-1);
    width_even_upper=width_even_upper(1:minn-1);
    width_odd_lower=width_odd_lower(1:minn-1);
    width_odd_upper=width_odd_upper(1:minn-1);
    end
  df2=(width_even_lower+width_odd_upper+width_even_upper+width_odd_lower)/2*2;


  K2=-K*btz.*energy_diff/temperature;
  width_even_lower(i_even)=-0.5*sum(K2(i_even,i_even));    % k00+k02+k04...
  width_odd_lower(i_odd)=-0.5*sum(K2(i_odd,i_odd));    % k00+k02+k04...
  width_even_upper(i_even)=-0.5*sum(K2(i_even,i_even));    % k00+k02+k04...
  width_odd_upper(i_odd)=-0.5*sum(K2(i_odd,i_odd));    % k00+k02+k04...

  minn=min([length(width_even_lower) length(width_odd_lower)]);
  if (prb == 'R')
    width_even_lower=width_even_lower(1:minn);
    width_even_upper=width_even_upper(1:minn);
    width_odd_lower=width_odd_lower(1:minn);
    width_odd_upper=width_odd_upper(1:minn);
  else
    width_even_lower=width_even_lower(1:minn-1);
    width_even_upper=width_even_upper(1:minn-1);
    width_odd_lower=width_odd_lower(1:minn-1);
    width_odd_upper=width_odd_upper(1:minn-1);
    end
  df3=(width_even_lower+width_odd_upper+width_even_upper+width_odd_lower)/2*2;

  df=[df1' df2' df3']';
  end
