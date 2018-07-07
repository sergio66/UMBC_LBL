function df=wfun(xx,elower,wr_tape,jall,temperature,stuff) 

prb=stuff.prb; 
B0=stuff.B0; 
btz=stuff.btz; 
beta=stuff.beta;

a1=xx(1);a2=xx(2);a3=xx(3);
no_lines=length(elower);	

energy_diff=ones(no_lines,1)*[elower]'-[elower]*ones(1,no_lines);
energy_diff=abs(energy_diff)+eye(no_lines,no_lines);
  
K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature).*...
  (1-eye(no_lines,no_lines));
J=ones(no_lines,1)*(2.*jall+1)'; % This is NOT J, it is 2J+1
K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature).*J'./J;

i_even=find(rem(jall,2)==0);   % i=2,4,6...   J=2,4,6...
i_odd=find(rem(jall,2)~=0);    % i=1,3,5...   J=1,3,5...
i_even1=i_even-1;
i_odd2=i_even;
i_odd3=i_odd(2:length(i_odd)); % i=3,5,7...   J=3,5,7...






width_even_lower(i_even)=-0.5*beta*sum(K(i_even,i_even));   
width_even_lower(i_odd)=-0.5*(1-beta)*sum(K(i_even,i_odd)); 
width_odd_lower(i_even)=-0.5*(1-beta)*sum(K(i_odd,i_even));
width_odd_lower(i_odd)=-0.5*beta*sum(K(i_odd,i_odd));     

width_even_upper(i_even)=-0.5*beta*sum(K(i_even,i_even));   
width_even_upper(i_odd)=-0.5*(1-beta)*sum(K(i_even,i_odd));
width_odd_upper(i_even)=-0.5*(1-beta)*sum(K(i_odd,i_even)); 
width_odd_upper(i_odd)=-0.5*beta*sum(K(i_odd,i_odd));      

if prb=='P'
  width_even_lower=width_even_lower(2:length(width_even_lower));
  width_even_upper=width_even_upper(2:length(width_even_upper));
  width_odd_lower=width_odd_lower(2:length(width_odd_lower));
  width_odd_upper=width_odd_upper(2:length(width_odd_upper));
  end

df1=(width_even_lower+width_odd_lower+width_even_upper+width_odd_upper)/a1;
df1=df1(1:length(wr_tape));




K1=-K.*log(energy_diff/B0);
width_even_lower(i_even)=-0.5*beta*sum(K1(i_even,i_even));   
width_even_lower(i_odd)=-0.5*(1-beta)*sum(K1(i_even,i_odd)); 
width_odd_lower(i_even)=-0.5*(1-beta)*sum(K1(i_odd,i_even));
width_odd_lower(i_odd)=-0.5*beta*sum(K1(i_odd,i_odd));     

width_even_upper(i_even)=-0.5*beta*sum(K1(i_even,i_even));   
width_even_upper(i_odd)=-0.5*(1-beta)*sum(K1(i_even,i_odd));
width_odd_upper(i_even)=-0.5*(1-beta)*sum(K1(i_odd,i_even)); 
width_odd_upper(i_odd)=-0.5*beta*sum(K1(i_odd,i_odd));      

if prb=='P'
  width_even_lower=width_even_lower(2:length(width_even_lower));
  width_even_upper=width_even_upper(2:length(width_even_upper));
  width_odd_lower=width_odd_lower(2:length(width_odd_lower));
  width_odd_upper=width_odd_upper(2:length(width_odd_upper));
  end

df2=(width_even_lower+width_odd_lower+width_even_upper+width_odd_upper)/2*2;
df2=df2(1:length(wr_tape));





K2=-K*btz.*energy_diff/temperature;
width_even_lower(i_even)=-0.5*beta*sum(K2(i_even,i_even));   
width_even_lower(i_odd)=-0.5*(1-beta)*sum(K2(i_even,i_odd)); 
width_odd_lower(i_even)=-0.5*(1-beta)*sum(K2(i_odd,i_even));
width_odd_lower(i_odd)=-0.5*beta*sum(K2(i_odd,i_odd));     

width_even_upper(i_even)=-0.5*beta*sum(K2(i_even,i_even));   
width_even_upper(i_odd)=-0.5*(1-beta)*sum(K2(i_even,i_odd));
width_odd_upper(i_even)=-0.5*(1-beta)*sum(K2(i_odd,i_even)); 
width_odd_upper(i_odd)=-0.5*beta*sum(K2(i_odd,i_odd));      

if prb=='P'
  width_even_lower=width_even_lower(2:length(width_even_lower));
  width_even_upper=width_even_upper(2:length(width_even_upper));
  width_odd_lower=width_odd_lower(2:length(width_odd_lower));
  width_odd_upper=width_odd_upper(2:length(width_odd_upper));
  end

df3=(width_even_lower+width_odd_lower+width_even_upper+width_odd_upper)/2*2;
df3=df3(1:length(wr_tape));





df=[df1' df2' df3']';
