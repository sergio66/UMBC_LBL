%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wgradco2.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This function simply speeds up the optimization done in wfunco2.m by
%    including the gradient in the calculation.
%
function df=wgradco2(xx,elower,wq,jall,prb,num)

global beta beta_pure beta_for bsm duration frequency_shift
global pressure_self pressure_for btz B0 temperature temperature_ref; 
global density Boltzmann mass_CO2 speed_light path_length pressure_ref;
global K_scale_mixing K_scale_voigt K_scale_lor;
global trans_ampl population population_t t_rawdata;
global voi_back
global K_back strenrt w_forr w_selfr

a1=xx(1);a2=xx(2);a3=xx(3);
B0=0.4;  btz=1.4387863;
no_lines=length(elower);

% Here, the energy difference between levels is calculated

energy_diff=ones(no_lines,1)*[elower]'-[elower]*ones(1,no_lines);
energy_diff=abs(energy_diff)+eye(no_lines,no_lines);
  
K=a1*((energy_diff/B0).^(-a2)).*exp(-a3*btz*energy_diff/temperature(num)).*...
  (1-eye(no_lines,no_lines));
J=ones(no_lines,1)*(2*jall+1)'; % This is NOT J, it is 2J+1
K=tril(K,-1)+triu(K,1).*exp(btz*energy_diff/temperature(num)).*J'./J;

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

K2=-K*btz.*energy_diff/temperature(num);
        width_even_lower=-0.5*sum(K2(i_even,i_even));    % k00+k02+k04...
	width_odd_upper=-0.5*sum(K2(i_odd,i_odd));       % k11+k13+k15...
if prb=='P'
	width_odd_upper=width_odd_upper(1:length(width_odd_upper)-1);
end

df3=(width_even_lower + width_odd_upper)/2*2;


df=[df1' df2' df3']';
