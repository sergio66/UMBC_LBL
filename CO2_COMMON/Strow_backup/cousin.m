function [chi]=cousin(f,freqq,w_tot,temperature,pressure_for,pressure_self)

%%%%%%%%%%%%%%%%%%%%%%%%%%% cousin.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This file computes Cousin Chi function abs. coeff.
%
% from /salsify/scratch4/Strow/Tobin_home/tobin/Co2pr/Cousin
%%%%%%%%%%Chi Function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   N2 broadened
chi_n2_238 = x_n2_238(f,freqq);
if temperature >= 238
  chi_n2_296 = x_n2_296(f,freqq);
  chi_n2=chi_n2_238+(chi_n2_296-chi_n2_238)/(296-238)*(temperature-238);
else
  chi_n2_193 = x_n2_193(f,freqq);
  chi_n2=chi_n2_193+(chi_n2_238-chi_n2_193)/(238-193)*(temperature-193);
  end

%    CO2 broadened
chi_co2_218 = x_co2_218(f,freqq);
chi_co2_296 = x_co2_296(f,freqq);
chi_co2=chi_co2_218+(chi_co2_296-chi_co2_218)/(296-218)*(temperature-218);

chi=(pressure_self*chi_co2+pressure_for*chi_n2)/(pressure_self+pressure_for);

plot(f,chi); pause(0.1);


