function [qfcn]=q(A,B,C,D,E,lines,T); 
% function [qfcn]=q(A,B,C,D,E,lines,T); 
% initialize coefficients vectors for qtips coefficients 

a1 = ones(length(lines.iso),1);b1 = ones(length(lines.iso),1); 
c1 = ones(length(lines.iso),1);d1 = ones(length(lines.iso),1); 

% Assign coefficients according to isotope 
no_isotopes = max(lines.iso); 
for i = 1: no_isotopes 
  ind = find(lines.iso == i); 
  if (length(ind) > 0)
    a1(ind) = a1(ind)*A(i); 
    b1(ind) = b1(ind)*B(i); 
    c1(ind) = c1(ind)*C(i); 
    d1(ind) = d1(ind)*D(i); 
    end
end 

% Evaluate partition functions at desired temperature and 296K 
Qt   = a1 + b1*T   + c1*T^2   + d1*T^3; 
Q296 = a1 + b1*296.0 + c1*(296.0^2) + d1*(296.0^3); 

qfcn = Q296./Qt;
