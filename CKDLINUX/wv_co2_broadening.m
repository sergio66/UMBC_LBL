function y = wv_co2_broadening(lines);

% Fourier transform infrared spectroscopy measurements of H2O-broadened half-widths of CO2 at 4.3 mm1
% Keeyoon Sung, Linda R. Brown, Robert A. Toth, and Timothy J. Crawford
% Can. J. Phys. 87: 469–484 (2009) doi:10.1139/P08-130

p(1) = 0.137;                 a = -0.005017;     a1 = -0.03563;
p(2) = 3.0e-2;                b = 1.085;         b1 = 0.06318;
p(3) = 2.63e-3;               c = 35.67;         d1 = 0.001498;
p(4) = -8.17e-6;              q = 48.84;         e1 = -1.848e-5;
p(5) = 1.11e-7;                                  f1 = 4.924e-8;
p(6) = -5.51e-8;                                 c1 = 0.1093;

%{
m1 = 0:61;
F = zeros(size(m1));
for ii = 1 : 6
  F = F + p(ii) * (abs(m1).^(ii-1));
end

m2 = 62:121;
G = (a*m2.^2 + b*abs(m2) + c) ./ (abs(m2) + q);

m3 = 0:121;
H = a1./m3./m3 + b1./abs(m3) + d1*abs(m3) + e1*m3.^2 + f1*(abs(m3)).^3 + c1;

figure(1); plot(m1,F,'b',m2,G,'g',m3,H,'r')
figure(2); plot(m3,H,'r')
y = H;
%}

if ~isfield(lines,'jlower')
  disp('warning : no jlower : setting it using bslq')
  if ~isfield(lines,'bslq')
    disp('warning : no bslq .. .just using 0.13 cm-1 atm-')
    lines.jlower = ones(size(lines.stren));
  else
    lines.jlower = str2num(lines.bslq(:,6:8));
  end
end  
m3 = lines.jlower;
y = a1./m3./m3 + b1./abs(m3) + d1*abs(m3) + e1*m3.^2 + f1*(abs(m3)).^3 + c1;
bad = find(isnan(y) | isinf(y));
y(bad) = 0.137;
%figure(1); plot(m3,y,'.');