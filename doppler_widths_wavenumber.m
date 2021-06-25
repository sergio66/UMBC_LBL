function y = doppler_widths_wavenumber(v0,T,m)

k       = 1.380658e-23; 
c_light = 2.99792458e8;         %ms-1 
amu     = 1.6605402e-27;        %nucleon mass/kg 

if nargin == 0
  v0 = 0.1:0.1:40000;   %% new upto VIS/UV
  v0 = 605 : 1 : 2830;  %% IR
  v0 = 200 : 1 : 3200;  %% FIR-NIR
end

if nargin <= 1
  T  = [100 150 200 250 300 350];
  T  = [100 : 10 : 350];
  T  = [200 : 20 : 320];
end

if nargin <= 2
  m  = 18;  %% WV
  m  = 48;  %% O3
  m  = 44;  %% Co2
end

mass    = m*amu;

for ii = 1 : length(T)
  y(ii,:) =  v0.*sqrt(2*log(2)*k*T(ii)/mass/c_light/c_light);
end

figure(1); plot(v0,y);        xlabel('v0 cm-1'); ylabel('\delta \nu (doppler)');
  hl = legend(num2str(T'),'location','northwest'); set(hl,'fontsize',10); grid
  xlim([min(v0) max(v0)]);

figure(2); plot(v0,y/0.0025); xlabel('v0 cm-1'); ylabel('\delta \nu (doppler)/0.0025');
  hl = legend(num2str(T'),'location','northwest'); set(hl,'fontsize',10); grid
  xlim([min(v0) max(v0)]);

if m == 44
  boo = find(v0 >= 670,1); 
  junk = [min(y(:,boo)) max(y(:,boo)) mean(y(:,boo)) std(y(:,boo))];
  fprintf(1,'CO2 667 cm-1   min/max mean/std dv =  %8.6f %8.6f %8.6f %8.6f cm-1 \n',junk);

  boo = find(v0 >= 2300,1); 
  junk = [min(y(:,boo)) max(y(:,boo)) mean(y(:,boo)) std(y(:,boo))];
  fprintf(1,'CO2 2300 cm-1   min/max mean/std dv =  %8.6f %8.6f %8.6f %8.6f cm-1 \n',junk);

elseif m == 18
  boo = find(v0 >= 1500,1); 
  junk = [min(y(:,boo)) max(y(:,boo)) mean(y(:,boo)) std(y(:,boo))];
  fprintf(1,'WV 1500 cm-1   min/max mean/std dv =  %8.6f %8.6f %8.6f %8.6f cm-1 \n',junk);

elseif m == 48
  boo = find(v0 >= 1050,1); 
  junk = [min(y(:,boo)) max(y(:,boo)) mean(y(:,boo)) std(y(:,boo))];
  fprintf(1,'O3 1050 cm-1   min/max mean/std dv =  %8.6f %8.6f %8.6f %8.6f cm-1 \n',junk);
end

