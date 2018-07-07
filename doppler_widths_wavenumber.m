k       = 1.380658e-23; 
c_light = 2.99792458e8;         %ms-1 
amu     = 1.6605402e-27;        %nucleon mass/kg 

v0 = 0.1:0.1:40000;   %% new upto VIS/UV
v0 = 605 : 1 : 2830;  %% IR
v0 = 200 : 1 : 3200;  %% FIR-NIR

m  = 48;
m  = 44;
mass    = m*amu;
T  = [100 150 200 250 300 350]; 
for ii = 1 : length(T)
  y(ii,:) =  v0.*sqrt(2*log(2)*k*T(ii)/mass/c_light/c_light);
end

figure(1); plot(v0,y);        xlabel('v0 cm-1'); ylabel('\delta \nu (doppler)');
  hl = legend(num2str(T'),'location','northwest'); set(hl,'fontsize',10); grid

figure(2); plot(v0,y/0.0025); xlabel('v0 cm-1'); ylabel('\delta \nu (doppler)/0.0025');
  hl = legend(num2str(T'),'location','northwest'); set(hl,'fontsize',10); grid

