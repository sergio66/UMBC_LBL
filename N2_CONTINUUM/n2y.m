ipfile = 'IPFILES/n2one';

%% see Lafferty 1996 paper  vol 35, no 30 : applied optics pg 5911
MGC = 8.314674269981136  ;
prof = load(ipfile);
pp = prof(:,3); T = prof(:,4); q = prof(:,5);
p = prof(:,2);

L = q * 1000 * 10000;   %% kmoles/cm2 --> moles/m2
L = L * MGC * T ./(pp * 101325)  %% gas cell length in m

am = p * 101325 * 273/(T*101325);  %% gas density in amagats

xx = 2300 + [0:10:140];
B   = [0.111 0.139 0.176 0.197 0.175 0.142 0.132 ...
       0.129 0.130 0.133 0.135 0.131 0.124 0.116 0.104]*1e-5;  
beta =[-247 -150 -160 -181 -195 -205 -211 -210 ...
       -205 -190 -168 -143 -108 -63  +1];;
abscoeff = (p * 273/T)^2*(0.8387 - 0.0754*T/296);
abscoeff = (0.8387 - 0.0754*T/296);
abscoeff = (p * 273/T)^1 * (0.8387 - 0.0754*T/296);
abscoeff = abscoeff .* B .* exp(beta*(1/296 - 1/T));    %% in am-1 cm-1 

odx = abscoeff * (L*100) * am;  %% multiply abscoeff by length(in cm) * den(am)

plot(xx,odx,'r')