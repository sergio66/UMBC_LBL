ipfile = 'IPFILES/n2one';

f1 = 2305; f2 = 2430;

topts.O2O3N2continuumVers = 4;
[w,dx4] = run8(22,f1,f2,ipfile,topts);      %% new lbl
disp('done 4')

topts.O2O3N2continuumVers = 1;
[w,dx1] = run8(22,f1,f2,ipfile,topts);      %% new lbl
disp('done 1')

topts.O2O3N2continuumVers = 3;
[w,dx3] = run8(22,f1,f2,ipfile,topts);      %% new lbl
disp('done 3')

topts.O2O3N2continuumVers = -1;
[w,dNOCON] = run8(22,f1,f2,ipfile,topts);      %% new lbl
disp('done NOCON')


%[w,d] = run8(22,2305,2430,ipfile);      %% new lbl
%[w,dorig] = run8(22,2305,2430,ipfile);  %% old lbl
%[w,d] = run8(22,2305,2310,ipfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see Lafferty 1996 paper  vol 35, no 30 : applied optics pg 5911
MGC = 8.314674269981136  ;
prof = load(ipfile);
pp = prof(:,3); T = prof(:,4); q = prof(:,5);
p = prof(:,2);

L = q * 1000 * 10000;   %% kmoles/cm2 --> moles/m2
L = L * MGC * T ./(pp * 101325)  %% gas cell length in m

am = p * 101325 * 273/(T*101325);  %% gas density in amagats

x = 2300 + [0:10:140];
B   = [0.111 0.139 0.176 0.197 0.175 0.142 0.132 ...
       0.129 0.130 0.133 0.135 0.131 0.124 0.116 0.104]*1e-5;  
beta =[-247 -150 -160 -181 -195 -205 -211 -210 ...
       -205 -190 -168 -143 -108 -63  +1];;
%abscoeff = (p * 273/T)^2*(0.8387 - 0.0754*T/296);
abscoeff = (p*273/T) * (0.8387 - 0.0754*T/296);
abscoeff = abscoeff .* B .* exp(beta*(1/296 - 1/T));    %% in am-1 cm-1 

od = abscoeff * (L*100) * am;  %% multiply abscoeff by length(in cm) * den(am)

plot(x,od)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is for a path length of 1 m
genln2 = load('/asl/packages/Genln2/Glab/lite.tau');
genln2 = load('genln2_n2.tau');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dUMBC = dx3;  dORIG = dx1; dNEW = dx4;

odx = interp1(x,od,w);
clf
plot(w,dNEW,w,dORIG,'o-',w,dUMBC,...
     w,odx+dNOCON,'mo-',genln2(:,1),genln2(:,2)*L,'k','Linewidth',2);  
hl = ...
  legend('new LBLRTM','ORIG UMBC LBL','CURRENT UMBC LBL','Lafferty','Genln2');
set(hl,'fontsize',10); ylabel('OD for 235m'); xlabel('cm-1')
grid 

%% so ORIG and GENLN2 are the same ....


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
