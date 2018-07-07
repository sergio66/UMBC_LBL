press = load('IPFILES/oco_p');
vmr = load('IPFILES/oco_vmr_o2');
temperature = load('IPFILES/oco_t');
z = load('IPFILES/oco_zhgt');

%% pV = nRT ==> q = L n/V = L p/RT
MGC = 8.314674269981136  ;

press = press /101325;     %% need atm
partpress = press .* vmr;  %% need atm
z = z * 100;               %% need cm

GasAmt = 101325 * z.* partpress/1e9/MGC./temperature; %change to kmoles/cm2

xstartup

std_o2 = load('IPFILES/std_o2');
figure(1); semilogy(std_o2(:,4),std_o2(:,2),temperature,press,'r'); set(gca,'ydir','reverse')
figure(2); semilogy(std_o2(:,5),std_o2(:,2),GasAmt,press,'r');      set(gca,'ydir','reverse')
figure(3); semilogy(p2h(std_o2(:,2)*1013.25)*100,std_o2(:,2),z,press,'r');  
  set(gca,'ydir','reverse')

fid = fopen('IPFILES/oco_o2_profile','w');
data = [(1:length(press))' press partpress temperature GasAmt];
whos data
fprintf(fid,'%3i %8.6e %8.6e  %8.6f  %8.6e \n',data');
fclose(fid);

c = 3e8;
lambda = [396750000000000 388500000000000];  %% start stop 
vu = c./lambda;
vu = vu *1e6;
vu = 10000./vu;

%%[w,d] = run8(7,14000,15000,'IPFILES/n2one');  
%% plot(w,d)
%% [w,d] = run8(7,14000,15000,'IPFILES/oco_o2_profile');

topts.ffin = 0.01/5;
%[w,d] = run8(7,14200,14800,'IPFILES/co2two',topts);  
%[wx,dx] = run8(7,14250,14750,'IPFILES/oco_o2_profile',topts);
%save OCO_PWang/oco_o2_vis.mat wx dx

%[w,d] = run8(7,12950,13250,'IPFILES/co2two',topts);  
[wx2,dx2] = run8(7,12950,13250,'IPFILES/oco_o2_profile',topts);
save OCO_PWang/oco_o2_vis2.mat wx2 dx2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error('lll')

clear all
cd OCO_PWang
abs_read

load oco_o2_vis2

plot(10000./Wvnm,sum(TAU_TG'),10000./wx2,sum(dx2),'r'); 
set(gca,'xdir','reverse')

plot(10000./Wvnm,exp(-sum(TAU_TG')),...
     10000./wx2,exp(-sum(dx2)),'r'); set(gca,'xdir','reverse')

plot(10000./Wvnm,exp(-sum(TAU_TG')),...
     10000./wx2,exp(-sum(dx2)*0.075),'r'); set(gca,'xdir','reverse')
