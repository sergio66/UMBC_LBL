function y = doppler_lorentz_widths_wavenumber_usstd(m)

%% can easily edit to read in trp_ mls_ mlw_ sas_ saw_ std_water

if m == 18
  gasID = 1;
elseif m == 44
  gasID = 2;
elseif m == 48
  gasID = 3;
end
HITRAN        = '/asl/rta/hitran/h16.by.gas';
fnamePRE = [HITRAN '/g'];
fnamePOST = '.dat';
fnameIN   = int2str(gasID);
fname     = [fnamePRE fnameIN fnamePOST];
[lineORIG,hitran_version,hlist_qtips] = hitread(500,3000,0,gasID,fname,-1);

profile = load('/home/sergio/SPECTRA/IPFILES/std_water');
        press      = profile(:,2); %% atm
        partpress  = profile(:,3); %% atm
        tempr      = profile(:,4); %% K

for jj = 1 : 100
        pwr       = lineORIG.abcoef;
        for_brd   =  lineORIG.abroad;
        self_brd  = lineORIG.sbroad;
	freq      = lineORIG.wnum+press(jj)*lineORIG.tsp;  %freq shift
        energy    = lineORIG.els;
	s0        = lineORIG.stren;
        brd(:,jj) = broad(press(jj),partpress(jj),1.0,for_brd,self_brd,pwr,tempr(jj),gasID);
end

v0 = 500 : 25 : 3000;
yd = doppler_widths_wavenumber(v0,tempr,m);

figure(2)
i700 = find(v0 >= 700,1);
semilogx(max(brd),1:100,yd(:,i700),1:100,'linewidth',2); hl = legend('Lorenz broad','Doppler broad','location','best','fontsize',10);
  xlabel('Broadening cm-1'); ylabel('layer num');

loglog(max(brd),press*1013,yd(:,i700),press*1013,'linewidth',2); 
  xlabel('Broadening cm-1'); ylabel('p(mb)');
  set(gca,'ydir','reverse'); ylim([0.01 1000])
  hl = legend('Lorenz broad','Doppler broad','location','best','fontsize',10);

semilogx(max(brd),p2h(press*1013)/1000,yd(:,i700),p2h(press*1013)/1000,'linewidth',2); 
  xlabel('Broadening cm-1'); ylabel('H(km)');
  %set(gca,'ydir','reverse'); ylim([0.01 1000])
  hl = legend('Lorenz broad','Doppler broad','location','best','fontsize',10);

figure(3);
  dn = 0:0.002:0.2; plot(dn,hist(lineORIG.abroad,dn),'b',dn,hist(lineORIG.sbroad,dn),'r'); 
  ylabel('histogram'); title('(b)abroad (r)sbroad'); xlabel('\gamma cm-1/atm');
