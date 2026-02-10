fname = [hitranpath '/h16.by.gas/'];  %% this is technically same as do_HITRAN_vers.m
fname = [hitranpath '/h20.by.gas/'];  %% this is technically same as do_HITRAN_vers.m
fname = [hitranpath '/h24.by.gas/'];  %% this is technically same as do_HITRAN_vers.m

fop = '/carrot/s1/sergio/pin_feb2002_sea_airsnadir_op.22deg.rtp';
[h,ha,p,pa] = rtpread(fop);
for ii = 1 : 30
  str = ['doink = p.gas_' num2str(ii,'%2d') '(4,49);']; eval(str);
  gas(ii) = doink;  %% in molecules/cm2 at surface
  end
lays = 1:p.nlevs(49)-1;
loglog(p.gas_1(lays,49),p.plays(lays,49),...
       p.gas_2(lays,49),p.plays(lays,49),...
       p.gas_3(lays,49),p.plays(lays,49),...
       p.gas_4(lays,49),p.plays(lays,49),...
       p.gas_5(lays,49),p.plays(lays,49),...
       p.gas_6(lays,49),p.plays(lays,49),...
       p.gas_9(lays,49),p.plays(lays,49),'Linewidth',2)
hold on
loglog(p.gas_12(lays,49),p.plays(lays,49),'b--','Linewidth',2)
set(gca,'ydir','reverse'); grid on
xlabel('amount (molecules/cm2)'); ylabel('p(mb)'); title('US Std Gas Profiles')
h1 = ...
legend('Water','CO2','O3','N2O','CO','CH4','HNO3','SO2','Location','NorthEast')
set(h1,'Fontsize',8)
axis([1e9 1e22 0 1000])

iWhich = input(' +1 +2 +3 +4 +5 +6 : MW IR NIR VIS (FIR-NIR) (IR-VIS) : ');

if iWhich == 6   %% all
  for gid = 1 : 30
    [lines,hitran_version] = hitread_simple(00001,30000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(10000./lines.wnum,lines.stren); title(num2str(gid));
      axis([10000/30000 10000/01 1e-30 (max(lines.stren))]);
      xlabel('Wavelength(um)'); ylabel('log10(linestren)')
      pause
      end
    end

elseif iWhich == 5   %% fir-nir
  for gid = 1 : 30
    [lines,hitran_version] = hitread_simple(00300,10000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(10000./lines.wnum,lines.stren); title(num2str(gid));
      axis([10000/10000 10000/00300 1e-30 (max(lines.stren))]);
      xlabel('Wavelength(um)'); ylabel('log10(linestren)')
      printer = ['print -dpng /home/sergio/PAPERS/SPECIALTALKS/IOWA08/FIGS/'];
      printer = [printer 'gas_stren' num2str(gid)];
      eval(printer);

      clf; 
      semilogy(10000./lines.wnum,lines.stren*gas(gid)); title(num2str(gid));
      axis([10000/10000 10000/00300 1e-30*gas(gid) max(lines.stren*gas(gid))]);
      axis([10000/10000 10000/00300 1e-10 1e-0]); grid on;
      axis([10000/10000 20          1e-10 1e-0]); grid on;
      line([3.5 3.5],[1e-10 1e-0],'color','k','linewidth',2)
      line([15  15],[1e-10 1e-0],'color','k','linewidth',2)
      xlabel('Wavelength(um)'); ylabel('log10(linestren) * q(surf)')
      printer = ['print -dpng /home/sergio/PAPERS/SPECIALTALKS/IOWA08/FIGS/'];
      printer = [printer 'gasamtXstren' num2str(gid)];
      eval(printer);
      pause
      end
    end

elseif iWhich == 4  %% vis
  for gid = 1 : 30
    [lines,hitran_version] = hitread_simple(05000,30000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(10000./lines.wnum,lines.stren); title(num2str(gid));
      axis([10000/30000 10000/05000 1e-30 (max(lines.stren))]);
      xlabel('Wavelength(um)'); ylabel('log10(linestren)')
      pause
      end
    end

elseif iWhich == 3 %% nir
  for gid = 1 : 30
    [lines,hitran_version] = hitread_simple(03000,05000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(10000./lines.wnum,lines.stren); title(num2str(gid));
      axis([10000/5000 10000/3000 1e-30 (max(lines.stren))]);
      xlabel('Wavelength(um)'); ylabel('log10(linestren)')
      pause
      end
    end

elseif iWhich == 2  %% ir
  for gid = 1 : 30
    [lines,hitran_version] = hitread_simple(00500,03000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(10000./lines.wnum,lines.stren); title(num2str(gid));
      axis([10000/3000 10000/0500 1e-30 (max(lines.stren))]);
      xlabel('Wavelength(um)'); ylabel('log10(linestren)')
      pause
      end
    end

elseif iWhich == 1  %% mw
  for gid = 1 : 30
    [lines,hitran_version] = hitread_simple(00001,00500,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(10000./lines.wnum,lines.stren); title(num2str(gid));
      axis([10000/500 10000/001 1e-30 (max(lines.stren))]);
      xlabel('Wavelength(um)'); ylabel('log10(linestren)')
      pause
      end
    end

 end
