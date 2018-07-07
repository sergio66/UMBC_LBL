function lines = show_vis_ir_lines_wavenumber(hitV,iWhich,gid0);

%% input
%%   hitV = hitran version
%%   iWHich = +1 +2 +3 +4 +5 +6 +7 : MW IR NIR VIS (FIR-NIR) (NIR-VSWIR) (IR-VIS) 
%%   gid0   = gasID

%hitV = input('enter HITRAN version (2004, 2008, 2012) : ');

if hitV == 1992
  fname = '/asl/data/hitran/h92.by.gas/';
elseif hitV == 1996
  fname = '/asl/data/hitran/h96.by.gas/';
elseif hitV == 1998
  fname = '/asl/data/hitran/h98.by.gas/';
elseif hitV == 2000
  fname = '/asl/data/hitran/h2k.by.gas/';
elseif hitV == 2004
  fname = '/asl/data/hitran/h04.by.gas/';
elseif hitV == 2008
  fname = '/asl/data/hitran/h08.by.gas/';
elseif hitV == 2012
  fname = '/asl/data/hitran/h12.by.gas/';
elseif hitV == 2016
  fname = '/asl/data/hitran/h16.by.gas/';
else
  error('1992, 1996, 1998, 2000, 2004, 2008, 2012, 2016 ...')
end

%{
fop = '/strowdata1/s1/sergio/pin_feb2002_sea_airsnadir_op.22deg.rtp';
[h,ha,p,pa] = rtpread(fop);
for ii = 1 : 30
  str = ['doink = p.gas_' num2str(ii,'%2d') '(4,49);']; eval(str);
  gas(ii) = doink;  %% in molecules/cm2 at surface
end
%}

% iWhich = input(' +1 +2 +3 +4 +5 +6 +7 : MW IR NIR VIS (FIR-NIR) (NIR-VSWIR) (IR-VIS) : ');
if iWhich == 7   %% all
  for gid = gid0 : gid0    %% 1 : 30
    %[lines,hitran_version] = hitread_simple(00001,30000, 0,gid,fname);
    [lines,hitran_version] = hitread_simple(00000,300000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([01 30000 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')
      pause(0.1)
    end
  end

elseif iWhich == 6   %% nir-vswir
  for gid = gid0 : gid0    %% 1 : 30
    [lines,hitran_version] = hitread_simple(02500,05000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([02500 05000 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')

      clf; semilogy(lines.wnum,lines.stren*gas(gid)); title(num2str(gid));
      axis([2500 5000 1e-30*gas(gid) max(lines.stren*gas(gid))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren) * q(surf)')

      pause(0.1)
    end
  end

elseif iWhich == 5   %% fir-nir
  for gid = gid0 : gid0    %% 1 : 30
    [lines,hitran_version] = hitread_simple(00300,10000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([00300 10000 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')

      clf; semilogy(lines.wnum,lines.stren*gas(gid)); title(num2str(gid));
      axis([300 10000 1e-30*gas(gid) max(lines.stren*gas(gid))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren) * q(surf)')

      pause(0.1)
    end
  end

elseif iWhich == 4  %% vis
  for gid = gid0 : gid0    %% 1 : 30
    [lines,hitran_version] = hitread_simple(05000,30000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([5000 30000 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')
      pause(0.1)
    end
  end

elseif iWhich == 3 %% nir
  for gid = gid0 : gid0    %% 1 : 30
    [lines,hitran_version] = hitread_simple(02800,06000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([3000 5000 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')
      pause(0.1)
    end
  end

elseif iWhich == 2  %% ir
  for gid = gid0 : gid0    %% 1 : 30
    [lines,hitran_version] = hitread_simple(00500,03000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([500 3000 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')
      pause(0.1)
    end
  end

elseif iWhich == 1  %% mw
  for gid = gid0 : gid0    %% 1 : 30
    [lines,hitran_version] = hitread_simple(00001,00500,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([01 500 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')
      pause(0.1)
    end
  end

end