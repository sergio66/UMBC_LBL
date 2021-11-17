function lines = show_vis_ir_lines_wavenumber_water(hitV,iWhich,gasid,gasisotopes);

if gasid == 32
  error('do not have gas_32 in pin_feb2002_sea_airsnadir_op.22deg.rtp');
end

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
elseif hitV == 2020
  fname = '/asl/data/hitran/h20.by.gas/';
else
  hitV
  error('1992, 1996, 1998, 2000, 2004, 2008, 2012, 2016, 2020 ..')
end

lines = [];

%fop = '/carrot/s1/sergio/pin_feb2002_sea_airsnadir_op.22deg.rtp';
%[h,ha,p,pa] = rtpread(fop);
%for ii = 1 : 31
%  str = ['doink = p.gas_' num2str(ii,'%2d') '(4,49);']; eval(str);
%  gas(ii) = doink;  %% in molecules/cm2 at surface
%end

%iWhich = input(' +1 +2 +3 +4 +5 +6 +7 : MW IR NIR VIS (FIR-NIR) (NIR-VSWIR) (IR-VIS) : ');

if iWhich ~= 7
  error('can only handle iWhich = 7');
end

if iWhich == 7   %% all
  for gid = gasid : gasid
    [lines,hitran_version] = hitread_simple(00001,30000,0,gid,fname);
    shoop = [];
    for jj = 1 : length(gasisotopes)
      woop = find(lines.iso == gasisotopes(jj));
      if lines.linct > 0 & length(woop) > 0
        shoop = [shoop woop];
      end
    end
    woop = shoop;
    clf; semilogy(lines.wnum(woop),lines.stren(woop)); 
    title(['gasid = ' num2str(gid) '; iso = ' num2str(gasisotopes)]); pause(0.1)
    axis([01 50000 1e-30 (max(lines.stren(woop)))]);
    xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')

    lines.wnum = lines.wnum(woop);
    lines.stren = lines.stren(woop);
    lines.iso = lines.iso(woop);
  end

elseif iWhich == 6   %% nir-vswir
  for gid = gasid : gasid
    [lines,hitran_version] = hitread_simple(02500,06000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([02500 05000 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')

      clf; semilogy(lines.wnum,lines.stren*gas(gid)); title(num2str(gid));
      axis([2500 5000 1e-30*gas(gid) max(lines.stren*gas(gid))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren) * q(surf)')

    end
  end

elseif iWhich == 5   %% fir-nir
  for gid = gasid : gasid
    [lines,hitran_version] = hitread_simple(00300,10000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([00300 10000 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')

      clf; semilogy(lines.wnum,lines.stren*gas(gid)); title(num2str(gid));
      axis([300 10000 1e-30*gas(gid) max(lines.stren*gas(gid))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren) * q(surf)')

    end
  end

elseif iWhich == 4  %% vis
  for gid = gasid : gasid
    [lines,hitran_version] = hitread_simple(05000,30000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([5000 30000 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')
    end
  end

elseif iWhich == 3 %% nir
  for gid = gasid : gasid
    [lines,hitran_version] = hitread_simple(03000,05000,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([3000 5000 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')
    end
  end

elseif iWhich == 2  %% ir
  for gid = gasid : gasid
    [lines,hitran_version] = hitread_simple(00500,03000,0,gid,fname);
%    if lines.linct > 0
%      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
%      axis([500 3000 1e-30 (max(lines.stren))]);
%      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')
%    end
%  end
    shoop = [];
    for jj = 1 : length(gasisotopes)
      woop = find(lines.iso == gasisotopes(jj));
      if lines.linct > 0 & length(woop) > 0
        shoop = [shoop woop];
      end
    end
    woop = shoop;
    clf; semilogy(lines.wnum(woop),lines.stren(woop)); 
    title(['gasid = ' num2str(gid) '; iso = ' num2str(gasisotopes)]); pause(0.1)
    axis([01 50000 1e-30 (max(lines.stren(woop)))]);
    xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')

    lines.wnum = lines.wnum(woop);
    lines.stren = lines.stren(woop);
    lines.iso = lines.iso(woop);
  end

elseif iWhich == 1  %% mw
  for gid = gasid : gasid    
    [lines,hitran_version] = hitread_simple(00001,00500,0,gid,fname);
    if lines.linct > 0
      clf; semilogy(lines.wnum,lines.stren); title(num2str(gid));
      axis([01 500 1e-30 (max(lines.stren))]);
      xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')
    end
  end

end
