function lines = show_vis_ir_lines_wavenumber_individualgas(iWhich,gasids,fplots,hit_version);

% lines = show_vis_ir_lines_wavenumber_individualgas(iWhich,gasids,fplots,hitversion);
% where
% iWhich = +1 +2 +3  +4     +5           +6        +7     +8       +9 
%          MW IR NIR VIS (FIR-NIR) (NIR-VSWIR) (IR-VIS)  (VIS-UV) [f1 f2]
% fplots = [f1 f2] if iWhich == 9, not needed otherwise
%  optional hit)version = 2012 (default), else eg 2004, 2008 

lines = [];

hitranpath
if nargin == 3
  fname = '/asl/data/hitran/h24.by.gas/';
  fname = [HITRAN 'h24.by.gas/'];  
elseif hit_version == 2020
  fname = '/asl/data/hitran/h20.by.gas/';
  fname = [HITRAN 'h20.by.gas/'];    
elseif hit_version == 2016
  fname = '/asl/data/hitran/h16.by.gas/';
  fname = [HITRAN 'h16.by.gas/'];    
elseif hit_version == 2012
  fname = '/asl/data/hitran/h12.by.gas/';
  fname = [HITRAN 'h12.by.gas/'];    
elseif hit_version == 2008
  fname = '/asl/data/hitran/h08.by.gas/';
  fname = [HITRAN 'h08.by.gas/'];    
elseif hit_version == 2004
  fname = '/asl/data/hitran/h04.by.gas/';
  fname = [HITRAN 'h04.by.gas/'];    
elseif hit_version == 2000
  fname = '/asl/data/hitran/h2k.by.gas/';
  fname = [HITRAN 'h2k.by.gas/'];    
elseif hit_version == 1996
  fname = '/asl/data/hitran/h96.by.gas/';
  fname = [HITRAN 'h96.by.gas/'];    
else
  error('only have H1996,2000,2004,2008,2012,2016,2020,2024')
end



iNew = +1;
if iNew < 0
  if intersect(gasids,32)
    error('do not have p.gas_32 in pin_feb2002_sea_airsnadir_op.22deg.rtp');
    end

  fop = '/carrot/s1/sergio/pin_feb2002_sea_airsnadir_op.22deg.rtp';
  fop = '/asl/s1/sergio/pin_feb2002_sea_airsnadir_op.22deg.rtp';
  [h,ha,p,pa] = rtpread(fop);
  %% recall prof 49 = US Std
  for ii = 1 : 31
    str = ['doink = p.gas_' num2str(ii,'%2d') '(4,49);']; eval(str);
    gas(ii) = doink;  %% in molecules/cm2 at surface
    end
  end

if iNew > 0
  %load /asl/packages/klayersV205/Data/refprof_usstd2010.mat
  load /asl/packages/klayersV205/Data/refprof_usstd16Aug2010.mat
  glistall = [[1:42] [51:81]];
  for ii = 1 : length(glistall)
    gas = glistall(ii);
    str = ['doink = prof.gas_' num2str(gas) '(4);']; eval(str);
    gas(ii) = doink;  %% in molecules/cm2 at surface
    end
  end

%iWhich = input(' +1 +2 +3 +4 +5 +6 +7 : MW IR NIR VIS (FIR-NIR) (NIR-VSWIR) (IR-VIS) : ');

switch iWhich
  case  9
    f1 = fplots(1);     f2 = fplots(2);   %%% generic special
  case  8
    f1 = 12500; f2 = 50000; %% all bands, from VIS -- UV
  case  7
    f1 = 00001; f2 = 50000; %% all bands, from FIR -- UV
  case  6
    f1 = 02500; f2 = 06000; %% from nir-vswir
  case  5
    f1 = 00300; f2 = 10000; %% from fir - nir
  case  4
    f1 = 05000; f2 = 30000; %% vis
  case  3
    f1 = 03000; f2 = 05000; %% nir
  case  2
    f1 = 00500; f2 = 03000; %% ir
  case  1
    f1 = 00001; f2 = 00500; %% fir
  end

for gid = 1 : length(gasids)
  gasid = gasids(gid);
  [lines,hitran_version] = hitread_simple(f1,f2,0,gasid,fname);
  if lines.linct > 0
    clf; semilogy(lines.wnum,lines.stren); title(num2str(gasid));
    if (max(lines.stren)) > 1e-30
      axis([f1 f2 1e-30 (max(lines.stren))]);
    else
      axis([f1 f2 max(lines.stren)/10 (max(lines.stren))]);
    end
    xlabel('Wavenumber(cm-1)'); ylabel('log10(linestren)')
    pause(0.1)
    end
  end

