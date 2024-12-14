addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools

fname = input('Enter layers rtp file : ');
%% fname = '/asl/s1/sergio/pin_feb2002_sea_airsnadir_op.sun.rtp';

ee = exist(fname);
if ee > 0
  [h,ha,p,pa] = rtpread(fname);
else
  fprintf(1,'%s dne \n',fname);
end

gunit = h.gunit;
if gunit(1) ~= 1 & h.ptype < 1
  error('does not look like LATYERS file!!!');
else
  fprintf(1,'there are %6i profiles in %s \n',length(p.stemp),fname);
end

iX = input('enter which profile : ');

nlays = p.nlevs(iX) - 1;

a0 = load('/home/sergio/HITRAN2UMBCLBL/REFPROF/refproTRUE.mat');

refpro.glist = h.glist;
refpro.mpres = a0.refpro.mpres;   %%% assuming we are using DEFAULT AIRS pressure layers
refpro.plev  = a0.refpro.plev;

plevs = p.plevs(1:nlays+1,iX)/1013.25;  %% change to atm
refpro.plev = plevs;

mpresN = plevs(1:end-1) - plevs(2:end);
mpresD = log(plevs(1:end-1) ./ plevs(2:end));
refpro.mpres = mpresN ./ mpresD;

L = p.palts(1:nlays+1,iX); L = abs(diff(L));   %% in meters

refpro.mtemp = p.ptemp(1:nlays,iX);
for gg = 1 : length(h.glist)
  xer = ['x = p.gas_' num2str(h.glist(gg)) ';'];
  eval(xer);
  x = x(1:nlays,iX);                 %% molecules/cm2
  refpro.gamnt(:,gg) = x/6.023e26;   %% change to kmol/cm2

  %% how to find partial pressure???
  %% pV = nRT ==> n/V = p/RT ==> gas amount q = L n/V = p/RT L
  %% p = q R T /l       

  x = x*1e4/6.023e23;  %% change from molecules/cm2 to mol/m2
  pp = x * 8.31 .* refpro.mtemp ./L;  %% in N/m2
  refpro.gpart(:,gg) = pp/1.01325e5;
end

if length(a0.refpro.glist) > length(h.glist)
  haha = [length(a0.refpro.glist)  length(h.glist)];
  fprintf(1,'  The reference profile had %2i gases while your rtp file had %2i gases \n',haha);
  iX = input('  Add on the extra gases from the default profile??? (-1 = no +1 = yes) : ');
  if iX > 0
    %% first index the pressure layers
    for ll = 1 : length(refpro.mpres)
      pdiff = abs(a0.refpro.mpres - refpro.mpres(ll));
      layer_index(ll) = find(pdiff == min(pdiff));
    end
    figure(4);
    plot(refpro.mpres,a0.refpro.mpres(layer_index),refpro.mpres,refpro.mpres-a0.refpro.mpres(layer_index))
    xlabel('rtp players'); ylabel('reference pro players');

    %% then go through the different gases
    iCnt = length(refpro.glist);
    [theset,I] = setdiff(a0.refpro.glist,h.glist);
    for gg = 1 : length(theset)
      iCnt = iCnt + 1;
      refpro.glist(iCnt) = theset(gg);
      refpro.gamnt(:,iCnt) = a0.refpro.gamnt(layer_index,I(gg));
      refpro.gpart(:,iCnt) = a0.refpro.gpart(layer_index,I(gg));
    end

    %% now sort the gas indices
    [Y,I] = sort(refpro.glist);
    refpro.glist = refpro.glist(I);
    refpro.gamnt = refpro.gamnt(:,I);
    refpro.gpart = refpro.gpart(:,I);

  end
end

save oneprof.mat refpro

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); semilogy(a0.refpro.mtemp,a0.refpro.mpres,'b',refpro.mtemp,refpro.mpres,'r')
  set(gca,'ydir','reverse'); xlabel('T(K)'); ylabel('pres (atm)'); 
  title('(b) reference (r) selected profile'); grid

figure(2); loglog(a0.refpro.gamnt(:,1),a0.refpro.mpres,'b',refpro.gamnt(:,1),refpro.mpres,'r')
  set(gca,'ydir','reverse'); xlabel('gas amt kmol/cm2'); ylabel('pres (atm)'); 
  title('WV (b) reference (r) selected profile'); grid

figure(3); loglog(a0.refpro.gpart(:,1),a0.refpro.mpres,'b',refpro.gpart(:,1),refpro.mpres,'r')
  set(gca,'ydir','reverse'); xlabel('gas partial pres (atm)'); ylabel('pres (atm)'); 
  title('WV (b) reference (r) selected profile'); grid

