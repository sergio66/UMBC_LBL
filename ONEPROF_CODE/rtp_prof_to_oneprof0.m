addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools

%% fname = '/asl/s1/sergio/pin_feb2002_sea_airsnadir_op.sun.rtp';
%% fname = '/home/chepplew/projects/klayers_wrk/regr49_pbl.op.rtp';

set_file_names
%fname = input('Enter layers rtp file : ');
%iProf = input('enter which profile : ');

fname = frtp; 
fprintf(1,'will look at profile %4i from rtp file = %s \n',iProf,fname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ee = exist(fname);
if ee > 0
  [h,ha,p,pa] = rtpread(fname);
else
  fprintf(1,'%s dne \n',fname);
end

addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
if ~isfield(p,'plays')
  disp('adding p.plays')
  p = make_rtp_plays(p);
end

gunit = h.gunit;
if gunit(1) ~= 1 & h.ptype < 1
  error('does not look like LAYERS file!!!');
else
  fprintf(1,'there are %6i profiles in %s \n',length(p.stemp),fname);
end

nlays = p.nlevs(iProf) - 1;

a0 = load('/home/sergio/HITRAN2UMBCLBL/REFPROF/refproTRUE.mat');

rtpProf.glist = h.glist;

%% assuming we are using DEFAULT AIRS pressure layers
%% rtpProf.mpres = a0.refpro.mpres;   
%% rtpProf.plev  = a0.refpro.plev;

plevs = p.plevs(1:nlays+1,iProf)/1013.25;  %% change to atm
rtpProf.plev = plevs;                      %% going to have one more level than mpres
mpresN = plevs(1:end-1) - plevs(2:end);
mpresD = log(plevs(1:end-1) ./ plevs(2:end));
rtpProf.mpres = mpresN ./ mpresD;
rtpProf.mpres = p.plays(1:nlays,iProf)/1013.25;

L = p.palts(1:nlays+1,iProf); L = abs(diff(L));   %% in meters

rtpProf.mtemp = p.ptemp(1:nlays,iProf);
for gg = 1 : length(h.glist)
  xer = ['x = p.gas_' num2str(h.glist(gg)) ';'];
  eval(xer);
  x = x(1:nlays,iProf);                 %% molecules/cm2
  rtpProf.gamnt(:,gg) = x/6.023e26;   %% change to kmol/cm2

  %% how to find partial pressure???
  %% pV = nRT ==> n/V = p/RT ==> gas amount q = L n/V = p/RT L
  %% p = q R T /l       

  x = x*1e4/6.023e23;  %% change from molecules/cm2 to mol/m2
  pp = x * 8.31 .* rtpProf.mtemp ./L;  %% in N/m2
  rtpProf.gpart(:,gg) = pp/1.01325e5;
end

if length(a0.refpro.glist) > length(h.glist)
  haha = [length(a0.refpro.glist)  length(h.glist)];
  fprintf(1,'  The reference profile had %2i gases while your rtp file had %2i gases \n',haha);
  iX = input('  Add on the extra gases from the default profile??? (-1/DEFAULT = no +1 = yes) : ');
  if length(iX) == 0
    iX = -1;
  end
  if iX > 0
    %% first index the pressure layers
    for ll = 1 : length(refpro.mpres)
      pdiff = abs(a0.refpro.mpres - refpro.mpres(ll));
      layer_index(ll) = find(pdiff == min(pdiff));
    end
    figure(4);
    plot(rtpProf.mpres,a0.refpro.mpres(layer_index),rtpProf.mpres,rtpProf.mpres-a0.refpro.mpres(layer_index))
    xlabel('rtp players'); ylabel('reference pro players');

    %% then go through the different gases
    iCnt = length(rtpProf.glist);
    [theset,I] = setdiff(a0.refpro.glist,h.glist);
    for gg = 1 : length(theset)
      iCnt = iCnt + 1;
      rtpProf.glist(iCnt) = theset(gg);
      rtpProf.gamnt(:,iCnt) = a0.refpro.gamnt(layer_index,I(gg));
      rtpProf.gpart(:,iCnt) = a0.refpro.gpart(layer_index,I(gg));
    end

    %% now sort the gas indices
    [Y,I] = sort(rtpProf.glist);
    rtpProf.glist = rtpProf.glist(I);
    rtpProf.gamnt = rtpProf.gamnt(:,I);
    rtpProf.gpart = rtpProf.gpart(:,I);

  end
end

rtpProf.rtp = fname;
rtpProf.iProf = iProf;

junk = findstr(fname,'/');
junk = junk(end)+1;
junkname = fname(junk:end);
junk = findstr(junkname,'.');
junk= junk(1);
junkname = junkname(1:junk-1);

outfile = ['PROFILES/oneprof_' junkname '_' num2str(iProf) '.mat'];
saver = ['save ' outfile ' rtpProf' ];                      eval(saver)
lner  = ['!rm oneprof.mat; ln -s ' outfile ' oneprof.mat']; eval(lner)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); semilogy(a0.refpro.mtemp,a0.refpro.mpres,'b',rtpProf.mtemp,rtpProf.mpres,'r')
  set(gca,'ydir','reverse'); xlabel('T(K)'); ylabel('pres (atm)'); 
  title('(b) reference (r) selected profile'); grid

figure(2); loglog(a0.refpro.gamnt(:,1),a0.refpro.mpres,'b',rtpProf.gamnt(:,1),rtpProf.mpres,'r')
  set(gca,'ydir','reverse'); xlabel('gas amt kmol/cm2'); ylabel('pres (atm)'); 
  title('WV (b) reference (r) selected profile'); grid

figure(3); loglog(a0.refpro.gpart(:,1),a0.refpro.mpres,'b',rtpProf.gpart(:,1),rtpProf.mpres,'r')
  set(gca,'ydir','reverse'); xlabel('gas partial pres (atm)'); ylabel('pres (atm)'); 
  title('WV (b) reference (r) selected profile'); grid

