function [theprofs] = rtp_prof_to_oneprof1(fname,iY);

addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools

%% fname = '/asl/s1/sergio/RTP_pin_feb2002/pin_feb2002_sea_airsnadir_op.sun.rtp';
%% fname = input('Enter layers rtp file : ');

ee = exist(fname);
if ee > 0
  [h,ha,p,pa] = rtpread(fname);
else
  fprintf(1,'%s dne \n',fname);
end

gunit = h.gunit;
if gunit(1) ~= 1 & h.ptype < 1
  error('does not look like LAYERS file!!!');
else
  fprintf(1,'there are %6i profiles in %s \n',length(p.stemp),fname);
end

%% iY = input('enter which profile : ');

nlays = p.nlevs(iY) - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%
a0 = load('/home/sergio/HITRAN2UMBCLBL/REFPROF/refproTRUE.mat');

theprofs.glist = h.glist;
theprofs.mpres = a0.refpro.mpres;   %%% assuming we are using DEFAULT AIRS pressure layers
theprofs.plev  = a0.refpro.plev;

plevs = p.plevs(1:nlays+1,iY)/1013.25;  %% change to atm
theprofs.plev = plevs;

mpresN = plevs(1:end-1) - plevs(2:end);
mpresD = log(plevs(1:end-1) ./ plevs(2:end));
theprofs.mpres = mpresN ./ mpresD;

L = p.palts(1:nlays+1,iY); L = abs(diff(L));   %% in meters

theprofs.mtemp = p.ptemp(1:nlays,iY);
for gg = 1 : length(h.glist)
  xer = ['x = p.gas_' num2str(h.glist(gg)) ';'];
  eval(xer);
  x = x(1:nlays,iY);                 %% molecules/cm2
  theprofs.gamnt(:,gg) = x/6.023e26;   %% change to kmol/cm2

  %% how to find partial pressure???
  %% pV = nRT ==> n/V = p/RT ==> gas amount q = L n/V = p/RT L
  %% p = q R T /l       

  x = x*1e4/6.023e23;  %% change from molecules/cm2 to mol/m2
  pp = x * 8.31 .* theprofs.mtemp ./L;  %% in N/m2
  theprofs.gpart(:,gg) = pp/1.01325e5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%

if length(a0.refpro.glist) > length(h.glist)
  haha = [length(a0.refpro.glist)  length(h.glist)];
  fprintf(1,'  The reference profile had %2i gases while your rtp file had %2i gases \n',haha);
  %iX = input('  Add on the extra gases from the default profile??? (-1 = no +1 = yes) : ');
  iX = +1;
  if iX > 0
    disp('adding on extra gases from STD')
    %% first index the pressure layers
    for ll = 1 : length(theprofs.mpres)
      pdiff = abs(a0.refpro.mpres - theprofs.mpres(ll));
      layer_index(ll) = find(pdiff == min(pdiff));
    end
    figure(4);
    plot(theprofs.mpres,a0.refpro.mpres(layer_index),theprofs.mpres,theprofs.mpres-a0.refpro.mpres(layer_index))
    xlabel('rtp players'); ylabel('reference pro players');

    %% then go through the different gases
    iCnt = length(theprofs.glist);
    [theset,I] = setdiff(a0.refpro.glist,h.glist);
    for gg = 1 : length(theset)
      iCnt = iCnt + 1;
      theprofs.glist(iCnt) = theset(gg);
      theprofs.gamnt(:,iCnt) = a0.refpro.gamnt(layer_index,I(gg));
      theprofs.gpart(:,iCnt) = a0.refpro.gpart(layer_index,I(gg));
    end

    %% now sort the gas indices
    [Y,I] = sort(theprofs.glist);
    theprofs.glist = theprofs.glist(I);
    theprofs.gamnt = theprofs.gamnt(:,I);
    theprofs.gpart = theprofs.gpart(:,I);

  end
end

%%% save oneprof.mat theprofs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); semilogy(a0.refpro.mtemp,a0.refpro.mpres,'b',theprofs.mtemp,theprofs.mpres,'r')
  set(gca,'ydir','reverse'); xlabel('T(K)'); ylabel('pres (atm)'); 
  legend('reference',['profile ' num2str(iY)],'location','best','fontsize',10); grid
  title('T(p)')

figure(2); loglog(a0.refpro.gamnt(:,1),a0.refpro.mpres,'b',theprofs.gamnt(:,1),theprofs.mpres,'r')
  set(gca,'ydir','reverse'); xlabel('gas amt kmol/cm2'); ylabel('pres (atm)'); 
  legend('reference',['profile ' num2str(iY)],'location','best','fontsize',10); grid
  title('WV amount (p)')

figure(3); loglog(a0.refpro.gpart(:,1),a0.refpro.mpres,'b',theprofs.gpart(:,1),theprofs.mpres,'r')
  set(gca,'ydir','reverse'); xlabel('gas partial pres (atm)'); ylabel('pres (atm)'); 
  legend('reference',['profile ' num2str(iY)],'location','best','fontsize',10); grid
  title('WV partial press (p)')
