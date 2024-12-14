clear all
gidlist = 51 : 81;
for gidxx = 1 : length(gidlist)
  gid = gidlist(gidxx);

  [freqchunk, numchunk, iCnt(gid), iDone(gid)] = xsec_done_already(gid);
  %disp('ret to continue : '); pause

  figure(3);
  plot(1:gid,iCnt,'bo-',1:gid,iDone,'rx-'); 
    xlabel('GasID'); ylabel('Number chunks needed/done'); grid on

  pause(0.1);
end

figure(3);
plot(gidlist,iCnt(gidlist),'bo-',gidlist,iDone(gidlist),'rx-'); 
    xlabel('GasID'); ylabel('Number chunks needed/done'); grid on
