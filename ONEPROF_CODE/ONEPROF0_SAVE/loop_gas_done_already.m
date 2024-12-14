if exist('gidlist') & exist('iCnt') & exist('iDone')
  gidlist0 = gidlist;
  iCnt0    = iCnt;
  iDone0   = iDone;
  iBroken0 = iBroken;
  disp('  saving off gidlist iCnt iDone iBroken')
  clear gidlist iCnt iDone iX iY gid gidxx iBroken
end

gidlist = setdiff(1:47,[30 35 42]);   %% those last 3 are NOT in the breakout, and 
                                      %%% OOPS no profiles for 43,44,45,46,47
gidlist = setdiff(1:42,[30 35 42]);   %% those last 3 are NOT in the breakout

iX = input('Enter gas list (-1 for all from 1 to 42, [a b] to choose, and +N to specify ONE) : ');
if length(iX) == 1 & iX == -1
  gidlist = 1 : 42;
elseif length(iX) == 1 & iX > 0
  gidlist = iX;
else
  gidlist = iX(1) : iX(2);
end

gidlist = setdiff(gidlist,[30 35 42]);   %% last 3 are NOT in the breakout, 
                                         %% and OOPS no profiles for 43,44,45,46,47

iY = input('Enter year (-1) for this one, or 2008 or 2012 : ');

%%%%%%%%%%%%%%%%%%%%%%%%%

for gidxx = 1 : length(gidlist)
  gid = gidlist(gidxx);

  [freqchunk, numchunk, iCnt(gid), iDone(gid), iZeroCnt(gid), iNumLines(gid), iBroken(gid)] = gas_done_already(gid,iY);
  %disp('ret to continue : '); pause

  figure(3);
  plot(1:gid,iCnt,'bo-',1:gid,iDone,'rx-',1:gid,iBroken,'k'); 
    xlabel('GasID'); 
    ylabel('Number chunks needed(b)/done(r), potentially done(k)'); grid on

  pause(0.1);
  %disp('<RET> to continue'); pause
end

figure(2);
plot(gidlist,iNumLines(gidlist),'bo-'); title('NumLines');

figure(3);
plot(gidlist,iCnt(gidlist),'bo-',gidlist,iDone(gidlist),'rx-',gidlist,iZeroCnt(gidlist)/11,'gd-',...
     gidlist,iBroken(gidlist),'k','linewidth',3); 
hold on
plot(gidlist,iDone(gidlist) + iBroken(gidlist),'rd--','linewidth',4)
hold off
   xlabel('GasID')
   ylabel('Number chunks needed(b)/done(r), potentially done(k)'); grid on
ax = axis;
axis([ax(1) ax(2) ax(3) ax(4)+1])
title('(g) number of chunks being worked on');

fprintf(1,'require total of %4i chunks (summed over all gases) to be done \n',sum(iCnt))
fprintf(1,'  chunks done so far (summed over all gases) = %4i \n',sum(iDone))
fprintf(1,'  empty chunks/empty files = %6.2f/%4i \n',sum(iZeroCnt(gidlist)/11),sum(iZeroCnt(gidlist)))

if exist('gidlist0') & exist('iCnt0') & exist('iDone0')
  hold on
  plot(gidlist0,iCnt0(gidlist0),'co-',gidlist0,iDone0(gidlist0),'ms-',gidlist,iBroken0(gidlist),'y'); 
    xlabel('GasID'); ylabel('Number chunks needed/done'); grid on
     ylabel('Number chunks needed(b)/done(r), potentially done(k)'); grid on
  hold off
  fprintf(1,'  chunks previsouly done (summed over all gases) = %4i \n',sum(iDone0))
  title('(g) number of chunks being worked on \newline (c/m) old chunks needed/done');
end


disp(' ')
disp('Figure 3 : blue  = chunks needed');
disp('         : red   = chunks done (red/blue overlap when everything done)');
disp('         : green = empty files found (divided by 11), implies being worked on')
disp('         : black = broken chunks done ie unfilled chunks, but files already done)')
disp(' ')
disp('         : cyan    = previous iteration, chunks needed')
disp('         : magenta = previous iteration, chunks done')
disp('         : yellow  = previous iteration, broken chunks done');
disp('  ')
disp('So magenta --> red, yellow --> black = signs of progress!!!')
