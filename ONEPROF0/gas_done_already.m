function  [freqchunk, numchunk, iCnt, iDone, iZeroCnt,iNumLines,iBroken] = gas_done_already(gasid,iWhichFreqBoundaries);

%% input
%%   iWhichFreqBoundaries = default (-1 for this database, or eg 2008 or 2012)
%%
%% output
%%   iCnt  = how many FULL chunks (each 11 files) need to be done
%%   iDone = how many FULL chunks (each 11 files) HAVE been  done
%%   iZeroCnt = how many empty files found
%%   iNumLines = number of lines in this band
%%   iBroken = how many BROKEN chunks (summed over all unfilled chunks, nonzero files) HAVE been  done

if nargin == 1
  iWhichFreqBoundaries = -1;  %% this one
end

iHitLinesFound = 0;
iNumLines = 0;
iBroken = 0;

figure(1); clf

addpath /home/sergio/SPECTRA
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools
addpath /asl/matlib/science
addpath /asl/matlib/aslutil

iZeroCnt = 0;
oo = [];

homedir = pwd;

nbox = 5; pointsPerChunk = 10000;

gid = gasid;
if iWhichFreqBoundaries == -1
  freq_boundaries
elseif iWhichFreqBoundaries == 2012
  freq_boundaries
  lala = strfind(dirout,'H2012_RUN8_NIRDATABASE');
  lenx = length('H2012_RUN8_NIRDATABASE');
  diroutX = dirout;
  diroutX(lala:lala+lenx-1) = 'H2012_RUN8_NIRDATABASE';
  dirout = diroutX;
elseif iWhichFreqBoundaries == 2008
  freq_boundaries
  lala = strfind(dirout,'H2012_RUN8_NIRDATABASE');
  lenx = length('H2012_RUN8_NIRDATABASE');
  diroutX = dirout;
  diroutX(lala:lala+lenx-1) = 'H2008_RUN8_NIRDATABASE';
  dirout = diroutX(1:end-5);    %% they do not have eg g2.dat but only g2
elseif iWhichFreqBoundaries == 2004
  freq_boundaries
  lala = strfind(dirout,'H2012_RUN8_NIRDATABASE');
  lenx = length('H2012_RUN8_NIRDATABASE');
  diroutX = dirout;
  diroutX(lala:lala+lenx-1) = 'H2004_RUN8_NIRDATABASE';
  dirout = diroutX(1:end-5);    %% they do not have eg g2.dat but only g2
end

fA = wn1; fB = wn2; df = dv;
dir0 = dirout;

%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/SPECTRA
[iYes,line] = findlines_plot(wn1-25,wn2+25,gid);

iaChunk = [];
iCnt = 0;
for wn = wn1 : dv : wn2
  woo = find(line.wnum >= wn-25 & line.wnum <= wn+dv+25);
  if length(woo) >= 1
    iCnt = iCnt + 1;
    iaChunk(iCnt) = wn;
    iaNumLinesInChunk(iCnt) = length(woo);
  end
end
iNumLines = length(line.wnum);

if length(iaChunk) == 0
  numchunk = 0;
  iCnt = 0;
  iDone = 0;
  freqchunk = [];
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf
%gasid = input('enter gas id : ');

%if gasid == 1
%  error('not for WV')
%else
%  iso = [-1];
%end

lines = show_vis_ir_lines_wavenumber(2012,7,gasid);

[mm,nn] = size(lines.iso);
if nn > 0
  figure(1);
  semilogy(lines.wnum,lines.stren);
  pause(0.1)
  dfx = 25;
  oo  = find(lines.wnum >= fA-df    &  lines.wnum <= fB+df);
  oo  = find(lines.wnum >= fA-dfx   &  lines.wnum <= fB+df+dfx);
  if length(oo) > 0
    semilogy(lines.wnum(oo),lines.stren(oo));
    if length(oo) > 1
      axis([fA-dfx fB+df+dfx min(lines.stren(oo)) max(lines.stren(oo))])
    else
      plot(lines.wnum(oo),lines.stren(oo),'o');
    end

    fprintf(1,'found %5i lines between %8.2f & %8.2f cm-1 \n',length(oo),fA,fB)

    fchunk = 0;
    for ff = fA : df : fB
      fchunk = fchunk + 1;  freqchunk(fchunk) = ff;
      hitlinesfound(fchunk) = 0;
      boof = find(lines.wnum >= ff &  lines.wnum <= ff+df);
      hitlinesfound(fchunk) = length(boof);
      if length(boof) == 0
        boof2 = find(lines.wnum >= ff-25 &  lines.wnum <= ff+df+25);  
        hitlinesfound(fchunk) = -length(boof2);  %%% << indicates this could be a chunk ie on the edge
      end
    end
  else
    fchunk = 0;
    for ff = fA : df : fB
      fchunk = fchunk + 1;  freqchunk(fchunk) = ff;
      hitlinesfound(fchunk) = 0;
    end  
  end
else
  fchunk = 0;
  for ff = fA : df : fB
    fchunk = fchunk + 1;  freqchunk(fchunk) = ff;
    hitlinesfound(fchunk) = 0;
  end  
end

cder = ['cd ' dir0]; eval(cder);

clear foundname
jj = 0;
jjx = 0;

if gasid >= 1 
  fchunk = 0;
  for ff = fA : df : fB
    fchunk = fchunk + 1;  freqchunk(fchunk) = ff;
    ifound = 0;
    for kk = 6
      if gid > 1
        fname = ['std' num2str(ff) '_' num2str(gasid) '_' num2str(kk) '.mat'];
      else
        fname = ['stdH2O' num2str(ff) '_' num2str(gasid) '_' num2str(kk) '_2.mat'];
      end
      dirr = dir(fname);
      if length(dirr) == 1 & dirr.bytes > 0
        jjx = +1;
        jj = jj + 1;
        ifound = ifound + 1;
        foundname(jj).name = dirr.name;
        foundname(jj).date = dirr.date;
      elseif length(dirr) == 1 & dirr.bytes == 0
        % comment out this ZERO file size warning
        % fprintf(1,'warning found file %s but size = 0 bytes! \n',fname)
        iZeroCnt = iZeroCnt + 1;
      end
    end
    numchunk(fchunk) = ifound;
  end 
end

gasidx = gid;

yesno = [];
cder = ['cd ' homedir]; eval(cder);
fprintf(1,'found %5i files for gasid %3i \n',jj,gasid);
if  jj > 0
  % comment out this tedious print of every file/date/time
  % for kk = 1 : jj
  %   fprintf(1,'%4i %3i %s %s\n',kk,gasid,foundname(kk).name,foundname(kk).date);
  % end
  for mm = 1 : length(freqchunk)
    if abs(hitlinesfound(mm)) > 0 & numchunk(mm) == 1
      yesno(mm) = +1;  %% 
    else
      yesno(mm) = 0;
    end
  end

  %% now print the results
  disp(' ')
  disp('gid  freqchunk hitlinesfound numchunk chunkdone(Y/N/M)')
  disp('------------------------------------------------------')

  for ii = 1 : iCnt
    woo = find(freqchunk == iaChunk(ii));  
    common(ii) = woo;
  end
  bwop = setdiff(1:length(freqchunk), common);
  if sum(abs(yesno(bwop))) > 0
    error(' >>>>>>>>>> WARNING : somehow found different chunk files made by run8, than what we expect from findlines_plot');
  end

  %% [gasid*ones(1,length(freqchunk)); freqchunk; hitlinesfound; numchunk; yesno]'
  wah = [gasid*ones(1,length(common)); freqchunk(common); hitlinesfound(common); numchunk(common); yesno(common)];
  fprintf(1,'%4i %8.2f     %6i        %3i    %3i \n',wah);
  fprintf(1,'found %5i lines between fA & fB cm-1 \n',length(oo))
  figure(1); title(num2str(gasid));
end

if iZeroCnt > 0
  fprintf(1,'  ---> WARNING found %4i files of size 0 bytes .... \n',iZeroCnt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(oo) > 0 & jjx == 0
  %% lines have been found in eg 2830 - 3530 cm-1 for this gas, but so far NOTHING completed
  for kk = 1 : jj
    fprintf(1,'%4i  %3i %s %s \n',kk,gasid,foundname(kk).name,foundname(kk).date);
  end

  disp('gid  freqchunk hitlinesfound numchunk')
  disp('--------------------------------------')
  %% [gasid*ones(1,length(freqchunk)); freqchunk; hitlinesfound; numchunk]'
  wah = [gasid*ones(1,length(freqchunk)); freqchunk; hitlinesfound; numchunk];
  fprintf(1,'%4i %8.2f     %6i        %3i    \n',wah);
  fprintf(1,'found %5i lines between fA & fB cm-1 \n',length(oo))
  figure(1); title(num2str(gasid));

  fprintf(1,' WARNING : gasid %3i has %5i lines w/in %8.2f,%8.2f cm-1 \n',gasid,length(oo),fA,fB)
  fprintf(1,'           but ZERO files found \n')

  for ii = 1 : iCnt
    woo = find(freqchunk == iaChunk(ii));
    iaChunkNum(ii) = numchunk(woo);
    if hitlinesfound(woo) > 0 | hitlinesfound(woo) < 0
      iHitLinesFound = iHitLinesFound + 1;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if jjx > 0
  %% lines have been found in eg 2830 - 3530 cm-1 for this gas, and SOMETHINGs already completed
  iDone = 0;
  iHitLinesFound = 0;
  for ii = 1 : iCnt
    woo = find(freqchunk == iaChunk(ii));
    iaChunkNum(ii) = numchunk(woo);
    if yesno(woo) == 1
      iDone = iDone + 1;
      iaChunkDone(ii) = +1;
    else
      iaChunkDone(ii) = -1;
    end

    if hitlinesfound(woo) > 0 | hitlinesfound(woo) < 0
      iHitLinesFound = iHitLinesFound + 1;
    end

  end
  figure(2); clf; bar(iaChunk,iaChunkNum);
  figure(2); title(num2str(gasid));
else
  iDone = 0;
end

iCnt = iHitLinesFound;   %% !!

%%   iCnt  = how many chunks need to be done
%%   iDone = how many chunks HAVE been  done
%%   iZeroCnt = how many empty files found
%%   iNumLines = number of lines in this band

NumFilesFullChunks = iDone * 11;
iBrokenChunks = jj - NumFilesFullChunks - iZeroCnt; %% jj = numfilesfound
iBrokenChunks = max(iBrokenChunks/11,0);
fprintf(1,'%2i chunks in this band; potentially %2i for gasid %2i, [completed full, completed broken] = %2i %5.3f \n',iCnt,iHitLinesFound,gasid,iDone,iBrokenChunks)

iBroken = iBrokenChunks;