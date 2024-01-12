function  lineNEW = translate2oldHITparams(lineOLD);

heldplot = -1;

if length(lineOLD.igas) == 0
  lineOLD
  error('pah : lineOLD looks empty')
end

lala = lineOLD.ilsgq;
[m,n] = size(lala);
if n ~= lineOLD.linct
  error('in translate2oldHITparams,  n ~= lineOLD.linct');
end

if m == 1
  disp('    >>> in translate2oldHITparams, NO NEED TO change iusgq and ilsgq into old formats for CO2!!! <<<< ')
  lineNEW = lineOLD;
  return
end

disp('    >>> in translate2oldHITparams, changing iusgq and ilsgq into old formats for CO2!!! <<<< ')

if (lineOLD.igas(1) ~= 2 & lineOLD.igas(1) ~= 5)
  error('expecting gasID == 2 or 5')
end

%{
%% for example for CO2
aa = lineOLD.iusgq;
size(aa) = 15 x linescnt
and they look like
  15x3 char array

    '   '
    '   '
    '   '
    '   '
    '   '
    '   '
    '   '
    '462'
    '   '
    '456'
    '   '
    '456'
    '   '
    '100'
    '342'
So you see only 5 of the columns should have any quantum IDs
%}

%% same co2vibs for both H2004 and H2008 so this should be fine
%%see eg /home/sergio/SPECTRA/Global_Data_HITRAN2004/BD_Vibs.for
if lineOLD.igas(1) == 2
  %xvibs = load('/home/sergio/SPECTRA/Global_Data_HITRAN2008/co2vibs');
  xvibs = load('/home/sergio/SPECTRA/Global_Data_HITRAN2004/co2vibs');
  lenID = 5;
elseif lineOLD.igas(1) == 5
  xvibs = load('/home/sergio/SPECTRA/Global_Data_HITRAN2004/covibs');
  lenID = 1;  
end

lineNEW.igas   = lineOLD.igas;
lineNEW.iso    = lineOLD.iso;
lineNEW.wnum   = lineOLD.wnum;
lineNEW.stren  = lineOLD.stren;
lineNEW.tprob  = lineOLD.tprob;
lineNEW.abroad = lineOLD.abroad;
lineNEW.sbroad = lineOLD.sbroad;
lineNEW.els    = lineOLD.els;
lineNEW.abcoef = lineOLD.abcoef;
lineNEW.tsp    = lineOLD.tsp;
lineNEW.ai     = lineOLD.ai;
lineNEW.ref    = lineOLD.ref;
lineNEW.linct  = lineOLD.linct;
lineNEW.gasid  = lineOLD.gasid;

%lineNEW.uslq = [lineOLD.uslq(:,2:9) lineOLD.uslq(:,1)];
%lineNEW.bslq = [lineOLD.bslq(:,2:9) lineOLD.bslq(:,1)];
%lineNEW.uslq = [lineOLD.uslq(:,3:9) lineOLD.uslq(:,1:2)];
%lineNEW.bslq = [lineOLD.bslq(:,3:9) lineOLD.bslq(:,1:2)];

lineNEW.uslq = [lineOLD.uslq(:,2:10)];
lineNEW.bslq = [lineOLD.bslq(:,2:10)];

aa = lineOLD.iusgq;
popoU = aa;
lala    = ~isspace(popoU); popoU = popoU(lala); 
%% >>>>>>>>>>> this is new, Feb 2, 2020 %%
if (lenID*lineOLD.linct ~= length(popoU))
  disp('>>>>> WARNING!!! translate2oldHITparams.m lenID*lineOLD.linct ~= length(popoU)')
  [mm,nn] = size(aa);
  for ii = 1 : mm
    junk = aa(ii,:);
    junk2 = strfind(junk,' ');
    numcount(ii) = length(junk2);
  end
  ratioA = numcount/length(junk);
  plot(1:15,numcount/length(junk),'bx-'); title(['iusgq GID ' num2str(lineOLD.igas(1))]);
  xlabel('character index'); ylabel('frac blanks'); grid
  hold on
  heldplot = +1;

  junk = [ lenID lineOLD.linct  lenID*lineOLD.linct  length(popoU)];
  fprintf(1,'maybe problem for HITEMP CO2??? lenID*lineOLD.linct  length(popoU %2i %10i %10i %10i \n',junk);

  iMethod = 2;
  if lineOLD.igas(1) == 2 & iMethod == 1
    fprintf(1,'old iusqg fix for GID %2i \n',lineOLD.igas(1))
    popoU = popoU(1:lenID*lineOLD.linct);
  elseif lineOLD.igas(1) == 2 & iMethod == 2
    fprintf(1,'new iusqg fix for GID %2i \n',lineOLD.igas(1))
    aa(7,:) = ' ';
    aa(9,:) = ' ';
    aa(11,:) = ' ';
    iiX = 08; junk = aa(iiX,:); junkB = strfind(junk,' '); fprintf(1,'iUSGQ %6i bad for index %2i \n',length(junkB),iiX); junk(junkB) = '9'; aa(iiX,:) = junk;
    iiX = 10; junk = aa(iiX,:); junkB = strfind(junk,' '); fprintf(1,'iUSGQ %6i bad for index %2i \n',length(junkB),iiX); junk(junkB) = '9'; aa(iiX,:) = junk;
    iiX = 12; junk = aa(iiX,:); junkB = strfind(junk,' '); fprintf(1,'iUSGQ %6i bad for index %2i \n',length(junkB),iiX); junk(junkB) = '9'; aa(iiX,:) = junk;
    iiX = 14; junk = aa(iiX,:); junkB = strfind(junk,' '); fprintf(1,'iUSGQ %6i bad for index %2i \n',length(junkB),iiX); junk(junkB) = '9'; aa(iiX,:) = junk;
    iiX = 15; junk = aa(iiX,:); junkB = strfind(junk,' '); fprintf(1,'iUSGQ %6i bad for index %2i \n',length(junkB),iiX); junk(junkB) = '9'; aa(iiX,:) = junk;
    popoU = aa;
    lala    = ~isspace(popoU); popoU = popoU(lala); 
  end
  %whos popoU aa
  %[lenID lineOLD.linct]
  %aa(:,1:3)
  %error(';lkgs')

  [mm,nn] = size(aa);
  for ii = 1 : mm
    junk = aa(ii,:);
    junk2 = strfind(junk,' ');
    numcount(ii) = length(junk2);
  end
  plot(1:15,numcount/length(junk),'ro-'); title(['iusgq GID ' num2str(lineOLD.igas(1))]);
  xlabel('character index'); ylabel('frac blanks'); grid
  hold off
  ratioB = numcount/length(junk);

  %[(1:15); ratioA; ratioB; ratioA-ratioB]'
end
%% >>>>>>>>>>> this is new, Feb 2, 2020 %%

popoU   = reshape(popoU,lenID,lineOLD.linct);
popoU   = popoU'; popoU = str2num(popoU);
popoUU  = ones(1,lineOLD.linct) * -1;    % to make sure all arrays have same size
foundUU = ones(1,lineOLD.linct) * -1;
for ii = 1 : length(xvibs)
  xx = xvibs(ii,2);
  kk = find(popoU == xx);
  popoUU(kk) = ii;
  foundUU(kk) = +1;
end

aa = lineOLD.ilsgq;
popoL = aa;
lala    = ~isspace(popoL); popoL = popoL(lala); 
%% >>>>>>>>>>> this is new, Feb 2, 2020 %%
if (lenID*lineOLD.linct ~= length(popoL))
  disp('>>>>> WARNING!!! translate2oldHITparams.m lenID*lineOLD.linct ~= length(popoL)')
  [mm,nn] = size(aa);
  for ii = 1 : mm
    junk = aa(ii,:);
    junk2 = strfind(junk,' ');
    numcount(ii) = length(junk2);
  end
  plot(1:15,numcount/length(junk),'x-'); title(['ilsgq GID ' num2str(lineOLD.igas(1))]);
  xlabel('character index'); ylabel('frac blanks'); grid

  junk = [ lenID lineOLD.linct  lenID*lineOLD.linct  length(popoL)];
  fprintf(1,'maybe problem for HITEMP CO2??? lenID*lineOLD.linct  length(popoL %2i %10i %10i %10i \n',junk);

  iMethod = 2;
  if lineOLD.igas(1) == 2 & iMethod == 1
    %% old
    popoL = popoL(1:lenID*lineOLD.linct);
  elseif lineOLD.igas(1) == 2 & iMethod == 2
    %% new
    aa(7,:) = ' ';
    aa(9,:) = ' ';
    aa(11,:) = ' ';
    iiX = 08; junk = aa(iiX,:); junkB = strfind(junk,' '); fprintf(1,'iLSGQ %6i bad for index %2i \n',length(junkB),iiX); junk(junkB) = '9'; aa(iiX,:) = junk;
    iiX = 10; junk = aa(iiX,:); junkB = strfind(junk,' '); fprintf(1,'iLSGQ %6i bad for index %2i \n',length(junkB),iiX); junk(junkB) = '9'; aa(iiX,:) = junk;
    iiX = 12; junk = aa(iiX,:); junkB = strfind(junk,' '); fprintf(1,'iLSGQ %6i bad for index %2i \n',length(junkB),iiX); junk(junkB) = '9'; aa(iiX,:) = junk;
    iiX = 14; junk = aa(iiX,:); junkB = strfind(junk,' '); fprintf(1,'iLSGQ %6i bad for index %2i \n',length(junkB),iiX); junk(junkB) = '9'; aa(iiX,:) = junk;
    iiX = 15; junk = aa(iiX,:); junkB = strfind(junk,' '); fprintf(1,'iLSGQ %6i bad for index %2i \n',length(junkB),iiX); junk(junkB) = '9'; aa(iiX,:) = junk;
    popoL = aa;
    lala    = ~isspace(popoL); popoL = popoL(lala); 
  end
  %whos popoL aa
  %[lenID lineOLD.linct]
  %aa(:,1:3)
end
%% >>>>>>>>>>> this is new, Feb 2, 2020 %%

popoL   = reshape(popoL,lenID,lineOLD.linct);
popoL   = popoL'; popoL = str2num(popoL);
popoLL  = ones(1,lineOLD.linct) * -1; % to make sure all arrays have same size
foundLL = ones(1,lineOLD.linct) * -1;
for ii = 1 : length(xvibs)
  xx = xvibs(ii,2);
  kk = find(popoL == xx);
  popoLL(kk) = ii;
  foundLL(kk) = +1;
end

lineNEW.iusgq = popoUU;
lineNEW.ilsgq = popoLL;

%lineOLD
%lineNEW

%semilogy(lineNEW.wnum,lineNEW.stren); 
if heldplot > 0
  hold off
end
