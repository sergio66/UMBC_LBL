function  lineNEW = translate2oldHITparams(lineOLD);

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