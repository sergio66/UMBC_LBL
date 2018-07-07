function [newline] = findUnionNew_oneLayer(...
                                 line,GasAmt,temperature,press,partpress,...
                                 A,B,C,D,G,hitran_version,mass_info,...
                                 gasID,mass_iso,optdepthMIN,isotope,TheLayers);
  
% same as findUnionNew except it is for one layer
% function [newline]=findUnion(line,GasAmt,Temp,Press,PartPress,...
%                              A,B,C,D,G,hitran_version,mass_info,...
%                              gasID,mass_iso,optdepthMIN,isotope,TheLayers);
% do a UNION of all lines so that lines do not turn on and off, turning  
% both SVD and Fast Models crazy 
% same as findUnion except is accounts for ABCDG not existing in later HITRAN
% and so uses info : hitran_version from hitread.m
%                    mass_info      from load_mass_dat.m

% also the code only selects the lines belonging to the chosen isotopes!

line = should_I_translate2oldHITparams(line,gasID,hitran_version); 

numlines=line.linct;

if (isotope ~= 0)  %have to choose the lines belonging to particular isotope
  fprintf(1,'\n findUnionNew_oneLayer : all isotopes have %6i lines',numlines);
  number_of_isotopes = unique(line.iso);
  number_of_isotopes = max(line.iso);
  if (number_of_isotopes < isotope)
    error('you have selected nonexisting isotope!');
  end
  theindex = [];
  for oo = 1 : length(isotope)
    theindexx = find(line.iso == isotope(oo));
    wah = [gasID isotope(oo)  length(theindexx)];
    fprintf(1,'  gasID   isotope   numlines = %3i %3i %6i \n',wah);    
    theindex = [theindex theindexx];
  end
  if (length(theindex) == 0) 
    error('sorry! this isotope(s) has no lines in freq region chosen!');
  else
    line.linct  = length(theindex); 
    line.wnum   = line.wnum(theindex); 
    line.stren  = line.stren(theindex); 
    line.tprob  = line.tprob(theindex); 
    line.abroad = line.abroad(theindex); 
    line.sbroad = line.sbroad(theindex); 
    line.els    = line.els(theindex); 
    line.abcoef = line.abcoef(theindex); 
    line.tsp    = line.tsp(theindex); 
    line.iusgq  = line.iusgq(theindex); 
    line.ilsgq  = line.ilsgq(theindex); 
    line.iso    = line.iso(theindex);  
    line.gasid  = line.gasid(theindex);  
    numlines    = line.linct;
    line.uslq   = line.uslq(theindex,:); 
    line.bslq   = line.bslq(theindex,:); 
  end
  fprintf(1,'\n chosen isotope, number of lines = %6i',numlines);
end

if (line.linct > 0)
  fprintf(1,' >> number of lines = %6i \n',line.linct);
  fprintf(1,' >> min/max wavennumber of lines = %8.6f %8.6f \n',min(line.wnum),max(line.wnum));
  scatter(line.wnum,log10(line.stren),20,line.iso,'filled'); colorbar
  grid; xlabel('wavenumber cm-1'); ylabel('log10(stren)'); title('isotope')
else
  fprintf(1,'sorry : no lines found!!!! \n');
end

%unique(line.iso)
%whos mass_iso 
%line

%% hmm very oldly, esp for CO2 in H2012, we seem to get sotope 0
%% 20 4599.238816 4.171E-30 9.805E-02.06810.081  568.52390.76-.005947
boo = find(line.iso == 0);
if length(boo > 0)
  fprintf(1,'WOW WOW : have isotope 0 for gas %2i !!! ... removing those %6i out of %6i lines \n',...
          gasID,length(boo),length(line.wnum))
  disp('you may get a "WARNING .... findUnionNew somehow missed some lines (strengthM == 0)" msg')
  disp('can ignore that ...')
  theindex = find(line.iso > 0);
  line.linct  = length(theindex); 
  line.wnum   = line.wnum(theindex); 
  line.stren  = line.stren(theindex); 
  line.tprob  = line.tprob(theindex); 
  line.abroad = line.abroad(theindex); 
  line.sbroad = line.sbroad(theindex); 
  line.els    = line.els(theindex); 
  line.abcoef = line.abcoef(theindex); 
  line.tsp    = line.tsp(theindex); 
  line.iusgq  = line.iusgq(theindex); 
  line.ilsgq  = line.ilsgq(theindex); 
  line.iso    = line.iso(theindex);  
  line.gasid  = line.gasid(theindex);  
  numlines    = line.linct;
  line.uslq   = line.uslq(theindex,:); 
  line.bslq   = line.bslq(theindex,:); 
  line.ai     = line.ai(theindex,:); 
end

mass = mass_iso(line.iso); 
mass = mass';

newline = line;

theindex = 1 : length(line.wnum);
  %now set "newline" with all these chosen lines
  newline.linct  = length(theindex); 
  newline.wnum   = line.wnum(theindex); 
  newline.stren  = line.stren(theindex); 
  newline.tprob  = line.tprob(theindex); 
  newline.abroad = line.abroad(theindex); 
  newline.sbroad = line.sbroad(theindex); 
  newline.els    = line.els(theindex); 
  newline.abcoef = line.abcoef(theindex); 
  newline.tsp    = line.tsp(theindex); 
  if gasID == 2
    %% have done the translation from 15xN char to Nx1 integrers
    newline.iusgq  = line.iusgq(theindex);    %% before Feb 2014
    newline.ilsgq  = line.ilsgq(theindex);    %% before Feb 2014
  else
    %% keep as 15xN char
    newline.iusgq  = line.iusgq(:,theindex);    %% after Feb 2014
    newline.ilsgq  = line.ilsgq(:,theindex);    %% after Feb 2014
  end
  newline.iso    = line.iso(theindex); 
  newline.gasid  = line.gasid(theindex); 
  newline.uslq   = line.uslq(theindex,:); 
  newline.bslq   = line.bslq(theindex,:); 
  newline.ref    = line.ref(theindex,:); 

newline
