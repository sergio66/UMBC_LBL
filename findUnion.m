function [newline]=findUnion(line,GasAmt,temperature,press,partpress,...
                             A,B,C,D,G,...
                             gasID,mass_iso,optdepthMIN,isotope,TheLayers);

% function [newline]=findUnion(line,GasAmt,Temp,Press,PartPress,...
%                              A,B,C,D,G,...
%                              gasID,mass_iso,optdepthMIN,isotope,TheLayers);
% do a UNION of all lines so that lines do not turn on and off, turning  
% both SVD and Fast Models crazy 

% also the code only selects the lines belonging to the chosen isotopes!

numlines=line.linct;

if (isotope ~= 0)  %have to choose the lines belonging to particular isotope
  fprintf(1,'\n all isotopes together, number of lines = %6i',numlines);
  number_of_isotopes = max(line.iso);
  if (number_of_isotopes < isotope)
    error('you have selected nonexisting isotope!');
    end
  theindex=find(line.iso == isotope);
  if (length(theindex) == 0) 
    error('sorry! this isotope has no lines in freq region chosen!');
  else
    line.linct   = length(theindex); 
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
    numlines=line.linct;
    end
  fprintf(1,'\n chosen isotope, number of lines = %6i',numlines);
  end

if (line.linct > 0)
  plot(line.iso,'+')
else
  fprintf(1,'sorry : no lines found!!!! \n');
  end

mass=mass_iso(line.iso); 
mass=mass';

%%%%%%%%%%% only do layer with MAX gas amount
amtMax = max(GasAmt);
amtMin = min(GasAmt);
jj     = find(abs((GasAmt-amtMax)/amtMin) <= 10*eps);
if (length(jj) > 1)
  jj = jj(1);
  end
fprintf(1,'max gas amount is in layer %3i \n',jj);

tempr = temperature(jj); 
p     = press(jj); 
ps    = partpress(jj); 
amt   = GasAmt(jj);

if (line.linct > 0)
  qfcn     = q(A,B,C,D,G,line,tempr);
  pwr      = line.abcoef; 
  for_brd  = line.abroad; 
  self_brd = line.sbroad; 
  freq     = line.wnum+press(jj)*line.tsp;  %freq shift 
  energy   = line.els; 
  s0       = line.stren; 
  brd      = broad(p,ps,1.0,for_brd,self_brd,pwr,tempr,gasID); 
  strength = find_stren(qfcn,freq,tempr,energy,s0,amt); 

  maxVHH   = vhhOPTDEP(freq,tempr,mass,brd);  %max of VHH
  %optical depth > some value 
  %tempUNIONold=find(strength(length(mass))*maxVHH > optdepthMIN)
  tempUNION = find(strength.*maxVHH > optdepthMIN);
  plot(freq,strength.*maxVHH); %%pause(0.1)
  theindex = tempUNION;
  len0     = length(theindex);

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
  newline.iusgq  = line.iusgq(theindex); 
  newline.ilsgq  = line.ilsgq(theindex); 
  newline.iso    = line.iso(theindex); 

  if (numlines ~= length(theindex))
    fprintf(1,'\n in findUnion, making sure all levels have same lines \n');
    fprintf(1,'\n orig num = %5i new num = %5i \n',numlines,length(theindex));
    end
  end

