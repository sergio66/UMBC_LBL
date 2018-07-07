function lineNEW = should_I_translate2oldHITparams(line,gasID,hitran_version);

lineNEW = line;

bwah = 0;
if strfind(hitran_version,'h04')
  bwah = 1;
elseif strfind(hitran_version,'h08')
   bwah = 1;
elseif strfind(hitran_version,'h12') 
  bwah = 1; 
elseif strfind(hitran_version,'h16') 
  bwah = 1; 
end

if bwah > 0 & gasID == 2
  disp('for CO2, need to translate H04,H08,H12,H16 char(ilsgq,iusgq) to integer')
  lineNEW = translate2oldHITparams(line);
elseif bwah > 0 & gasID == 5
  disp('for CO, need to translate H04,H08,H12,H16 char(ilsgq,iusgq) to integer')
  lineNEW = translate2oldHITparams(line);
end
