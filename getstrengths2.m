fnamePRE = '/asl/data/hitran/h98.by.gas/g';
do_HITRAN_vers;
fnamePRE = [HITRAN '/g'];; 

fnamePOST = '.dat';                        

maxii = 13;

for ii = 1:maxii
  gasID = ii;
  fnameIN = int2str(gasID);
  fname = [fnamePRE fnameIN fnamePOST];
  [line] = hitread(400,500,1e-28,gasID,fname);
  maxstren(ii) = 0.0;
  if line.linct > 1
    tt = num2str(ii);
    eval(['load RefProf/refgas' tt]);
    eval(['prof = refgas' tt ';']);
    maxq = max(prof(:,5));
    plot(line.wnum,line.stren); title(tt); pause(5);
    maxstren(ii) = max(line.stren*maxq*6.02e26);
    end
  end

semilogy(1:maxii,maxstren,'+'); grid on
