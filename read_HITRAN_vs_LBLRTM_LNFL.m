function [linesU,linesL] = read_HITRAN_vs_LBLRTM_LNFL(f1,f2,gid);

iVersL = 12.4;
iVersL = 12.8;
[iYesL,linesL] = read_LBLRTM_LNFL(f1,f2,gid,iVersL);

iVersU = 2012;
iVersU = 2016;
[iYesU,linesU] = findlines_plot(f1,f2,gid,iVersU);

figure(1); semilogy(linesU.wnum,linesU.stren,'b.',linesL.wnum,linesL.stren,'r.'); title('b : HITRAN r:AER');

isoU = unique(linesU.iso);
fprintf(1,'HITRAN : found %6i total lines %2i isotopes \n',length(linesU.iso),length(isoU));
for ii = 1 : length(isoU)
  bahU = find(linesU.iso == ii);
  fprintf(1,'     iso %2i has %6i lines \n',ii,length(bahU));
end

disp(' ')
isoL = unique(linesL.iso);
fprintf(1,'LBLRTM : found %6i total lines %2i isotopes \n',length(linesL.iso),length(isoL));
for ii = 1 : length(isoL)
  bahL = find(linesL.iso == ii);
  fprintf(1,'     iso %2i has %6i lines \n',ii,length(bahL));
end

ii=1;
  bahU = find(linesU.iso == ii);
  bahL = find(linesL.iso == ii);
  figure(2); semilogy(linesU.wnum(bahU),linesU.stren(bahU),'b.',linesL.wnum(bahL),linesL.stren(bahL),'r.'); title('b : HITRAN r:AER ISO = 1');

  [Y,IU,IL] = intersect(linesU.wnum(bahU),linesL.wnum(bahL));
  fprintf(1,'  found %5i common ISO=1 lines \n',length(Y));
  figure(3); plot(linesU.wnum(bahU(IU)),linesU.stren(bahU(IU)) ./ linesL.stren(bahL(IL))','r.'); title('HITRAN/AER ISO = 1');  

ix = 0;
for ff = f1 : 0.0625 : f2
  ix = ix + 1;
  fjunk(ix) = ff;
  
  boo = find(linesU.wnum >= ff-0.25 & linesU.wnum <= ff+0.25);
  if length(boo) > 0
    maxU(ix) = max(linesU.stren(boo));
    xah = find(linesU.stren(boo) > 1e-26);    
    meanU(ix) = mean(linesU.stren(boo(xah)));
  else
    maxU(ix) = NaN;;
    meanU(ix) = NaN;
    
  end

  boo = find(linesL.wnum >= ff-0.25 & linesL.wnum <= ff+0.25);
  if length(boo) > 0
    maxL(ix) = max(linesL.stren(boo));
    xah = find(linesL.stren(boo) > 1e-26);        
    meanL(ix) = mean(linesL.stren(boo(xah)));        
  else
    maxL(ix) = NaN;;
    meanL(ix) = NaN;    
  end
end
figure(4); semilogy(fjunk,maxU,'b.',fjunk,maxL,'r',fjunk,meanU,'c.',fjunk,meanL,'m');
  hl = legend('max HITRAN','max LBLRTM','mean HITRAN','mean LBLRTM');
  
figure(5); semilogy(fjunk,maxU ./ maxL,'r', fjunk,meanU ./ meanL,'b'); title('(r) maxHITRAN/maxAER (b)meanHITRAN/meanAER'); 

monk = maxU ./ maxL;   monk(isnan(monk)) = 1;
bonk = meanU ./ meanL; bonk(isnan(bonk)) = 1;
figure(6); plot(fjunk,sgolayfilt(monk,3,101),'r',fjunk,sgolayfilt(bonk,3,101),'b','linewidth',2); title('smoothed (r) maxHITRAN/maxAER (b)meanHITRAN/meanAER'); 
axis([f1 f2 0.8 1.2]); grid