function [iYes1,line1,iYes2,line2] = findlines_plot_compareGEISA_HITRAN(wv1,wv2,gid,GEISA,HITRAN);

addpath /home/sergio/MATLABCODE
figure(1)

if nargin < 4
  GEISA = 2015;
  HITRAN = 2016;
  HITRAN = 2020;
end
if nargin < 5
  HITRAN = 2016;
  HITRAN = 2020;
end

if GEISA ~= 2015
  error('only have GEISA 2015');
end

if mod(HITRAN,4) ~= 0
  error('HITRAN should be 2002,2004,....')
end  

[iYes1,line1] = findlines_plot(wv1,wv2,gid,GEISA,-1);
disp(' ')
[iYes2,line2] = findlines_plot(wv1,wv2,gid,HITRAN);

figure(1); clf; semilogy(line1.wnum,line1.stren,'b.'); title(['GEISA ' num2str(GEISA)]); 
  fprintf(1,'GEISA  %4i has %6i lines \n',GEISA,length(line1.wnum))
ax1 = axis;
figure(2); clf; semilogy(line2.wnum,line2.stren,'r.'); title(['HITRAN ' num2str(HITRAN)]); 
  fprintf(1,'HITRAN %4i has %6i lines \n',HITRAN,length(line2.wnum))
ax2 = axis;
figure(3); clf; semilogy(line1.wnum,line1.stren,'b.',line2.wnum,line2.stren,'r.');
  title('(b) GEISA (r) HITRAN')

figure(1); axis([ax1(1) ax1(2) min(ax1(3),ax2(3)) max(ax1(4),ax2(4))]); grid on
figure(2); axis([ax1(1) ax1(2) min(ax1(3),ax2(3)) max(ax1(4),ax2(4))]); grid on
figure(3); axis([ax1(1) ax1(2) min(ax1(3),ax2(3)) max(ax1(4),ax2(4))]); grid on

[I,i1,i2] = intersect(line1.wnum,line2.wnum);
oo1 = find(line1.iso == 1);
oo2 = find(line2.iso == 1);
for ii = 1 : length(oo1)
  bwop = abs(line1.wnum(oo1(ii))-line2.wnum(oo2));
  bwop = find(bwop == min(bwop),1);
  if abs(line1.wnum(oo1(ii))-line2.wnum(oo2(bwop))) < 0.005
    oo21(ii) = bwop;
  else
    oo21(ii) = NaN;
  end
end
wah = find(isfinite(oo21));
figure(1); plot(line1.wnum(oo1(wah)),line2.stren(oo2(oo21(wah)))./line1.stren(oo1(wah))); grid
    title(['STREN iso=1 ' num2str(HITRAN) '/' num2str(GEISA)]); 
figure(2); semilogy(line2.wnum(oo2(oo21(wah))),line2.stren(oo2(oo21(wah))),'rx',line1.wnum(oo1(wah)),line1.stren(oo1(wah)),'b.'); grid
    title(['STREN iso=1 (r)H' num2str(HITRAN) '(b)G' num2str(GEISA)]); 
disp('ret')

if length(I) > 0
  fprintf(1,'found %8i out of %8i common wavenumbers, plotting the strens and linewidths \n',length(I),length(line1.wnum))
  figure(3);
  semilogy(line2.wnum,line2.stren,'r.',line1.wnum,line1.stren,'bo'); grid
    title(['STREN (b)G' num2str(GEISA) ' (r)H' num2str(HITRAN)]);
  figure(4); clf
    wah = min(max(line1.iso),max(line2.iso));
    oo1 = find(line1.iso == 1);     oo2 = find(line2.iso == 1);
      semilogy(line2.wnum(oo2),line2.stren(oo2),'ro',line1.wnum(oo1),line1.stren(oo1),'m.')
      hold on
    oo1 = find(line1.iso == 2);     oo2 = find(line2.iso == 2);      
      semilogy(line2.wnum(oo2),line2.stren(oo2),'bo',line1.wnum(oo1),line1.stren(oo1),'c.')      
    oo1 = find(line1.iso == 3);     oo2 = find(line2.iso == 3);      
      semilogy(line2.wnum(oo2),line2.stren(oo2),'ko',line1.wnum(oo1),line1.stren(oo1),'g.')
      hold off
    str1a = ['H' num2str(HITRAN) ' iso1'];     str1b = ['G' num2str(GEISA) ' iso1'];
    str2a = ['H' num2str(HITRAN) ' iso2'];     str2b = ['G' num2str(GEISA) ' iso2'];
    str3a = ['H' num2str(HITRAN) ' iso3'];     str3b = ['G' num2str(GEISA) ' iso3'];    
    hl = legend(str1a,str1b,str2a,str2b,str3a,str3b,'location','best'); grid
    
  iDo = input('Do you want them plotted (-1) log(strn) or (+1) difference in stren (0) percent diff in stren ? ');
  if iDo < 0
    figure(4); clf; semilogy(line1.wnum(i1),line1.stren(i1),'b.',line2.wnum(i2),line2.stren(i2),'r+'); grid
    title(['STREN (b)G' num2str(GEISA) ' (r)H' num2str(HITRAN)]); 
    grid on
  elseif iDo > 0
    figure(4); clf; plot(line1.wnum(i1),line2.stren(i2) - line1.stren(i1),'k+'); grid
    title(['STREN G' num2str(GEISA) '-H' num2str(HITRAN)]); 
    grid on
  elseif iDo == 0
    figure(4); clf; plot(line1.wnum(i1),100*(line2.stren(i2) - line1.stren(i1))./line2.stren(i2),'k+'); grid
    title(['STREN percent diff ' num2str(HITRAN) '-' num2str(GEISA)]); 
    grid on
  end

  figure(5);
  semilogy(line2.wnum,line2.abroad,'r.',line1.wnum,line1.abroad,'bo'); grid
    title(['ABROAD (b)G' num2str(GEISA) ' (r)H' num2str(HITRAN)]); 
  iDo = input('Do you want them plotted (-1) broadening or (+1) difference in air broadening ? ');
  if iDo <= 0
    figure(6); clf; plot(line1.wnum(i1),line1.abroad(i1),'b.',line2.wnum(i2),line2.abroad(i2),'r+'); grid
    title(['ABROAD (b)G' num2str(GEISA) ' (r)H' num2str(HITRAN)]); 
    grid on
  elseif iDo > 0
    figure(6); clf; plot(line1.wnum(i1),line2.abroad(i2) - line1.abroad(i1),'k+'); grid
    title(['ABROAD ' num2str(HITRAN) '-' num2str(GEISA)]); 
    grid on
  end

  figure(7);
  semilogy(line2.wnum,line2.sbroad,'r.',line1.wnum,line1.sbroad,'bo'); grid
    title(['SBROAD (b)G' num2str(GEISA) ' (r)H' num2str(HITRAN)]); 
  iDo = input('Do you want them plotted (-1) broadening or (+1) difference in self broadening ? ');
  if iDo <= 0
    figure(8); clf; plot(line1.wnum(i1),line1.sbroad(i1),'b.',line2.wnum(i2),line2.sbroad(i2),'r+'); grid
    title(['SBROAD (b)G' num2str(GEISA) ' (r)H' num2str(HITRAN)]); 
    grid on
  elseif iDo > 0
    figure(8); clf; plot(line1.wnum(i1),line2.sbroad(i2) - line1.sbroad(i1),'k+'); grid
    title(['SBROAD ' num2str(HITRAN) '-' num2str(GEISA)]); 
    grid on
  end

  figure(9);
  semilogx(line1.wnum(i1),line1.abroad(i1)-line2.abroad(i2),'b.',line1.wnum(i1),line1.sbroad(i1)-line2.sbroad(i2),'r.'); grid
  semilogx(line1.stren(i1),line1.abroad(i1)-line2.abroad(i2),'b.',line1.stren(i1),line1.sbroad(i1)-line2.sbroad(i2),'r.'); grid
  xlabel('strength'); ylabel('\delta broad'); title('(b) air (r) self')

  %keyboard_nowindow
  %{
  iso1 = find(line1.iso == 1);
  iso2 = find(line2.iso == 1);
  junk2 = struct;
  for ii = 1 : length(iso1)
    boo = abs(line1.wnum(iso1(ii)) - line2.wnum(iso2));
    boo = find(boo == min(boo),1);
    junk2.iso(ii) = line2.iso(iso2(boo));
    junk2.wnum(ii) = line2.wnum(iso2(boo));
    junk2.stren(ii) = line2.stren(iso2(boo));
    junk2.abroad(ii) = line2.abroad(iso2(boo));    
    junk2.sbroad(ii) = line2.sbroad(iso2(boo));
    junk2.tsp(ii) = line2.tsp(iso2(boo));
  end
  figure(10); plot(line1.wnum(iso1)-junk2.wnum)
  figure(10); plot(line1.wnum(iso1),line1.stren(iso1)-junk2.stren)
  figure(10);
    semilogx(line1.stren(iso1),line1.abroad(iso1)-junk2.abroad,'b.',...
             line1.stren(iso1),line1.sbroad(iso1)-junk2.sbroad,'r.')
    xlabel('strength'); ylabel('\delta broad'); title('(b) air (r) self')  
  %}
  
  figure(4); figure(6); figure(8); figure(9);
end
