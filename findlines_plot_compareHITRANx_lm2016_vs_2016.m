function [iYes1,line1,iYes2,line2] = findlines_plot_compareHITRANx_lm2016_vs_2016(wv1,wv2,gid,HITRAN1,HITRAN2);

addpath /home/sergio/MATLABCODE

if nargin < 4
  %% looking at the linemix database
  HITRAN2 = 2017;  %% right now symbolic link to ../h16.by.gas/CO2/new_lm_g2.dat_Mar20_10.43am    which is based on orig LM linemix database downloaded from HITRAN in 2017
  HITRAN2 = 2018;  %% right now symbolic link to ../h16.by.gas/CO2/new_lm_g3_0.dat_Oct06_08.20am  which is based on new database sent to me by Iouli in Aug 2018
end

if HITRAN1 == HITRAN2
  error('HITRAN1 == HITRAN2')
end

[iYes1,line1] = findlines_plot(wv1,wv2,gid,HITRAN1);
disp(' ')
[iYes2,line2] = findlines_plot(wv1,wv2,gid,HITRAN2);

figure(1); clf; semilogy(line1.wnum,line1.stren,'b.'); title(num2str(HITRAN1)); 
  fprintf(1,'%4i has %6i lines \n',HITRAN1,length(line1.wnum))
ax1 = axis;
figure(2); clf; semilogy(line2.wnum,line2.stren,'r.'); title(num2str(HITRAN2)); 
  fprintf(1,'%4i has %6i lines \n',HITRAN2,length(line2.wnum))
ax2 = axis;

figure(1); axis([ax1(1) ax1(2) min(ax1(3),ax2(3)) max(ax1(4),ax2(4))]); grid on
figure(2); axis([ax1(1) ax1(2) min(ax1(3),ax2(3)) max(ax1(4),ax2(4))]); grid on

[I,i1,i2] = intersect(line1.wnum,line2.wnum);

if length(I) > 0
  fprintf(1,'found %8i out of %8i common wavenumbers, plotting the strens and linewidths \n',length(I),length(line1.wnum))
  figure(3);
  semilogy(line2.wnum,line2.stren,'r.',line1.wnum,line1.stren,'bo');
    title(['STREN (b) HITRAN2016 (r) LM2016']); 
  iDo = input('Do you want them plotted (-1) log(strn) or (+1) difference in stren (0) ratio of stren ? ');
  if iDo < 0
    figure(4); clf; semilogy(line1.wnum(i1),line1.stren(i1),'b.',line2.wnum(i2),line2.stren(i2),'r+')
    title(['STREN (b) HITRAN2016 (r) LM2016']); 
    grid on
  elseif iDo > 0
    figure(4); clf; plot(line1.wnum(i1),line1.stren(i1) - line2.stren(i2),'k+')
    title(['STREN HITRAN2016-LM2016']);         
    grid on
  elseif iDo == 0
    figure(4); clf; plot(line1.wnum(i1),line1.stren(i1) ./ line2.stren(i2),'k+')
    title(['STREN HITRAN2016/LM2016']);         
    grid on
  end

  figure(5);
  semilogy(line2.wnum,line2.abroad,'r.',line1.wnum,line1.abroad,'bo');
    title(['ABROAD (b) HITRAN2016 (r) LM2016']); 
  iDo = input('Do you want them plotted (-1) broadening or (+1) difference in air broadening ? ');
  if iDo < 0
    figure(6); clf; plot(line1.wnum(i1),line1.abroad(i1),'b.',line2.wnum(i2),line2.abroad(i2),'r+')
    title(['ABROAD (b) HITRAN2016 (r) LM2016']); 
    grid on
  elseif iDo > 0
    figure(6); clf; plot(line1.wnum(i1),line1.abroad(i1) - line2.abroad(i2),'k+')
    title(['ABROAD HITRAN2016-LM2016']); 
    grid on
  elseif iDo == 0
    figure(6); clf; plot(line1.wnum(i1),line1.abroad(i1) ./ line2.abroad(i2),'k+')
    title(['ABROAD HITRAN2016./LM2016']); 
    grid on
  end

  figure(7);
  semilogy(line2.wnum,line2.sbroad,'r.',line1.wnum,line1.sbroad,'bo');
    title(['SBROAD (b) HITRAN2016 (r) LM2016']); 
  iDo = input('Do you want them plotted (-1) broadening or (+1) difference in self broadening ? ');
  if iDo < 0
    figure(8); clf; plot(line1.wnum(i1),line1.sbroad(i1),'b.',line2.wnum(i2),line2.sbroad(i2),'r+')
    title(['SBROAD (b) HITRAN2016 (r) LM2016']); 
    grid on
  elseif iDo > 0
    figure(8); clf; plot(line1.wnum(i1),line1.sbroad(i1) - line2.sbroad(i2),'k+')
    title(['SBROAD HITRAN2016-LM2016']);     
    grid on
  elseif iDo == 0
    figure(8); clf; plot(line1.wnum(i1),line1.sbroad(i1) ./ line2.sbroad(i2),'k+')
    title(['SBROAD HITRAN2016/LM2016']);     
    grid on
  end

  figure(9);
  semilogx(line1.wnum(i1),line1.abroad(i1)-line2.abroad(i2),'b.',line1.wnum(i1),line1.sbroad(i1)-line2.sbroad(i2),'r.')
  semilogx(line1.stren(i1),line1.abroad(i1)-line2.abroad(i2),'b.',line1.stren(i1),line1.sbroad(i1)-line2.sbroad(i2),'r.')
  semilogx(line1.stren(i1),line1.abroad(i1)./line2.abroad(i2),'b.',line1.stren(i1),line1.sbroad(i1)./line2.sbroad(i2),'r.')    
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
