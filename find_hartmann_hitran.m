function [ggout,lenHitr,lenHart] = find_hartmann_hitran(lineIN,fhartdir,xx);

%% finds all lines in current HITRAN database that are in the HARTMANN band
%% line IN = HITRAN info
%% xx from /home/sergio/SPECTRA/JMHARTMANN/LM_PQR_CO2_2.0/Data/BandInfo.dat
%%    where xx(1),xx(2),xx(3) = iso iusgq ilsgq

%%simple
gg = find(lineIN.iso == xx(1) & lineIN.iusgq == xx(2) & lineIN.ilsgq == xx(3));

lenHitr = 0;
if length(gg) > 0
  ggbslq = lineIN.bslq(gg,6:8); ggbslqPQR = lineIN.bslq(gg,5);
  gg_lower = str2num(ggbslq);
  pp = find(ggbslqPQR == 'P'); gg_upper(pp) = gg_lower(pp) - 1; poink(pp) = -1;
  qq = find(ggbslqPQR == 'Q'); gg_upper(qq) = gg_lower(qq);     poink(qq) = 0;
  rr = find(ggbslqPQR == 'R'); gg_upper(rr) = gg_lower(rr) + 1; poink(rr) = +1;
  crap = lineIN.wnum(gg)/1000; crap = crap'; poink = poink'; 
  gg_upper = gg_upper'; gg_lower = gg_lower; 
  hitran_lines = [crap gg_lower gg_upper poink];
  [lenHitr,n1] = size(hitran_lines);
  [Y,I] = sort(hitran_lines(:,2)); 
  hitran_lines = ...
     [hitran_lines(I,1)  hitran_lines(I,2) hitran_lines(I,3) poink(I)]; 
end

%% but let us actually look at HARTMAN info for this band
strA = num2str(xx(1));
if xx(2) < 10
  strB = ['00' num2str(xx(2))];
elseif xx(2) < 100
  strB = ['0'  num2str(xx(2))];
else
  strB = [     num2str(xx(2))];
end

if xx(3) < 10
  strC = ['00' num2str(xx(3))];
elseif xx(3) < 100
  strC = ['0'  num2str(xx(3))];
else
  strC = [     num2str(xx(3))];
end
hartname = [fhartdir 'S' strA strB strC '.dat'];
hartxx = load(hartname);
hartmann_lines = [hartxx(:,1)/1000 hartxx(:,8:9)];
[Y,I] = sort(hartmann_lines(:,2)); 
hartmann_lines = ...
  [hartmann_lines(I,1)  hartmann_lines(I,2) hartmann_lines(I,3)]; 
[lenHart,n2]  = size(hartmann_lines);

if length(gg) > 0
  figure(2); clf;
  %% now find the lines in common
  lophit = [hitran_lines(:,2)   hitran_lines(:,3)];   [m1,n1] = size(lophit);
  lophat = [hartmann_lines(:,2) hartmann_lines(:,3)]; [m2,n2] = size(lophat);
  ggcount = 0; ggout = []; f0 = []; shift = [];
  for ll = 1 : m2
    poof = ones(m1,1)*lophat(ll,:) - lophit; 
    poof = abs(poof); poof = sum(poof');
    ix = find(poof == 0);
    if length(ix) == 1
      ggcount = ggcount + 1;
      ggout(ggcount) = gg(ix);
      f0(ggcount) = 1000*hartmann_lines(ll,1);
      fshift(ggcount) = 1000*(hartmann_lines(ll,1) - hitran_lines(ix,1));
    end
  end
  str = ['iso '  num2str(xx(1)) ' u,l = ' num2str(xx(2)) ' ' num2str(xx(3))];
  if length(f0) > 0
    plot(f0,fshift,'+'); title(str); pause(0.1)
  end
else
  %fprintf(1,'no coincident lines iso U L = %3i %3i %3i \n',xx(1),xx(2),xx(3));
  ggout = [];
end
figure(1)