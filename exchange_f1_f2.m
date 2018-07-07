function [freq,str] = exchange_f1_f2(line,ind,iso,v_u,v_l);

%%    where xx(1),xx(2),xx(3) = iso iusgq ilsgq
xx = [iso v_u v_l];

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

str = [strA strB strC];
fread = ['CO2_MATFILES/EXCHANGE/hartmanS' strA strB strC '.mat'];
ee = exist(fread,'file');
if ee > 0
  hartxx = load(fread); 
  if hartxx.iusgq ~= xx(2)
    fprintf(1,'%s iusgq does not match %3i %3i \n',fread,xx(2),hartxx.iusgq);
    error('ugh');
    end 
  if hartxx.ilsgq ~= xx(3)
    fprintf(1,'%s ilsgq does not match %3i %3i \n',fread,xx(3),hartxx.ilsgq);
    error('ugh');
    end 
  if hartxx.iso ~= xx(1)
    fprintf(1,'%s iso does not match %3i %3i \n',fread,xx(1),hartxx.iso);
    error('ugh');
    end 
  poink = zeros(size(hartxx.f_hartmann));
    pp = find(hartxx.j_upper-hartxx.j_lower == -1); poink(pp) = -1;
    qq = find(hartxx.j_upper-hartxx.j_lower ==  0); poink(qq) =  0;
    rr = find(hartxx.j_upper-hartxx.j_lower == +1); poink(rr) = +1;
  hartmann_lines = [hartxx.f_hartmann; hartxx.j_lower; hartxx.j_upper; poink]';
  [lenHart,n2]  = size(hartmann_lines);

  gg = ind;
  ggbslq = line.bslq(gg,6:8); ggbslqPQR = line.bslq(gg,5);
  gg_lower = str2num(ggbslq);
  pp = find(ggbslqPQR == 'P'); gg_upper(pp) = gg_lower(pp) - 1; poink(pp) = -1;
  qq = find(ggbslqPQR == 'Q'); gg_upper(qq) = gg_lower(qq);     poink(qq) = 0;
  rr = find(ggbslqPQR == 'R'); gg_upper(rr) = gg_lower(rr) + 1; poink(rr) = +1;
  crap = line.wnum(gg); crap = crap'; poink = poink'; 
  gg_upper = gg_upper'; gg_lower = gg_lower; 
  hitran_lines = [crap gg_lower gg_upper poink];
  [lenHitr,n1] = size(hitran_lines); 

  linejj = line.bslq(ind,:);
  linePQR = linejj(:,5);
  linejj  = num2str(linejj(:,6:8));
  for ii = 1 : length(ind)
    iix = find(hartxx.j_lower == hitran_lines(ii,2) & ...
               hartxx.j_upper == hitran_lines(ii,3));
    if length(iix) == 1
      freq(ii) = hartxx.f_hartmann(iix);
    else
      freq(ii) = line.wnum(ind(ii));
      end
    end
else
  freq = line.wnum(ind);
  end
freq = freq';