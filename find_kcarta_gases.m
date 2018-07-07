fx = input('Enter start stop wavenumbers (in 25cm-1 chunks) : ');
f1 = fx(1); f2 = fx(2);
fdat = f1 : 25 : f2;
data1 = load('/home/sergio/KCARTA/INCLUDE/comp107.param');
data2 = load('/home/sergio/KCARTA/INCLUDE/xsec107.param');

data = [data1; data2];

gases = [];
for jj = 1 : length(fdat)-1
  f1a = fdat(jj);
  f2a = f1a + 25;
  for ii = 1 : length(data)
    if data(ii,2) <= f1a & data(ii,3) >= f2a
      gases = [gases data(ii,1)];
      end
    end
  end
gasesX = unique(gases);

for jj = 1 : length(gasesX)
  doinkX(jj) = length(find(gases == gasesX(jj)));
  end

disp('gasID num of occurences');
disp('-----------------------');
[gasesX; doinkX]'