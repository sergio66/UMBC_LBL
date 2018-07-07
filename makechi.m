function chifcn = makechi(outwave,iFudgeType);

airs_chi_dir = '/asl/data/kcarta/KCARTADATA/General/ChiFile/co2_4um_fudge_';

if (iFudgeType == 2)
  chichunks = [2255 2280 2380 2405];
elseif (iFudgeType == 3)
  chichunks = [630 655 680 705 730 755 2180 2205 2230 2255 2280];
  chichunks = [chichunks 2355 2380 2405 2430 2530 2555 2580];
  end

%%this uses the CO2 chi functions that Scott has made from the AIRS radiances

chifcn = ones(size(outwave));

for ii = 1 : length(chichunks)
  if (iFudgeType == 2)
    file = [airs_chi_dir num2str(chichunks(ii)) '_b.txt'];
  elseif (iFudgeType == 3)
    file = [airs_chi_dir num2str(chichunks(ii)) '_c.txt'];
    end
  chis = load(file);
  plot(chis(:,1),chis(:,2)); pause(0.01);

  zz = interp1(chis(:,1),chis(:,2),outwave);
  jj = find(outwave < min(chis(:,1))); zz(jj) = 1;
  jj = find(outwave > max(chis(:,1))); zz(jj) = 1;

  jj = find(outwave >= min(chis(:,1)) & outwave <= max(chis(:,1)));
  chifcn(jj) = zz(jj);

  end

