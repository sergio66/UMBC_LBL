figure(1); clf
figure(2); clf

iaType = input('Enter which CKD versions to compare : [a b] : ');

i1 = 8; i2 = 8;
i1 = 4; i2 = 9;

for iBand = i1: i2
  if iBand == 1
    thedir = '/asl/data/kcarta/H2008_otherbands.ieee-le/FIR15_30/etc.ieee-le/';
  elseif iBand == 2
    thedir = '/asl/data/kcarta/H2008_otherbands.ieee-le/FIR30_50/etc.ieee-le/';
  elseif iBand == 3
    thedir = '/asl/data/kcarta/H2008_otherbands.ieee-le/FIR50_80/etc.ieee-le/';
  elseif iBand == 4
    thedir = '/asl/data/kcarta/H2004_otherbands.ieee-le/FIR80_150/etc.ieee-le/';
  elseif iBand == 5
    thedir = '/asl/data/kcarta/H2004_otherbands.ieee-le/FIR140_310/etc.ieee-le/';
  elseif iBand == 6
    thedir = '/asl/data/kcarta/H2004_otherbands.ieee-le/FIR300_510/etc.ieee-le/';
  elseif iBand == 7
    thedir = '/asl/data/kcarta/H2004_otherbands.ieee-le/FIR500_605/etc.ieee-le/';
  elseif iBand == 8
    thedir = '/asl/data/kcarta/H2008_otherbands.ieee-le/FIR605_2830/etc.ieee-le/';
    thedir = '/asl/data/kcarta/KCARTADATA/General/CKDieee_le/';
  elseif iBand == 9
    thedir = '/asl/data/kcarta/H2004_otherbands.ieee-le/NIR2830_3330/etc.ieee-le/';
  elseif iBand == 10
    thedir = '/asl/data/kcarta/H2004_otherbands.ieee-le/NIR3550_5550/etc.ieee-le/';
  elseif iBand == 11
    thedir = '/asl/data/kcarta/H2004_otherbands.ieee-le/NIR5550_8200/etc.ieee-le/';
  elseif iBand == 12
    thedir = '/asl/data/kcarta/H2004_otherbands.ieee-le/NIR8250_12250/etc.ieee-le/';
  elseif iBand == 13
    thedir = '/asl/data/kcarta/H2004_otherbands.ieee-le/VIS12000_25000/etc.ieee-le/';
  elseif iBand == 14
    thedir = '/asl/data/kcarta/H2004_otherbands.ieee-le/UV25000_45000/etc.ieee-le/';
  end

  fname2S = [thedir '/CKDSelf' num2str(iaType(2)) '.bin'];
  fname2F = [thedir '/CKDFor' num2str(iaType(2)) '.bin'];
  fname1S = [thedir '/CKDSelf' num2str(iaType(1)) '.bin'];
  fname1F = [thedir '/CKDFor' num2str(iaType(1)) '.bin'];

  if exist(fname2S) & exist(fname2F) & exist(fname1S) & exist(fname1F) 
    [ks2, freq2, temp2] = contread(fname2S);
    [kf2, freq2, temp2] = contread(fname2F);
    [ks1, freq1, temp1] = contread(fname1S);
    [kf1, freq1, temp1] = contread(fname1F);
    fprintf(1,'reading %s \n',fname2S);

    ii1 = find(temp1 == 300,1);
    ii2 = find(temp2 == 300,1);
    if (ii1 ~= ii2)
      error('ii1 ~= ii2')
    end
    ii = ii1;

    figure(1); hold on
    semilogy(freq1,ks1(ii,:),freq2,ks2(ii,:),'LineWidth',2);
    set(gca,'YScale','log');

    figure(2); hold on
    semilogy(freq1,kf1(ii,:),freq2,kf2(ii,:),'LineWidth',2);
    set(gca,'YScale','log');
  else
    fprintf(1,'oops Band %2i one or more files DNE \n',iBand)
  end

  pause
end

figure(1)
  legend(num2str(iaType(1)),num2str(iaType(2))); xlabel('Wavenumber cm-1'); ylabel('CS(300K)');
  set(gca,'YScale','log'); grid

figure(2)
  legend(num2str(iaType(1)),num2str(iaType(2))); xlabel('Wavenumber cm-1'); ylabel('CF(300K)');
  set(gca,'YScale','log'); grid