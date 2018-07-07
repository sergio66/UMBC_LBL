%% this is us
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf6.bin'; [kx6s, freq, temp] = contread(fname);
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor6.bin';  [kx6f, freq, temp] = contread(fname);

%% this is from LBLRTM in around 2012
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf25.bin'; [kx25s, freq, temp] = contread(fname);
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor25.bin';  [kx25f, freq, temp] = contread(fname);

%% this is from LBLRTM in around 2016
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDSelf27.bin'; [kx27s, freq, temp] = contread(fname);
fname = '/asl/data/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor27.bin';  [kx27f, freq, temp] = contread(fname);

figure(1); semilogy(freq,kx25s ./ kx6s); title('self 25/6')
figure(2); semilogy(freq,kx25f ./ kx6f); title('forn 25/6')
figure(3); semilogy(freq,kx27s ./ kx6s); title('self 25/7')
figure(4); semilogy(freq,kx27f ./ kx6f); title('forn 25/7')