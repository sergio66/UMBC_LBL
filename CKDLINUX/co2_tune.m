%%%this is for subroutine Co2_4um_fudge in kcoeffMAIN.f of kCARTA

data=load('tunmlt_jan04deliv.dat');

kc_chunks = 605:25:2805;

kc_freqs = 605:0.0025:2830-0.0025;

tun_fr  = data(:,2);
tun_co2 = data(:,3);
[B,I,J] = unique(tun_fr);
tun_fr  = tun_fr(I);
tun_co2 = tun_co2(I);

kc_co2 = interp1(tun_fr,tun_co2,kc_freqs);

kc_co2_new0 = kc_co2;
ii = isnan(kc_co2);
kc_co2_new0(ii) = 1.0;
ii = find(kc_freqs >= 1613 & kc_freqs < 2191.5);
kc_co2_new0(ii) = 1.0;
ii = find(kc_freqs > 2589);
kc_co2_new0(ii) = 1.0;
plot(kc_freqs,[kc_co2; kc_co2_new0]);
kc_co2_new = kc_co2_new0;

kchifile = '/asl/data/kcarta/KCARTADATA/General/ChiFile/';
kc_current = [2255 2280 2380 2405];

newchunk = 0;
for ii = 1 : length(kc_chunks)-1
  jj = find(kc_freqs >= kc_chunks(ii) & ...
            kc_freqs <= kc_chunks(ii+1)-0.0025+0.0005);

  current_int = intersect(kc_chunks(ii),kc_current);
  if length(current_int) > 0
    disp('loading file ....')
    fname = [kchifile 'co2_4um_fudge_' num2str(kc_chunks(ii)) '_b.txt'];
    current_chi = load(fname);
    kc_co2_new(jj) = kc_co2_new0(jj).*current_chi(:,2)';
    plot(kc_freqs(jj),kc_co2_new0(jj),current_chi(:,1),current_chi(:,2),...
         kc_freqs(jj),kc_co2_new(jj)); pause;
    end

  mean_mult = mean(kc_co2_new(jj));
  chipoints = find(kc_co2_new(jj) ~= 1);

  if (mean_mult ~= 1.0)
    newchunk = newchunk + 1;
    fprintf(1,'freq ii newchunk %5i %3i %3i %3i \n',floor(kc_chunks(ii)),...
            ii,newchunk,length(chipoints));
    kc_new(newchunk) = kc_chunks(ii);
    plot(kc_freqs(jj),kc_co2_new(jj)); axis tight
    str = ['chunk = ' num2str(ii) ' start freq = ' num2str(kc_chunks(ii))];
    str = [str ' mean chi = ' num2str(mean_mult)]; 
    title(str); pause(0.5);

    if length(jj) ~= 10000
      error('ooooh!')
      end

    do_output = -1;
    if do_output > 0
      outputstuff = [kc_freqs(jj); kc_co2_new(jj)];
      fname = [kchifile 'co2_4um_fudge_' num2str(kc_chunks(ii)) '_c.txt'];
      fid = fopen(fname,'w');
      fprintf(fid,'  %12.6f  %12.6f \n',outputstuff);
      fclose(fid);
      end

  else
    fprintf(1,'freq ii newchunk %5i %3i \n',floor(kc_chunks(ii)),ii);
    end
  end