function theans = quick_continuum_rtp_calc(fname,profile)

%% input : fname = rtp with layers profiles
%%         profile = which of the profiles to process
%% output : theans.freq = 05:1:2805
%%          theans.continuum : continuum at 100 layers
%% calls : /home/sergio/SPECTRA/tobin5.m

addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

%{
%% testing
fname = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49/regr49_1100_385ppm.op.rtp';
profile = 1;
theans = quick_continuum_rtp_calc(fname,profile);

topts.CKD = 6;
topts.ffin = 1;
profname = 'IPFILES/trp_water';
[outwave,out_array] = run8watercontinuum(1,605,2830,profname,topts);

boo = 1:400:890000;
plot(outwave(boo),flipud(out_array(:,boo)) ./ theans.continuum(boo:,1:2225));

figure(1); semilogy(outwave(boo),flipud(out_array(:,boo)),'b.', outwave(boo), theans.continuum(:,1:2225),'r');
           title('b : ckd6 r : tobin5')
figure(2); plot(outwave(boo),flipud(out_array(:,boo)) ./ theans.continuum(:,1:2225)); title('ckd6 / tobin5')
%}

[h,ha,p,pa] = rtpread(fname);
[h,p] = subset_rtp(h,p,[],[],profile);

temperature = p.ptemp;

if ~isfield(p,'plays')
  pN = p.plevs(1:end-1)-p.plevs(2:end);
  pD = log(p.plevs(1:end-1) ./ p.plevs(2:end));  
  p.plays = zeros(101,1);
  p.plays(2:101) = pN' ./ pD';
end

[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(h,p,1,1);
MR = [ppmvLAY/1e6; 0];

temp        = p.ptemp;
press       = p.plays/1013;
partpress   = p.plays .* MR/1013;
GasAmt      = p.gas_1/6.023e26;

freq = 605 : 2830;
selfmult = 1;
formult = 1;

array = [ones(101,1) press partpress temp GasAmt];
array = array(1:100,:);
theprofname = '/home/sergio/SPECTRA/IPFILES/water_mst';
theprofname = './water_mst';
fid=fopen(theprofname,'w');
fprintf(fid,'%4i    %12.8e  %12.8e  %12.8f  %12.8e \n',array');
fclose(fid);

%for ii = 1 : 100
%  continuum(ii,:) = tobin5(freq,temperature,press,partpress,GasAmt,selfmult,formult,ii);
%end

topts.CKD = 6;
topts.ffin = 1;
profname = './water_mst';
[freq, continuum] = run8watercontinuum(1,605,2830,profname,topts);

theans.freq      = freq;
theans.continuum = continuum;

semilogy(theans.freq,theans.continuum)