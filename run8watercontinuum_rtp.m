function [outwave,all_out_array] = run8watercontinuum_rtp(gasID,fmin,fmax,rtpfname,topts,iProfID);

%% WARNING : expects the wavenumber grid to be 0.0025 cm-1 when it pre-defines all_out_array
%% WARNING : expects the wavenumber grid to be 0.0025 cm-1 when it pre-defines all_out_array
%% WARNING : expects the wavenumber grid to be 0.0025 cm-1 when it pre-defines all_out_array


[h,ha,p,pa] = rtpread(rtpfname);
if h.pfields ~= 1
  error('expect layers profiles')
end

if nargin == 4
  topts = [];
  iProfID = 1 : length(p.stemp);
elseif nargin == 5
  iProfID = 1 : length(p.stemp);
end

if length(iProfID) > 1
  oomax = max(p.nlevs(iProfID)-1);
  fx = (fmax-fmin)/0.0025;
  all_out_array = zeros(length(iProfID),oomax,fx);
elseif length(iProfID) == 1
  oomax = max(p.nlevs(iProfID)-1);
  fx = (fmax-fmin)/0.0025;
  all_out_array = zeros(1,oomax,fx);
end

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(h,p,1:length(p.stemp),1);
mr_all = ppmvLAY/1e6;
[mm,nn] = size(mr_all);
mr_all(101,:) = 0;

% matrix   profname        this contains N x 5 matrix : layernum,p,pp,t,x
%                          where column 1 = layer number 
%                                column 2 = layer pressure (atm)
%                                column 3 = layer partial pressure (atm)
%                                column 4 = layer temperature (k)
%                                column 5 = layer gas amnt (kilomolecules/cm2)

topts

for ii = 1 : length(iProfID)
  figure(1); clf
  fprintf(1,'profile %5i of %5i \n',ii,length(iProfID));
  mr = mr_all(:,iProfID(ii));
  [hx,px] = subset_rtp(h,p,[],[],iProfID(ii));
  fid = fopen('IPFILES/junk_wc.txt','w');
  fopen(fid);
  nn = px.nlevs-1;
  junk = [ones(nn,1) px.plevs(1:nn)/1013.25 px.plevs(1:nn).*mr(1:nn)/1013.25 px.ptemp(1:nn) px.gas_1(1:nn)/6.023e26];
  fprintf(fid,'%3i %8.5e %8.5e %8.5f %8.5e \n',junk');
  fclose(fid);
  profname = 'IPFILES/junk_wc.txt';
  if length(topts) == 0
    [outwave,out_array] = run8watercontinuum(gasID,fmin,fmax,profname);
  else
    [outwave,out_array] = run8watercontinuum(gasID,fmin,fmax,profname,topts);
  end
  [mmm,nnn] = size(out_array);
  whos outwave out_array
  woof = (1:nn) + (oomax-nn);  %% assumes lay 1 = gnd, lay nn = toa
  woof = (1:nn);               %% assumes lay 1 = toa, lay nn = gnd
  all_out_array(ii,woof,:) = out_array;
  figure(2); clf; semilogy(squeeze(all_out_array(:,[1 oomax-5],1)),'o-');
  pause(0.1)
end

if length(iProfID) == 1
  all_out_array = squeeze(all_out_array(1,:,:));
end  
