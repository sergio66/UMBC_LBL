function [outwave,all_out_array] = run8watercontinuum_rtp(gasID,fmin,fmax,rtpfname,topts,iProfID,iSortPtotal_Hi2Lo);

%% WARNING : expects the wavenumber grid to be 0.0025 cm-1 when it pre-defines all_out_array
%% WARNING : expects the wavenumber grid to be 0.0025 cm-1 when it pre-defines all_out_array
%% WARNING : expects the wavenumber grid to be 0.0025 cm-1 when it pre-defines all_out_array

%% see run8watercontinuum.m : if topts.devide == -1, then you get the correct OD

[h,ha,p,pa] = rtpread(rtpfname);
if h.ptype ~= 1
  error('expect layers profiles')
end

addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/

if ~isfield(p,'plays')
  disp('adding p.plays')
  p = make_rtp_plays(p);
end

if nargin == 4
  topts = [];
  iProfID = 1 : length(p.stemp);
  iSortPtotal_Hi2Lo = -1;
elseif nargin == 5
  iProfID = 1 : length(p.stemp);
  iSortPtotal_Hi2Lo = -1;
elseif nargin == 6
  iSortPtotal_Hi2Lo = -1;
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

for ii = 1 : length(iProfID)
  figure(1); clf
  fprintf(1,'profile %5i of %5i \n',ii,length(iProfID));
  mr = mr_all(:,iProfID(ii));
  [hx,px] = subset_rtp_allcloudfields(h,p,[],[],iProfID(ii));
  fid = fopen('IPFILES/junk_wc.txt','w');
  fopen(fid);
  nn = px.nlevs-1;

  junk = [ones(nn,1) px.plays(1:nn)/1013.25 px.plays(1:nn).*mr(1:nn)/1013.25 px.ptemp(1:nn) px.gas_1(1:nn)/6.023e26];
  if iSortPtotal_Hi2Lo > 0
    %% see eg SPECTRA/ONEPROF_CODE/driver_rta.m
    boo = junk(:,2);
    [Y,I] = sort(boo,'descend');
    junk = junk(I,:);
  end

  fprintf(fid,'%3i %8.5e %8.5e %8.5f %8.5e \n',junk');
  fclose(fid);
  profname = 'IPFILES/junk_wc.txt';
  if length(topts) == 0
    [outwave,out_array] = run8watercontinuum(gasID,fmin,fmax,profname);
  else
    [outwave,out_array] = run8watercontinuum(gasID,fmin,fmax,profname,topts);
  end
  [mmm,nnn] = size(out_array);
  %whos outwave out_array
  woof = (1:nn) + (oomax-nn);  %% assumes lay 1 = gnd, lay nn = toa
  woof = (1:nn);               %% assumes lay 1 = toa, lay nn = gnd
  all_out_array(ii,woof,:) = out_array;
  figure(2); clf; semilogy(squeeze(all_out_array(:,[1 oomax-5],1)),'o-');
  pause(0.1)
end

if length(iProfID) == 1
  all_out_array = squeeze(all_out_array(1,:,:));
end  
