function [outwave,all_out_array] = run8_rtp(gasID,fmin,fmax,rtpfname,topts,iProfID)

if isfield(topts,'tempname')
  tempname = topts.tempname;
  topts = rmfield(topts,'tempname');
else
  tempname = [];
end

outwave = [];
all_out_array = [];

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
[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(h,p,1:length(p.stemp),gasID);
mr_all = ppmvLAY/1e6;
[mm,nn] = size(mr_all);
mr_all(101,:) = 0;
plot(mr_all); title(['MR for gas ' num2str(gasID)]);

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
  plot(mr,1:101); ylim([1 101]); set(gca,'ydir','reverse'); title(['MR for gas ' num2str(gasID) ' profile ' num2str(iProfID(ii))]); pause(0.1); 

  [hx,px] = subset_rtp_allcloudfields(h,p,[],[],iProfID(ii));
  nn = px.nlevs-1;
  str = ['junkQ = px.gas_' num2str(gasID) '(1:nn);'];
  eval(str);
  fprintf(1,'max gas amount for gasID %2i is %12.6e molecules/cm2 \n',gasID,max(max(junkQ)));

  if length(tempname) == 0
    fid = fopen('IPFILES/junk_wc.txt','w');
    profname = 'IPFILES/junk_wc.txt';
  else
    fid = fopen(tempname,'w');
    profname = tempname;
  end
  fopen(fid);
  junk = [ones(nn,1) px.plays(1:nn)/1013.25 px.plays(1:nn).*mr(1:nn)/1013.25 px.ptemp(1:nn) junkQ/6.023e26];
  fprintf(fid,'%3i %8.5e %8.5e %8.5f %8.5e \n',junk');
  fclose(fid);

  if gasID == 1
    if length(topts) == 0
      [outwave,out_array] = run8water(gasID,fmin,fmax,profname);
      [outwave,cout_array] = run8watercontinuum(gasID,fmin,fmax,profname);
    else
      [outwave,out_array] = run8water(gasID,fmin,fmax,profname,topts);
      [outwave,cout_array] = run8watercontinuum(gasID,fmin,fmax,profname,topts);
    end
    out_array = out_array + cout_array;

  %elseif gasID == 2
  %  if length(topts) == 0
  %    [outwave,out_array] = run8co2(gasID,fmin,fmax,profname);
  %  else
  %    [outwave,out_array] = run8co2(gasID,fmin,fmax,profname,topts);
  %  end

  elseif gasID <= 40
    if gasID == 2
      disp('WARNING : should really be doing CO2 linemixing, but just doing usual voigt')
    end
    if length(topts) == 0
      [outwave,out_array] = run8(gasID,fmin,fmax,profname);
    else
      [outwave,out_array] = run8(gasID,fmin,fmax,profname,topts);
    end
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
