%%%this is a script to run the line-by-line code (with nothing fancy)

outdir = '/carrot/s1/sergio/PAGANO_CH4/';      %%results go here

gasIDs = [1 2 3 6 10 13 15 20 23 24 26];   %% all gases in this region

chunks = 2905:25:2905;
chunks = 2930:25:2930;
chunks = 2880:25:2930;

for fx = 1 : length(chunks)
  fmin = chunks(fx);
  fmax = fmin + 25;

  outSum = zeros(100,10000);
  outLowestLayer = zeros(length(gasIDs),10000);

  for iG = 1 : length(gasIDs)
    gasID = gasIDs(iG);
    profname = ['/home/sergio/SPECTRA/IPFILES/PAGANO/std_'  num2str(gasID)];
    clear topts

    if gasID ~= 1
      [outwave,out_array] = run7(gasID,fmin,fmax,profname);
    else
      [outwave,out_arrayWB] = run7water(gasID,fmin,fmax,profname);
      topts.CKD = 1;
      [outwave,out_arrayC]= run7watercontinuum(gasID,fmin,fmax,profname,topts);
      out_array = out_arrayWB + out_arrayC;
      clear out_arrayWB out_arrayC
      end

    outLowestLayer(iG,:) = out_array(1,:);
    outSum = outSum + out_array;

    %%%the just in case saver
    saver = ['save ' outdir '/indgas' num2str(gasID) '_'];
    saver = [saver num2str(fmin) '.mat out_array outwave'];
    eval(saver);
    
    end
  saver = ['save ' outdir '/sumgases' num2str(length(gasIDs)) '_'];
  saver = [saver num2str(fmin) '.mat outLowestLayer outwave outSum'];
  eval(saver);
  end

plot(outwave,outLowestLayer)
semilogy(outwave,outLowestLayer,outwave,outSum(1,:),'k.-') 

