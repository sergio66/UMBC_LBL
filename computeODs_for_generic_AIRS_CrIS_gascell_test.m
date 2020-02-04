%% 1 atm = 760 torr

% pV = nRT ==> q = (n/V)L = p/RT 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define file with gas cell params %% three coluns,Nrows : colums are [gasID ps(torr) T(degC)]

paramfile = 'params_jan27_2020';    %% name of file with gas cell params ... note that
                                    %% (a) expectign units of torr and degC
                                    %% (b) there is no .txt or .whatever at the end of the file name
                                    %%     since this code appends _ip and .mat to output files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN CODE

addpath /home/sergio/SPECTRA

gascellparams = load([paramfile]);  %% three coluns,Nrows : colums are [gasID ps(torr) T(degC)]
[mmg,nng] = size(gascellparams);
if nng ~= 3
  error('need 3 columns');
end

homedir = pwd;

%for ii = 1 : mmg
for ii = 1 : 3

  %% see /home/sergio/SPECTRA/gascell2ip.m
  %% p  = 5.93000000000000e-02
  %% ps = 5.93000000000000e-02
  %% T = 288.1500;
  %% L = 1;

  gid = gascellparams(ii,1);
  q = gascell2ip(gascellparams(ii,2)/760.0,gascellparams(ii,2)/760.0,gascellparams(ii,3)+273.15,1.0);

  % matrix   profname       this contains N x 5 matrix : layernum,p,pp,t,x
  %                          where column 1 = layer number
  %                                column 2 = layer pressure (atm)
  %                                column 3 = layer partial pressure (atm)
  %                                column 4 = layer temperature (k)
  %                                column 5 = layer gas amt (kilomoles/cm2)

  ipfile = [homedir '/' 'gascell_' num2str(ii) '_gid_' num2str(gid) '_' paramfile '_ip'];
  fid = fopen(ipfile,'w');
  array = [gid gascellparams(ii,2)/760.0 gascellparams(ii,2)/760.0 gascellparams(ii,3)+273.15 q];
  fprintf(fid,'%3i %8.6e %8.6e %8.6f %8.6e \n',array);
  fclose(fid);

  fprintf(1,'doing gascell %3i of %3i gid = %2i ipfile = %s \n',ii,mmg,gid,ipfile)

  if (gid == 2)
    f1 = 605;    f2 = 1105;
    f1 = 605;    f2 = 1230;  f2 = 1305;

    cd  /home/sergio/SPECTRA
    topts = struct;
    topts.mainloop = +1; 
    %topts.mainloop = -1; 

    tic
    for ii = 1 : floor((f2-f1)/100)
      f1x = f1 + (ii-1)*100;
      f2x = f1x + 100;
      ind = (1:40000) + (ii-1)*40000;
      [xoutwave,xout_array] = run8co2(2,f1x,f2x,ipfile,topts);   %% need length(w) <= kmaxlen = 200000
      outwave(ind) = xoutwave;
      out_array(ind) = xout_array;
    end
    toc
    
    cd  /home/sergio/HITRAN2UMBCLBL/LBLRTM
    [xw,xdglab,xdlblrtm] = driver_glab_lblrtm_self_MANYLAY(2,f1,f2,ipfile,-1);
    semilogy(outwave,out_array,'b',xw,xdlblrtm,'r')
    plot(outwave,out_array./xdlblrtm','r')
    
    cd ~/MATLABCODE/QUICKTASKS_TELECON/CrIS_gascell_Jan2020/
    saver = ['save ' ipfile '.mat out_array outwave xw xdglab xdlblrtm'];
    eval(saver)

  elseif (gid == 5)
    f1 = 1780;    f2 = 2280;
    f1 = 2105;    f2 = 2605;

    cd /home/sergio/SPECTRA
    [outwave,out_array] = run8(5,f1,f2,ipfile);

    cd  /home/sergio/HITRAN2UMBCLBL/LBLRTM
    [xw,xdglab,xdlblrtm] = driver_glab_lblrtm_self_MANYLAY(5,f1,f2,ipfile,-1);
    semilogy(outwave,out_array,xw,xdlblrtm)

    cd ~/MATLABCODE/QUICKTASKS_TELECON/CrIS_gascell_Jan2020/
    saver = ['save ' ipfile '.mat out_array outwave xw xdglab xdlblrtm'];
    eval(saver)

  elseif (gid == 6)
    f1 = 1180;    f2 = 1780;
    f1 = 1180;    f2 = 1780;

    cd /home/sergio/SPECTRA
    [outwave,out_array] = run8(6,f1,f2,ipfile);

    cd  /home/sergio/HITRAN2UMBCLBL/LBLRTM
    [xw,xdglab,xdlblrtm] = driver_glab_lblrtm_self_MANYLAY(6,f1,f2,ipfile,-1);
    semilogy(outwave,out_array,xw,xdlblrtm)

    cd ~/MATLABCODE/QUICKTASKS_TELECON/CrIS_gascell_Jan2020/
    saver = ['save ' ipfile '.mat out_array outwave xw xdglab xdlblrtm'];
    eval(saver)

  elseif (gid == 1)
    error('have not coded this')
  end
end    
