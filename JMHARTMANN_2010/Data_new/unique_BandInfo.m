
%% see /home/sergio/SPECTRA/Global_Data_HITRAN20XY/mapHITRAN_v1v2lv3r.f
%% see co2vibs
dd = load('BandInfo.dat');
dd = dd(:,1);
dd = num2str(dd);

dd_donk = dd(:,2:6);
dd_donk = str2num(dd_donk);

co2vibs = load('/home/sergio/SPECTRA/Global_Data_HITRAN2008/co2vibs');
[Y,I] = intersect(co2vibs(:,2),lala);