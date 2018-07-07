%% this checks differences between hartmann and hitran for 662 band

clear zz662 hh 
zz662 = load('~/SPECTRA/CO2_MATFILES/hit662');
hh    = load('~/SPECTRA/JMHARTMANN/LM_PQR_CO2_2.0/Data/S3002001.dat');

whos zz662 hh
plot(1:length(zz662.v_lower),zz662.v_lower,...
     1:length(zz662.v_lower),zz662.v_upper)

jvalue = str2num(zz662.j_lower(:,6:8));
plot(zz662.freq,jvalue,'+',hh(:,1),hh(:,8),'r+'); title('b=hitran, r=hartman')
