[fr4_1,cous4_1]=run6co2(2,2255,2355,0.0005,0.1,0.5,1,1,2,150,5,...
   1e-28,1e-28,'f','1','b','../SPECTRA/IPFILES/RAL_CO2/davethesis_co2_4um'); 

[fr4_2,cous4_2]=run6co2(2,2355,2455,0.0005,0.1,0.5,1,1,2,150,5,...
 1e-28,1e-28,'f','1','b','../SPECTRA/IPFILES/RAL_CO2/davethesis_co2_4um'); 

[fr4_3,cous4_3]=run6co2(2,2455,2555,0.0005,0.1,0.5,1,1,2,150,5,...
   1e-28,1e-28,'f','1','b','../SPECTRA/IPFILES/RAL_CO2/davethesis_co2_4um'); 

fr4      =[fr4_1       fr4_2       fr4_3];
cous4    =[cous4_1     cous4_2     cous4_3];
save /salsify/scratch3/Sergio/TestRun6CO2/dave4cous.mat fr4 cous4

