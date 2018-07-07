
fa=605;
fb=fa+25.0;
delta=0.0005;
spacing=delta*5;

len=ceil(((fb-spacing)-fa)/spacing)+1

freq=[];
k=[];
for ii=1:12
  index=(1:len)+(ii-1)*len;
  [fr,k6]=run6co2(2,fa,fb,delta,0.1,0.5,1,1,2,150,5,1e-28,1e-28,'f','1','b','../SPECTRA/IPFILES/RAL_CO2/ral15');
  freq=[freq fr];
  k=[k k6];
  fa=fa+25.0;
  fb=fb+25.0;
  end

ktotal=k;
ktotal=exp(-ktotal); 
 
cd /salsify/scratch3/Sergio/RAL_DATA/MATLAB 
[rch, wch] = fconvkc(ktotal', 'ralB33', 'nb', 6); 
