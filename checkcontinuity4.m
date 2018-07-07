 fa=1955;
fb=fa+25.0;
delta=0.0005;
spacing=delta*5;

len=ceil(((fb-spacing)-fa)/spacing)+1

freq=[];
k=[];
kn2=[];
kco=[];
for ii=1:22
  index=(1:len)+(ii-1)*len;
  [fr,k6]=run6co2(2,fa,fb,delta,0.1,0.5,1,1,2,150,5,1e-28,1e-28,'f','1','b','../SPECTRA/IPFILES/RAL_CO2/ral4some');
  [fr,kn]=run6(22,fa,fb,delta,0.1,0.5,1,1,2,25,5,1e-28,1e-28,'v',1,'../SPECTRA/IPFILES/RAL_CO2/ral4n2');
  [fr,kc]=run6(5, fa,fb,delta,0.1,0.5,1,1,2,25,5,1e-28,1e-28,'v',1,'../SPECTRA/IPFILES/RAL_CO2/ral4some');

  freq=[freq fr];
  k=[k k6];
  kn2=[kn2 kn];
  kco=[kco kc];

  fa=fa+25.0;
  fb=fb+25.0;

  end

knit=zeros(4,length(kn2)); knit(2:4,:)=kn2;

ktotal=k+knit+kco*5e-7;
ktotal=exp(-ktotal);

cd /salsify/scratch3/Sergio/RAL_DATA/MATLAB
[rch, wch] = fconvkc(ktotal', 'ralB22', 'nb', 6);
