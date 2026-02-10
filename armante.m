topts.HITRAN = '/asl/data/geisa/g15.by.gas/';;
topts.HITRAN = [geisapath '/g15.by.gas/'];

[w,dx] = run8(3,975,1225,'IPFILES/co2one',topts);
[w,d] = run8(3,975,1225,'IPFILES/co2one');
save armante.mat w d dx topts

subplot(211); plot(w,dx,w,d)
  title('p=1.1139 atm  pO3=4.0418E-03 atm  T=269.429K');; ylabel('OD')
  grid
subplot(212); plot(w,dx./d)
  ylabel('OD geisa/OD hitran')
  grid

data = [w; d; dx];
whos data
fid = fopen('armante.txt','w');
fprintf(fid,'%12.5f %12.7e %12.7e \n',data);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

topts.HITRAN = '/asl/data/geisa/g15.by.gas/';;
topts.HITRAN = [geisapath '/g15.by.gas/'];

[w,dx] = run8(3,1200,1275,'IPFILES/co2one',topts);
[w,d] = run8(3,1200,1275,'IPFILES/co2one');
save armante2.mat w d dx topts

subplot(211); plot(w,dx,w,d)
  title('p=1.1139 atm  pO3=4.0418E-03 atm  T=269.429K');; ylabel('OD')
  grid
subplot(212); plot(w,dx./d)
  ylabel('OD geisa/OD hitran')
  grid

data = [w; d; dx];
whos data
fid = fopen('armante2.txt','w');
fprintf(fid,'%12.5f %12.7e %12.7e \n',data);
fclose(fid);

