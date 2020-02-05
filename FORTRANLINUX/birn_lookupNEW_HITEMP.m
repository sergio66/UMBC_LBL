%script file to produce fortran look up table for LINUX
%same as birn_lookup.m except 1) has more tau
%                             2)writes to a Fortran DATA file

%this shows that if T is fixed, there is VERY WEAK dependance on w_tot
%for ii=1:10
%  pow=ii-1;                                                       
%  y150(ii,:)=birnbaumWORKS(000:0.05:3000,1500,10^(-pow),150,3e-3);    
%  y200(ii,:)=birnbaumWORKS(000:0.05:3000,1500,10^(-pow),200,3e-3);  
%  y250(ii,:)=birnbaumWORKS(000:0.05:3000,1500,10^(-pow),250,3e-3); 
%  y300(ii,:)=birnbaumWORKS(000:0.05:3000,1500,10^(-pow),300,3e-3); 
%  y350(ii,:)=birnbaumWORKS(000:0.05:3000,1500,10^(-pow),350,3e-3); 

%  yx150(ii,:)=birnbaumWORKS(000:0.05:3000,1500,10^(-pow),150,6e-3); 
%  yx200(ii,:)=birnbaumWORKS(000:0.05:3000,1500,10^(-pow),200,6e-3); 
%  yx250(ii,:)=birnbaumWORKS(000:0.05:3000,1500,10^(-pow),250,6e-3); 
%  yx300(ii,:)=birnbaumWORKS(000:0.05:3000,1500,10^(-pow),300,6e-3); 
%  yx350(ii,:)=birnbaumWORKS(000:0.05:3000,1500,10^(-pow),350,6e-3); 

%  end

%ll=1:length(y150);
%plot(ll,y150)

%since w_tot is immaterial
%this does a scan of tau2 for 10 docs
for ii=1:10                                                                  
  pow=1;                                                                 
  doc=1e-3*(1+(ii-1)*3);                                  
  y100(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),100,doc);
  y150(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),150,doc);
  y200(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),200,doc);
  y250(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),250,doc);
  y300(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),300,doc);
  y350(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),350,doc);
  y400(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),400,doc);

  y450(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),450,doc);
  y500(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),500,doc);
  y550(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),550,doc);
  y600(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),600,doc);
  y650(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),650,doc);
  y700(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),700,doc);
  y750(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),750,doc);
  y800(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),800,doc);
  y850(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),850,doc);
  y900(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),900,doc);
  y950(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),950,doc);

  y1000(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),1000,doc);
  y1050(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),1050,doc);
  y1100(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),1100,doc);
  y1150(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),1150,doc);
  y1200(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),1200,doc);
  y1250(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),1250,doc);
  y1300(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),1300,doc);
  y1350(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),1350,doc);
  y1400(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),1400,doc);
  y1450(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),1450,doc);
  y1500(ii,:)=birnbaumWORKS(000:0.25:500,250,10^(-pow),1500,doc);
end

clear ii pow doc
whos

[mm,nn]=size(y100);
mm=max(mm,nn);

filename='bbiirrnn_HITEMP.dat';
fid=fopen([filename],'w','native');

filemark4=8;
filemark_vect=mm*8;

for tt = 100 : 50 : 1500
  temp=tt;         %write the temperature
  fprintf(1,'writing out T = %4i \n',temp);
  fwrite(fid,filemark4,'integer*4');
  fwrite(fid,temp,'real*8');
  fwrite(fid,filemark4,'integer*4');

  %loop over the rows
  str = ['junk = y' num2str(temp) ';'];
  eval(str);
  for ii=1:10
    fwrite(fid,filemark_vect,'integer*4');
    fwrite(fid,junk(ii,:),'real*8');
    fwrite(fid,filemark_vect,'integer*4');
  end
end

fclose(fid);
