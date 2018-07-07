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
  end

clear ii pow doc
whos

[mm,nn]=size(y100);
mm=max(mm,nn);

filename='bbiirrnn.dat';
fid=fopen([filename],'w','native');

filemark4=8;
filemark_vect=mm*8;

temp=100;         %write the temperature
fwrite(fid,filemark4,'integer*4');
fwrite(fid,temp,'real*8');
fwrite(fid,filemark4,'integer*4');
%loop over the rows
for ii=1:10
  fwrite(fid,filemark_vect,'integer*4');
  fwrite(fid,y100(ii,:),'real*8');
  fwrite(fid,filemark_vect,'integer*4');
  end

temp=150;         %write the temperature
fwrite(fid,filemark4,'integer*4');
fwrite(fid,temp,'real*8');
fwrite(fid,filemark4,'integer*4');
%loop over the rows
for ii=1:10
  fwrite(fid,filemark_vect,'integer*4');
  fwrite(fid,y150(ii,:),'real*8');
  fwrite(fid,filemark_vect,'integer*4');
  end

temp=200;         %write the temperature
fwrite(fid,filemark4,'integer*4');
fwrite(fid,temp,'real*8');
fwrite(fid,filemark4,'integer*4');
%loop over the rows
for ii=1:10
  fwrite(fid,filemark_vect,'integer*4');
  fwrite(fid,y200(ii,:),'real*8');
  fwrite(fid,filemark_vect,'integer*4');
  end

temp=250;         %write the temperature
fwrite(fid,filemark4,'integer*4');
fwrite(fid,temp,'real*8');
fwrite(fid,filemark4,'integer*4');
%loop over the rows
for ii=1:10
  fwrite(fid,filemark_vect,'integer*4');
  fwrite(fid,y250(ii,:),'real*8');
  fwrite(fid,filemark_vect,'integer*4');
  end

temp=300;         %write the temperature
fwrite(fid,filemark4,'integer*4');
fwrite(fid,temp,'real*8');
fwrite(fid,filemark4,'integer*4');
%loop over the rows
for ii=1:10
  fwrite(fid,filemark_vect,'integer*4');
  fwrite(fid,y300(ii,:),'real*8');
  fwrite(fid,filemark_vect,'integer*4');
  end

temp=350;         %write the temperature
fwrite(fid,filemark4,'integer*4');
fwrite(fid,temp,'real*8');
fwrite(fid,filemark4,'integer*4');
%loop over the rows
for ii=1:10
  fwrite(fid,filemark_vect,'integer*4');
  fwrite(fid,y350(ii,:),'real*8');
  fwrite(fid,filemark_vect,'integer*4');
  end

temp=400;         %write the temperature
fwrite(fid,filemark4,'integer*4');
fwrite(fid,temp,'real*8');
fwrite(fid,filemark4,'integer*4');
%loop over the rows
for ii=1:10
  fwrite(fid,filemark_vect,'integer*4');
  fwrite(fid,y400(ii,:),'real*8');
  fwrite(fid,filemark_vect,'integer*4');
  end

fclose(fid);
