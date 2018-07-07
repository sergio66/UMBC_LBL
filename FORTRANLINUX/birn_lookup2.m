%script file to produce fortran look up table

%this shows that if T is fixed, there is VERY WEAK dependance on w_tot
%for ii=1:10
%  pow=ii-1;                                                       
%  y150(ii,:)=birnbaumWORKS2(000:0.05:3000,1500,10^(-pow),150,3e-3);    
%  y200(ii,:)=birnbaumWORKS2(000:0.05:3000,1500,10^(-pow),200,3e-3);  
%  y250(ii,:)=birnbaumWORKS2(000:0.05:3000,1500,10^(-pow),250,3e-3); 
%  y300(ii,:)=birnbaumWORKS2(000:0.05:3000,1500,10^(-pow),300,3e-3); 
%  y350(ii,:)=birnbaumWORKS2(000:0.05:3000,1500,10^(-pow),350,3e-3); 

%  yx150(ii,:)=birnbaumWORKS2(000:0.05:3000,1500,10^(-pow),150,6e-3); 
%  yx200(ii,:)=birnbaumWORKS2(000:0.05:3000,1500,10^(-pow),200,6e-3); 
%  yx250(ii,:)=birnbaumWORKS2(000:0.05:3000,1500,10^(-pow),250,6e-3); 
%  yx300(ii,:)=birnbaumWORKS2(000:0.05:3000,1500,10^(-pow),300,6e-3); 
%  yx350(ii,:)=birnbaumWORKS2(000:0.05:3000,1500,10^(-pow),350,6e-3); 

%  end

%ll=1:length(y150);
%plot(ll,y150)

%since w_tot is immaterial
%this does a scan of tau2 for 5 temps
for ii=1:5                                                                  
  pow=1;                                                                 
  doc=1e-3*(3+(ii-1)*2);                                  
  y100(ii,:)=birnbaumWORKS2(000:0.25:500,250,10^(-pow),100,doc);
  y150(ii,:)=birnbaumWORKS2(000:0.25:500,250,10^(-pow),150,doc);
  y200(ii,:)=birnbaumWORKS2(000:0.25:500,250,10^(-pow),200,doc);
  y250(ii,:)=birnbaumWORKS2(000:0.25:500,250,10^(-pow),250,doc);
  y300(ii,:)=birnbaumWORKS2(000:0.25:500,250,10^(-pow),300,doc);
  y350(ii,:)=birnbaumWORKS2(000:0.25:500,250,10^(-pow),350,doc);
  y400(ii,:)=birnbaumWORKS2(000:0.25:500,250,10^(-pow),400,doc);
  end

clear ii pow doc
%save birnbaum.txt y100 y150 y200 y250 y300 y350 y400 -ascii -double
%save birnbaum.mat
fid=fopen('birn_lookup2.dat','w');

%%%%%%%%%%%%%%%%%%%%
ughugh=y100;
fprintf(fid,'c this is for T = 100 \n');
for ii=1:99
  index=(1:20)+(ii-1)*20;
  index1=index(1):index(length(index))-1;
  index2=index(length(index));
  kk=index(1);
  ll=index(length(index));
  ugh=ughugh(:,index1);
  ugh1=ughugh(:,index2);
  fprintf(fid,'      data ((chi100(i,j),j=1,5),i=%4i,%4i)/ \n',kk,ll);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh1);
  end
fprintf(fid,'      data ((chi100(i,j),j=1,5),i=1981,1990)/ \n',kk);
index=1981:1989;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=1990:1990;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

fprintf(fid,'      data ((chi100(i,j),j=1,5),i=1991,2001)/ \n',kk);
index=1991:2000;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=2001:2001;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

%%%%%%%%%%%%%%%%%%%%
ughugh=y150;
fprintf(fid,'c this is for T = 150 \n');
for ii=1:99
  index=(1:20)+(ii-1)*20;
  index1=index(1):index(length(index))-1;
  index2=index(length(index));
  kk=index(1);
  ll=index(length(index));
  ugh=ughugh(:,index1);
  ugh1=ughugh(:,index2);
  fprintf(fid,'      data ((chi150(i,j),j=1,5),i=%4i,%4i)/ \n',kk,ll);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh1);
  end
fprintf(fid,'      data ((chi150(i,j),j=1,5),i=1981,1990)/ \n',kk);
index=1981:1989;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=1990:1990;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

fprintf(fid,'      data ((chi150(i,j),j=1,5),i=1991,2001)/ \n',kk);
index=1991:2000;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=2001:2001;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

%%%%%%%%%%%%%%%%%%%%
ughugh=y200;
fprintf(fid,'c this is for T = 200 \n');
for ii=1:99
  index=(1:20)+(ii-1)*20;
  index1=index(1):index(length(index))-1;
  index2=index(length(index));
  kk=index(1);
  ll=index(length(index));
  ugh=ughugh(:,index1);
  ugh1=ughugh(:,index2);
  fprintf(fid,'      data ((chi200(i,j),j=1,5),i=%4i,%4i)/ \n',kk,ll);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh1);
  end
fprintf(fid,'      data ((chi200(i,j),j=1,5),i=1981,1990)/ \n',kk);
index=1981:1989;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=1990:1990;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

fprintf(fid,'      data ((chi200(i,j),j=1,5),i=1991,2001)/ \n',kk);
index=1991:2000;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=2001:2001;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

%%%%%%%%%%%%%%%%%%%%
ughugh=y250;
fprintf(fid,'c this is for T = 250 \n');
for ii=1:99
  index=(1:20)+(ii-1)*20;
  index1=index(1):index(length(index))-1;
  index2=index(length(index));
  kk=index(1);
  ll=index(length(index));
  ugh=ughugh(:,index1);
  ugh1=ughugh(:,index2);
  fprintf(fid,'      data ((chi250(i,j),j=1,5),i=%4i,%4i)/ \n',kk,ll);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh1);
  end
fprintf(fid,'      data ((chi250(i,j),j=1,5),i=1981,1990)/ \n',kk);
index=1981:1989;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=1990:1990;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

fprintf(fid,'      data ((chi250(i,j),j=1,5),i=1991,2001)/ \n',kk);
index=1991:2000;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=2001:2001;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

%%%%%%%%%%%%%%%%%%%%
ughugh=y300;
fprintf(fid,'c this is for T = 300 \n');
for ii=1:99
  index=(1:20)+(ii-1)*20;
  index1=index(1):index(length(index))-1;
  index2=index(length(index));
  kk=index(1);
  ll=index(length(index));
  ugh=ughugh(:,index1);
  ugh1=ughugh(:,index2);
  fprintf(fid,'      data ((chi300(i,j),j=1,5),i=%4i,%4i)/ \n',kk,ll);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh1);
  end
fprintf(fid,'      data ((chi300(i,j),j=1,5),i=1981,1990)/ \n',kk);
index=1981:1989;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=1990:1990;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

fprintf(fid,'      data ((chi300(i,j),j=1,5),i=1991,2001)/ \n',kk);
index=1991:2000;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=2001:2001;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

%%%%%%%%%%%%%%%%%%%%
ughugh=y350;
fprintf(fid,'c this is for T = 350 \n');
for ii=1:99
  index=(1:20)+(ii-1)*20;
  index1=index(1):index(length(index))-1;
  index2=index(length(index));
  kk=index(1);
  ll=index(length(index));
  ugh=ughugh(:,index1);
  ugh1=ughugh(:,index2);
  fprintf(fid,'      data ((chi350(i,j),j=1,5),i=%4i,%4i)/ \n',kk,ll);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh1);
  end
fprintf(fid,'      data ((chi350(i,j),j=1,5),i=1981,1990)/ \n',kk);
index=1981:1989;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=1990:1990;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

fprintf(fid,'      data ((chi350(i,j),j=1,5),i=1991,2001)/ \n',kk);
index=1991:2000;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=2001:2001;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

%%%%%%%%%%%%%%%%%%%%
ughugh=y400;
fprintf(fid,'c this is for T = 400 \n');
for ii=1:99
  index=(1:20)+(ii-1)*20;
  index1=index(1):index(length(index))-1;
  index2=index(length(index));
  kk=index(1);
  ll=index(length(index));
  ugh=ughugh(:,index1);
  ugh1=ughugh(:,index2);
  fprintf(fid,'      data ((chi400(i,j),j=1,5),i=%4i,%4i)/ \n',kk,ll);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
  fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh1);
  end
fprintf(fid,'      data ((chi400(i,j),j=1,5),i=1981,1990)/ \n',kk);
index=1981:1989;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=1990:1990;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

fprintf(fid,'      data ((chi400(i,j),j=1,5),i=1991,2001)/ \n',kk);
index=1991:2000;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
index=2001:2001;
ugh=ughugh(:,index);
fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh);

%%%%%%%%%%%%%%%%%%%%

%ughugh=y400;
%ugh1=ughugh(:,2001);
%ugh=ughugh(:,1:2000);
%fprintf(fid,'c this is for T = 400 \n');
%fprintf(fid,'      data (((chi(i,j,k),k=7,7),j=1,5),i=1,2001)/ \n');
%fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e, \n',ugh);
%fprintf(fid,'     + %6.4e, %6.4e, %6.4e, %6.4e, %6.4e/ \n',ugh1);

fclose(fid);

