function matr=findratio(band);
%function matr=findratio(band);
%this finds full/lor for the specified band

matr=zeros(4,8);

if (band == 668)
  matr(1,1)=650; matr(2,1)=660; matr(3,1)=680; matr(4,1)=700;
elseif (band == 740)
  matr(1,1)=700; matr(2,1)=720; matr(3,1)=750; matr(4,1)=770;
elseif (band == 2093)
  matr(1,1)=2075; matr(2,1)=2080; matr(3,1)=2105; matr(4,1)=2120;
  end

f=matr(1,1)-20:0.01:matr(4,1)+20;

jj=1;
for T=100:50:400
  T
  jj=jj+1;
  pts=matr(:,1);
  [f,full,fc]=yrun_deltpiOLD(T,band,1,1e-2,1e-4,'F','0','n',f);
  [f,lor,fc]=yrun_deltpiOLD(T,band,1,1e-2,1e-4,'L','0','n',f);  
  rat=full./lor; 
  subplot(211); plot(f,rat);
  subplot(212); plot(f,full,f,lor);
  for ii=1:4
    kk=find(f == pts(ii)); matr(ii,jj)=rat(kk);
    end  
  matr
  pause(0.1)
  end

