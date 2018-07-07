function matr=findratio(band,prb);
%function matr=findratio(band,prb);
%this finds full/lor for the specified band

matr=zeros(4,8);

if ((prb == 'P') | (prb == 'p'))
  if (band == 2380)
    matr(1,1)=2220; matr(2,1)=2230; matr(3,1)=2380; matr(4,1)=2390; 
    end
elseif ((prb == 'R') | (prb == 'r'))
  if (band == 2380)
    matr(1,1)=2300; matr(2,1)=2305; matr(3,1)=2430; matr(4,1)=2435; 
    end
  end

f=matr(1,1)-20:0.01:matr(4,1)+20;

jj=1;
for T=100:50:400
  T
  jj=jj+1;
  pts=matr(:,1);
%  [f,full,fc]=yrun_sigpiOLD(T,band,1,1e-2,1e-4,prb,'V','1','n',f);
  [f,full,fc]=yrun_sigpiOLD(T,band,1,1e-2,1e-4,prb,'F','0','n',f);
  [f,lor,fc]=yrun_sigpiOLD(T,band,1,1e-2,1e-4,prb,'L','0','n',f);  
  rat=full./lor; 
  subplot(211); plot(f,rat);
  subplot(212); plot(f,full,f,lor);
  for ii=1:4
    kk=find(f == pts(ii)); matr(ii,jj)=rat(kk);
    end
  end

format short e 
matr 
