function matr=findratio(band);
%function matr=findratio(band);
%this finds full/lor for the specified band

matr=zeros(4,8);

if ((prb == 'P') | (prb == 'p'))
  if (band == 2310)
    matr(1,1)=2220; matr(2,1)=2225; matr(3,1)=2330; matr(4,1)=2335; 
    end
elseif ((prb == 'R') | (prb == 'r'))
  if (band == 2310)
    matr(1,1)=2320; matr(2,1)=2325; matr(3,1)=2375; matr(4,1)=2380; 
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

