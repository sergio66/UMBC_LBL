function matr=findratio(band,prb);
%function matr=findratio(band);
%this finds full/lor for the specified band

matr=zeros(4,8);

if ((prb == 'P') | (prb == 'p'))
  if (band == 2350)
    matr(1,1)=2105; matr(2,1)=2155; matr(3,1)=2480; matr(4,1)=2530; 
    end
elseif ((prb == 'R') | (prb == 'r'))
  if (band == 2350)
    matr(1,1)=2205; matr(2,1)=2255; matr(3,1)=2480; matr(4,1)=2530; 
    end
  end

f=matr(1,1)-20:0.01:matr(4,1)+20;

jj=1;
for T=100:50:400
  T
  jj=jj+1;
  pts=matr(:,1);
  p=1; ps=1e-2; q=1e-4;
  p=1; ps=p/(1000); q=1e-4;
%  [f,full,fc]=yrun_sigsigOLD(T,band,p,ps,q,,prb,'V','1','n',f);
  [f,full,fc,sss]=yrun_sigsigOLD(T,band,p,ps,q,prb,'F','0','n',f);
  [f,lor,fc,sss]=yrun_sigsigOLD(T,band,p,ps,q,prb,'L','0','n',f);  
  rat=full./lor;
  ss(jj)=sss; 
  subplot(211); plot(f,rat);
  subplot(212); plot(f,full,f,lor);
  for ii=1:4
    kk=find(f == pts(ii)); matr(ii,jj)=rat(kk);
    end  
  end


format short e
matr
ss
