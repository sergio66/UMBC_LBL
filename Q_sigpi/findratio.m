function matr=findratio(band);
%function matr=findratio(band);
%this finds full/lor for the specified band

matr=zeros(4,8);

if (band == 618) 
  matr(1,1)=580; matr(2,1)=605; matr(3,1)=625; matr(4,1)=680; 
elseif (band == 648) 
  matr(1,1)=605; matr(2,1)=640; matr(3,1)=665; matr(4,1)=705; 
elseif (band == 662) 
  matr(1,1)=605; matr(2,1)=660; matr(3,1)=680; matr(4,1)=705; 
elseif (band == 667) 
  matr(1,1)=630; matr(2,1)=660; matr(3,1)=680; matr(4,1)=705;  
elseif (band == 720) 
  matr(1,1)=680; matr(2,1)=705; matr(3,1)=730; matr(4,1)=755;  
elseif (band == 791) 
  matr(1,1)=755; matr(2,1)=780; matr(3,1)=800; matr(4,1)=830;  
elseif (band == 2080) 
  matr(1,1)=2030; matr(2,1)=2065; matr(3,1)=2095; matr(4,1)=2130;  
  end 

f=matr(1,1)-20:0.01:matr(4,1)+20;

jj=1;
for T=100:50:400
  T
  jj=jj+1;
  pts=matr(:,1);
  p=1; ps=p/(3030);  q=1e-4;
%  p=1; ps=p/(3030*1.05);  q=1e-4;
%  [f,full,fc]=yrun_sigpiOLD(T,band,p,ps,q,'V','1','n',f);
  [f,full,fc,ss(jj)]=yrun_sigpiOLD(T,band,p,ps,q,'F','0','n',f);
  [f,lor,fc,ss(jj)]=yrun_sigpiOLD(T,band,p,ps,q,'L','0','n',f);  
  rat=full./lor; 
  subplot(211); plot(f,rat);
  subplot(212); plot(f,full,f,lor);
  for ii=1:4
    kk=find(f == pts(ii)); matr(ii,jj)=rat(kk);
    end  
  end

matr
ss
