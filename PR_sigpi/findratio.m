function matr=findratio(band,prb);
%function matr=findratio(band,prb);
%this finds full/lor for the specified band

%%%%%%%%%%%%%%% ONLY DO THIS FOR 720 BAND!!!!!!!!!!
%%%%%%%%%%%%%%% ALL OTHER BANDS USE 0.5 MIX RATIO!!!

matr=zeros(4,8);


if ((prb == 'P') | (prb == 'p'))
  if (band == 720)
    matr(1,1)=600; matr(2,1)=640; matr(3,1)=750; matr(4,1)=780; 
    end
elseif ((prb == 'R') | (prb == 'r'))
  if (band == 720)
    matr(1,1)=660; matr(2,1)=680; matr(3,1)=800; matr(4,1)=830; 
    end
  end

f=matr(1,1)-100:0.01:matr(4,1)+100;

jj=1;
for T=100:50:400
  T
  jj=jj+1;
  pts=matr(:,1);
  [f,full,fc]=yrun_sigpiOLD(T,band,1,1e-2,1e-4,prb,'F','0','n',f);
  [f,lor,fc]=yrun_sigpiOLD(T,band,1,1e-2,1e-4,prb,'L','0','n',f);  
  rat=full./lor; 
  subplot(211); plot(f,rat);
  subplot(212); plot(f,full,f,lor);
  for ii=1:4
    kk=find(f == pts(ii)); matr(ii,jj)=rat(kk);
    end  
  matr
  pause(0.1)
  end

