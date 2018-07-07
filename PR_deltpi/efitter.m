function [elower,jall] = efitter(jq,band,elowerq,elower,prb);

global v12p1

%disp('fitting for missing odd energy levels')

bds = [668;  % Evib
       3.916693061746631e-01;  % B
       1.390937512910130e-07]; % D

if band == 740
  if v12p1 > 0
    options(1) = 0;
    ebd = leastsq('efit',bds,options,[],jq,elowerq);
  else
    clear options; options = optimset('MaxFunEvals',1000);
    %ebd = lsqnonlin(@efit,bds,[],[],options,jq,elowerq);
    ebd = bds;   
    params = [jq elowerq];   
    ebd = fsolve(@(ebd) efit1(ebd,params),bds,options);   
    end
  end

if (band == 668 | band == 2093)
  ebd = bds;
  end

if (rem(max(jq),2) == 0)
  jall = [1:max(jq)+1]'; 	
else
  jall = [1:max(jq)]'; 	
  end

evib = ebd(1);	
elower =  (ebd(2)*jall.*(jall+1) - ebd(3)*(jall.*(jall+1)).^2);
index_q = jq;
elower(index_q) = elowerq-evib;

eloweref = (ebd(2)*jq.*(jq+1) - ebd(3)*(jq.*(jq+1)).^2);

clear start2 bds eloweref
