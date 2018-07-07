%%%%%%%%%%%%%%%%%%%%%%%%%  efitter.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('fitting for missing odd energy levels')
	bds=[0;0.39;1.5e-7];
	ebd=leastsq('efit',bds,[],[],jr,elowerr)
	jallr=[0:max(jr)+1]';evib=ebd(1);
	elower_r= (ebd(2)*jallr.*(jallr+1) - ebd(3)*(jallr.*(jallr+1)).^2);
        index_e=find(rem(jallr,2)==0);
	elower_r(index_e)=elowerr-evib;

	bds=[0;0.39;1.5e-7];
	ebd=leastsq('efit',bds,[],[],jp,elowerp)
	jallp=[1:max(jp)+1]';evib=ebd(1);
	elower_p= (ebd(2)*jallp.*(jallp+1) - ebd(3)*(jallp.*(jallp+1)).^2);
        index_e=find(rem(jallp,2)==0);
	elower_p(index_e)=elowerp-evib;

clear start2 bds elowerrf
