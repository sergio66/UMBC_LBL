%%%%%%%%%%%%%%%%%%%%%% wfunco2er.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx=[ -6.202132765289518e-02   3.455577145107923e-01   1.009589099173803e+00];
options(1)=1; options(5)=0; options(7)=1; options(14)=60;

for i=1:length(pick_no)
  number=i;

  a1a2a3_for=leastsq('wfunco2',xx,options,'wgradco2',elower_r,...
	w_forr(:,i),jallr,'R',number)
  W_co2for=wfun1co2(a1a2a3_for,elowerr,w_forr(:,i),jr,number);
  a1a2a3_self=leastsq('wfunco2',xx,options,'wgradco2',elower_r,...
	w_selfr(:,i),jallr,'R',number)
  W_co2self=wfun1co2(a1a2a3_self,elowerr,w_selfr(:,i),jr,number);
  eval(['W_plus_r' num2str(i) '=(pressure_self(i)*W_co2self+pressure_for(i)*W_co2for)/pressure_ref;']);


  a1a2a3_for=leastsq('wfunco2',xx,options,'wgradco2',elower_p,...
	w_forp(:,i),jallp,'P',number)
  W_co2for=wfun1co2(a1a2a3_for,elowerp,w_forp(:,i),jp,number);
  a1a2a3_self=leastsq('wfunco2',xx,options,'wgradco2',elower_p,...
	w_selfp(:,i),jallp,'P',number)
  W_co2self=wfun1co2(a1a2a3_self,elowerp,w_selfp(:,i),jp,number);
  eval(['W_plus_p' num2str(i) '=(pressure_self(i)*W_co2self+pressure_for(i)*W_co2for)/pressure_ref;'])

  end

clear options














