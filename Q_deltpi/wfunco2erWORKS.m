function [W_co2for,W_co2self] = wfunco2er(jq,elower,elowerq,w_selfq,w_forq,...
       band,jall,temperature,stuff,beta_pi,beta_delt); 

%%%%%%%%%%%%%%%%%%%%%% wfunco2er.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global v12p1

xx=[ -6.202132765289518e-02   3.455577145107923e-01   1.009589099173803e+00];

if band ~= 2093
  if v12p1 == +1
    options(1) = 0; 
    options(5) = 0; options(7) = 1; options(14) = 60;
    a1a2a3_for = leastsq('wfunco2',xx,options,'wgradco2',elower,w_forq,jall,...
        stuff.beta_pi_air,stuff.beta_delt_air,temperature,stuff);
  elseif v12p1 == -1
    clear options; options = optimset('MaxFunEvals',1000);
    a1a2a3_for = lsqnonlin(@wfunco2,xx,[],[],options,elower,w_forq,jall,...
	stuff.beta_pi_air,stuff.beta_delt_air,temperature,stuff);
  elseif v12p1 == 0
    clear options; options = optimset('MaxFunEvals',1000);
    xdata = elower;
    ydata = wq_tape;
    a1a2a3_for = lsqcurvefit(@(a1a2a3_for,xdata) wfunco2x(xx,xdata,jall,...
	stuff.beta_pi_air,stuff.beta_delt_air,temperature,stuff)),xx,xdata,
        ydata,options);
    end

elseif band == 2093
  a1a2a3_for = [-4.672362662290493e-02
	       3.332699473532611e-01
	       1.063319601751937e+00];
  end

%a1a2a3_self= [-7.717999298395123e-02     
%               3.668221452704686e-01
%               1.156757636655336e+00]
%a1a2a3_for = [-4.654825090171447e-02  
%               3.316590519181703e-01
%               1.063498003168863e+00]
W_co2for = wfun1co2(a1a2a3_for,elowerq,w_forq,jq,temperature,stuff);

if band ~= 2093
  if v12p1 > 0
    a1a2a3_self = ...
        leastsq('wfunco2',xx,options,'wgradco2',elower,w_selfq,jall,...
               stuff.beta_pi_self,stuff.beta_delt_self,temperature,stuff);
  else
    a1a2a3_self = lsqnonlin(@wfunco2,xx,[],[],options,elower,w_selfq,jall,...
               stuff.beta_pi_self,stuff.beta_delt_self,temperature,stuff);
    end
elseif band == 2093
  a1a2a3_self = [-7.223394382342538e-02
	        3.601014807073620e-01
	        1.162424400557535e+00];
  end

%a1a2a3_self= [-7.717999298395123e-02
%               3.668221452704686e-01
%               1.156757636655336e+00]
%a1a2a3_for=[ -4.654825090171447e-02
%              3.316590519181703e-01
%              1.063498003168863e+00]
W_co2self = wfun1co2(a1a2a3_self,elowerq,w_selfq,jq,temperature,stuff);

