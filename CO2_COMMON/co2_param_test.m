%c -rw-r--r--  1 sergio users 8550 2002-08-09 15:43 co2_param_AIRS2002.m
%c -rw-r--r--  1 sergio users 7359 2002-08-08 13:26 co2_param_JJOHNS2002.m
%c -rw-r--r--  1 sergio users 7099 2002-07-30 15:49 co2_param_JULY20_2002.m
% -rw-r--r--  1 sergio users 5499 2002-07-04 12:47 co2_param_JJOHNS2002_JULY3.m
%c -rw-r--r--  1 sergio users 5103 2002-07-03 10:27 co2_param_JJOHNS2002_ORIG.m
%c -rw-r--r--  1 sergio users 5668 2000-08-03 16:41 co2_param_RAL_NEW.m
%c -rw-r--r--  1 sergio users 7917 2000-08-02 15:35 co2_param_RAL_NEWALL.m
% -rw-r--r--  1 sergio users 3658 1999-08-09 18:33 co2_param_JOHNJOHN.m
% -rw-r--r--  1 sergio users 3086 1999-07-06 17:07 co2_param_RAL.m

pt = 1.0712E+00/1013.25;
ps = 3.8885E-04/1013.25;

pt = 1.0712E+00;
ps = 3.8885E-04;

band = 668;
band = 720;
band = 740;

ii=1;
[duration_self(ii),duration_for(ii),...
 beta_self(ii),beta_for(ii),b2s(ii),b2f(ii)] = ...
  co2_param_AIRS2002(band,pt,ps);   %%% <---- co2_param.m was linked to this

ii=2;
[duration_self(ii),duration_for(ii),...
  beta_self(ii),beta_for(ii),b2s(ii),b2f(ii)] = ...
  co2_param_JULY20_2002(band,pt,ps);

ii=3;
[duration_self(ii),duration_for(ii),...
  beta_self(ii),beta_for(ii),b2s(ii),b2f(ii)] = ...
  co2_param_JJOHNS2002(band,pt,ps,1);

ii=4;
[duration_self(ii),duration_for(ii),...
  beta_self(ii),beta_for(ii),b2s(ii),b2f(ii)] = ...
  co2_param_JJOHNS2002_JULY3(band,pt,ps);

ii=5;
[duration_self(ii),duration_for(ii),...
  beta_self(ii),beta_for(ii),b2s(ii),b2f(ii)] = ...
  co2_param_JJOHNS2002_ORIG(band,pt,ps);

ii=6;
[duration_self(ii),duration_for(ii),...
  beta_self(ii),beta_for(ii),b2s(ii),b2f(ii)] = ...
  co2_param_RAL_NEW(band,pt,ps);

%********* these are the odd one out *****************
% for 668,740,2093
ii=7;
[duration_self(ii),duration_for(ii),...
  beta_self(ii),beta_for(ii),b2s(ii),b2f(ii)] = ...
  co2_param_RAL_NEWALL(band,pt,ps);

ii=8;
[duration_self(ii),duration_for(ii),...
  beta_self(ii),beta_for(ii),b2s(ii),b2f(ii)] = ...
  co2_param_JOHNJOHN(band,pt,ps);

ii=9;
[duration_self(ii),duration_for(ii),...
  beta_self(ii),beta_for(ii),b2s(ii),b2f(ii)] = ...
  co2_param_RAL(band);
%********* this is the oddone out *****************

plot(1:9,[beta_self; beta_for],'+-'); title('betaself betafor')
plot(1:9,[b2s      ; b2f],'+-'); title('b2s b2f')
plot(1:9,[duration_self; duration_for],'+-'); title('durationself durationfor')
