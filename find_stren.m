function  [strength]=find_stren(qfcn,v0,T,E_li,s0,amt)
% function  [strength]=find_stren(qfcn,v0,T,E_li,s0,amt)
% renormalises the strength based on eqn in pg 31 of Genln2 manual
% qfcn = Q(T)/Q(296)
% v0   = central wavenumber
% T    = temperature
% s0   = strength
% E_li = lower state energy
% amt  = gas amt (kilomolecules/ cm2)

s00 = s0*6.022045e26;                           % or could do amt=amt*6.022e26
c2  = 1.4387863;                                % K/ cm-1  from Genln2 manual
sb  = exp(-c2*E_li/T)./exp(-c2*E_li/296.0);     % boltzman factor (distribution at T)
se  = (1-exp(-c2*v0/T))./(1-exp(-c2*v0/296.0)); % adjust for detailed balance

strength = amt*(qfcn'.*s00.*sb.*se);

%eval(['!/bin/rm boo.mat']);
%save boo.mat qfcn strength