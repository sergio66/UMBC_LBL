function [brd]=broad_weakco2(press,press_self,press_ref,air,self,m,T,iGas)

%function [brd]=broad_weakco2(press,press_self,press_ref,air,self,m,T,iGas)
%for the run7co2_weakback.m code

%compute the broadening by combing air and self broadening
%remember units are in cm-1 per atm at 296 K, so we need the pressures
% press      = current AIRS pressure in atm
% press_self = current self pressure in atm
% press_ref  = current reference pressure in atm
% air        = air broadening cm-1/atm at 296 k
% self       = air broadening cm-1/atm at 296 k
% m          = power relationship to scale brd wrt temperature
% T          = temperature
%iGas        = GAS ID

%this eqn is from pg 31 of Genln2 manual

%%%brd=air*(press-press_self)/press_ref + self*press_self/press_ref;
%assume press_ref = 1.0 atm

if (iGas ~= 2)                        %all gases except CO2
  %this is the vectorised code
  dummysmall = (self < eps);
  dummybig   = (self >= eps);
  slfb=self.*dummybig;

  if (sum(dummysmall) > 0)
    if (iGas == 1)                    %for water, self width ~ 5 x air width
      slfs=(5*air).*dummysmall;
    else
      slfs=air.*dummysmall;
      end
  else
    slfs=zeros(size(self));
    end
  slf=slfs+slfb;

  brd=air*(press-press_self) + slf*press_self;
  brd=(296.0/T).^(m).*brd;

else                                %%%%%%%%%%for CO2
  %this is the vectorised code
  dummysmall = (self < eps);
  dummybig   = (self >= eps);
  slfb=self.*dummybig;

  if (sum(dummysmall) > 0)
    slfs=air.*dummysmall;
  else
    slfs=zeros(size(self));
    end
  slf=slfs+slfb;

  brdf=(press-press_self)*air.*(296.0/T).^(m);
  %brds=(press_self)*slf.*(296.0/T).^(0.685); %%%% <---- was 0.685 before
  brds=(press_self)*slf.*(296.0/T).^(m);      %%%% <---- was 0.685 before
  brd=brdf+brds;
  end

