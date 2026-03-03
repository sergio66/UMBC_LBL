%%% this is in orig code in run8 since 1996
c2 = 1.4387863;               %% second radiation constant
AVOG = 6.022045E+26;          %% kilomolecules/cm2
PREF = 101325.0;              %% N/m2
TREF = 296.0;                 %% K
T273 = 273.15;                %% K
MGC = 8.314674269981136;      %% J/mol/K
%%% this is in orig code in run8 since 1996

%% this is in MT_CKD_H2O-4.3/src/phys_consts.f90
AVOG = 6.02214199E+23 * 1000; %% kilomolecules/cm2
C2 = 1.4387752;               %% second radiation constant
PREF = 101300.0;              %% N/m2
TREF = 296.0;                 %% K
T273 = 273.15;                %% K
MGC = 8.314674269981136;      %% J/mol/K
%% this is in MT_CKD_H2O-4.3/src/phys_consts.f90

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tc2k = T273;

torr2mb  = (PREF/100) / 760; 
torr2atm = 1 / 760; 
mb2atm   = 1 / (PREF/100); 

