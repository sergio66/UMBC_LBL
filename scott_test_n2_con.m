% PROGRAM test_n2_con
%
% Tet N2 continuum
% Based on "...N2 near 4.3 um..." by W.Lafferty, A.Solodov, A.Weber,
% B.Olsen, JM.Hartmann.  Applied Optics, v35n30, 20 Oct 1996, pg 5911-5917

clear

TREF = 273.15;
N2fraction = 0.79;


% fig7a
T = TREF + 30.3; % K
P = 990 / 1013.25; % atm
L = 5.7 *1000*100; % cm
%
amount = L*(2.6867E+19*TREF/T)*P/6.022045E+23; % moles/cm^2
[wn,alpha] = equation8(P,T);
fudge = (T/273/P)^2;
od = alpha*L;
trans_7af = exp(-od*fudge);
trans_7a = exp(-od);


% fig7b
T = TREF - 21.4; % K
P = 998 / 1013.25; % atm
L = 5.7 *1000*100; % cm
%
amount = L*(2.6867E+19*TREF/T)*P/6.022045E+23; % moles/cm^2
[wn,alpha] = equation8(P,T);
fudge = (T/273/P)^2;
od = alpha*L;
trans_7bf = exp(-od*fudge);
trans_7b = exp(-od);

figure(1)
clf
plot(wn,trans_7a,'b', wn,trans_7af,'c', ...
     wn,trans_7b,'r', wn,trans_7bf,'m'),axis([2200 2500 0.3 0.9]),grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sergios test profile
T = 287.497;   % Kelvin
P = 0.98687;   % atmospheres
PP = 0.77075;  % atmospheres
amount = 7.6861E-4; % kilomoles/cm^2 (maybe)
%
L = amount/((2.6867E+19*TREF/T)*PP/6.022045E+26)
L = 23500      % cm, approximate
%
[wn,alpha] = equation8(P,T);
fudge = (T/273/P)^2;
od = alpha*L;

figure(2)
clf
plot(wn,od,'b',wn,od*fudge,'c'),axis([2300 2440 0.015 0.04]),grid

%%% end of program %%%

%************************************************************************

n2y   %% this is sergio

clf
plot(wn,od,'b',wn,od*fudge,'c',...
     xx,odx,'r',xx,odx/am,'ro'),axis([2300 2440 0.015 0.04]),grid
