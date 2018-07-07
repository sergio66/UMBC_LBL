function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
         co_param(band,ptotal,pself)
%function [duration_self,duration_for,beta_self,beta_for,b2s,b2f] = ...
%         co_param(band,ptotal,pself)


%b2s,b2f=0 for all bands except the PR_deltpi bands
b2s=0.0;
b2f=0.0;

%%from RAL fits : note all were done at T = 296 K
pt=[2.6800e+01   1.6140e+02   5.6030e+02   9.6170e+02];  pt=pt/1013.25;
ps=[2.6800e+01   2.6800e+01   2.6800e+01   2.6800e+01];  ps=ps/1013.25;
pf=pt-ps;
ratio=pf./pt;

%use this d.of.c parameter everywhere!!!!!!!  from RAL fits
%%%%%%%%%%%%%%%%%%% BEST BEST BEST FOR RAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
dur=[10.46423654878965 6.573033770334661 5.614703762280506 5.251501801911941];
dur = [0.6             0.6               0.55              0.55];

%%%%%%%%%%%%%%%%%%% BEST BEST BEST FOR RAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%AIRS FUDGE FACTOR
airs_fudge_dur = 1.000;  dur = dur * airs_fudge_dur;
dur=dur*1e-3;
duration_self = dur(1);

if ((ptotal-pself)/ptotal < ratio(3))
  duration     = interp1(ratio,dur,(ptotal-pself)/ptotal);
else
  duration     = spline(ratio,dur,(ptotal-pself)/ptotal);
end

if ((ptotal-pself)/ptotal > 0.00001)
  duration_for = (ptotal*duration - pself*duration_self)/(ptotal-pself);
else
  duration_for = duration_self;
end


%%%%%%%%%%%%%%%%%%%% set up beta parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (band < 2200)    %%%%this sets up the P/Q/R branch SigPi, DeltPi mixing
  par_beta_self =   0.5358;            %self broadened     %%%%%used to be 0.5
  par_beta_for  =   0.6002;            %for  broadened     %%%%%used to be 0.62

  par_beta_self =   0.6000;            %self broadened     %%%%%used to be 0.5
  par_beta_for  =   0.5500;            %for  broadened     %%%%%used to be 0.62

end 
beta_self = par_beta_self;    %self broadened  
beta_for  = par_beta_for;     %foreign broadened
b2s = beta_self;
b2f = beta_for;

if ((band == 2150))
%%%%%%%%%%%%%%%%%%% BEST BEST BEST FOR RAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
  besp=[6.831266673446708 8.942573853443219 ...
        9.369380338865858 9.502017378026286];
%%%%%%%%%%%%%%%%%%% BEST BEST BEST FOR RAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
end

