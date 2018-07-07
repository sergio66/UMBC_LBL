function [absc, vgrid] = calc_xsec_UV_li_zhu(gf, v1, v2, dv, tp, pL);

% function [absc, vgrid] = calc_xsec_UV_li_shu(gf, v1, v2, dv, tp, pL);
%
% calc_xsec returns ozone UV  absorptions calculated by
% interpolation or extrapolation of tabulated data
%
% inputs
%   gf  - gas id
%   v1  - lower wavenumber
%   v2  - upper wavenumber
%   dv  - wavenumber increment
%   tp  - temperature profile
%   pL  - pressure levels for tp .... assumes O3 PARTIAL PRESSURE in atm
%
% implicit inputs
%   nvpts = (v2-v1)/dv + 1, number of wavenumber points
%   nplev = length(pL), number of pressure levels
%
% outputs
%   absc  - nvpts x nplev array of layer absorptions, kmoles/cm^2
%   vgrid - nvpts vector of output frequencies
%

% Original Version:
% - H. Motteler, 15 Dec 98
% Revisions:
% - interpolation bug fix; fancier plotting; removed xsec data 
%   directory param; renamed "calc_xsec"; HM, 4 Apr 00

% physical extrapolation parameter, Ke = Km*(Tm/Te)^C1, where 
%   - Km is measured absorption at temperature Tm, and 
%   - Ke is estimated absorption at temperature Te.
%
C1 = 0.75;

% input defaults
if (nargin < 6)
  error('wrong number of arguments')
  end

%% see /home/sergio/SPECTRA/findxsec_plot_UV.m
%% supplied by Li Zhu 2008, gives cross sections in units of /(atm cm) 
hitlin_fname = '/home/sergio/SPECTRA/VISIBLE_OD/O3.csv'; 
%% recall pV = nRT ==> n/V = 2.4486e+19 molecules/cm3 at 296 K 
%%        pV = nRT ==> n/V = 2.6549e+19 molecules/cm3 at 273 K 
%% so this file says cross section sigmaLI = 267.063 (1/atm cm) at 245.1 nm 
%% sigmaLI(1/atm cm) = sigmaHIT(cm2/molecule) * (molecules/cm3)/atm 
%%                   = 1e-20 * 9.9079E2 * 2.6549e+19 
liO3 = load(hitlin_fname); 
liO3(:,1) = liO3(:,1)/10;     %% change from Angstroms to nm 
liO3(:,1) = liO3(:,1)/1000;   %% change from nm to um 
liO3(:,1) = 10000./liO3(:,1); %% change from um to wavenumber  
v1x = min(liO3(:,1)); 
v2x = max(liO3(:,1)); 

% initializations
nvpts = 1 + round((v2-v1)/dv);  % points in output grid
vgrid = v1 + (0:nvpts-1) * dv;  % output wavenumber grid
v2 = vgrid(nvpts);              % make v2 a multiple of dv
nplev = length(pL); 	        % number of pressure levels
absc = zeros(nvpts, nplev);     % return zeros by default

c0 = liO3(:,2); c0 = interp1(liO3(:,1),c0,vgrid);
c1 = liO3(:,3); c1 = interp1(liO3(:,1),c1,vgrid);
c2 = liO3(:,4); c2 = interp1(liO3(:,1),c2,vgrid);

boo = find(vgrid < v1x | vgrid > v2x);
if length(boo) > 0
  c0(boo) = 0.0;
  c1(boo) = 0.0;
  c2(boo) = 0.0;
  end

T0 = 273.15;
for ll = 1 : nplev
  Tc = tp(ll) - T0;
  %% abscoef is in 1/(atm cm) while pL is in atm, so output absc is in cm-1
  absc(:,ll) = ( c0 + c1*Tc + c2*Tc*Tc)*pL(ll);   
  end