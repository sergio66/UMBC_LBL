function q = gascell2ip(varargin)

%% see run8.m

set_c1_c2_avog_pref_tref

if nargin == 0
  press       = input('Enter total pressure (in atm) : ');
  partpress   = input('Enter gas partial pressure (in atm) : ');
  temperature = input('Enter temperature (in K) : ');
  GasAmt      = input('Enter path cell length (in cm) ');
elseif nargin == 4
  press = varargin{1};
  partpress = varargin{2};
  temperature = varargin{3};
  GasAmt = varargin{4};
else
  error('need 0 arguments (interactive) or 4 arguments (p[atm],ps[atm],T[K],L[cm])')
end

%change to kmoles cm-2
GasAmt = GasAmt * PREF * partpress/1e9/MGC/temperature; %change to kmoles/cm2
q = GasAmt;

fprintf(1,'GasAmt = %10.6e  kmoles/cm2 \n',GasAmt)
