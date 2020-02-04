function q = gascell2ip(varargin)
%% see run8.m

MGC = 8.314674269981136  ;

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
GasAmt = GasAmt*101325*partpress/1e9/MGC/temperature; %change to kmoles/cm2
q = GasAmt;

fprintf(1,'GasAmt = %10.6e  kmoles/cm2 \n',GasAmt)
