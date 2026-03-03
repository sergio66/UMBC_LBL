%% this is common to the run8* set of codes

%     if which_isotopes = 0            use all isotopes
%     if which_isotopes = [-1 iaList]  throw away isotopes in iaList, keep rest
%     if which_isotopes = [iaList]     keep isotopes in iaList, throw rest
which_isotope = input('Enter which isotope \n  0 for all \n  [-1 iaX] to throw away isotopes in iaList, keep rest, \n  [iaX to only keep isotopes in list) \nEnter choice (default 0) : ');
if length(which_isotope) == 0
  which_isotope = 0;
end

do_load   = 0;
MinLayer  = 1; MaxLayer = 1; Step = 1;        %use only ONE "layer"
NumLayers = 1;
TheLayers = MinLayer:Step:MaxLayer;

press       = input('Enter total pressure (in atm) : ');
partpress   = input('Enter gas partial pressure (in atm) : ');
temperature = input('Enter temperature (in K) : ');

GasAmt      = input('Enter path cell length (in cm) ');
%change to kilomoles cm-2 
GasAmt      = GasAmt *PREF *partpress/1e9/MGC/temperature; %change to kmoles/cm2 
