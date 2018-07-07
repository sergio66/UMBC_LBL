function q = p_pp_T_L_to_q(press,partpress,temperature,cell_length)

MGC = 8.314674269981136  ;

%press       = input('Enter total pressure (in atm) : ');
%partpress   = input('Enter gas partial pressure (in atm) : ');
%temperature = input('Enter temperature (in K) : ');
%cell_Length = input('Enter path cell length (in cm) ');

%change to kmoles cm-2
GasAmt = cell_length*101325*partpress/1e9/MGC/temperature; %change to kmoles/cm2

fprintf(1,'GasAmt = %8.6e molecules/cm2 \n',GasAmt*6.023e26)
