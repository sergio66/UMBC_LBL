function LenCell = convert_input_run8prof_to_gas_cell_length(Pall,PPall,Tall,Qall);

% so now use /home/sergio/SPECTRA/set_c1_c2_avog_pref_tref.m for MGC
% so now use /home/sergio/SPECTRA/gas_cell_others.m to set params
%   pV = n R T ==> n/V = p/(RT) ==> L n/V = pL/RT ==> q = (partpress) L/(RT)
%   GasAmt = GasCellLen * partpress/1e9/MGC/temperature; %change to kmoles/cm2


%{
example

%% [sergio@chip-login2 MAKEIR_CO2_O3_N2O_CO_CH4_othergases_LBLRTM_v12.17_lnflv3.8.1]$ more load_refprof.m
load /home/sergio/HITRAN2UMBCLBL/REFPROF/refproTRUE.mat

gasID = 1;
Pall  = refpro.mpres;            %% atm * 1.01e5 atm--> N/m2
Qall  = refpro.gamnt(:,gasID);   %% kmol/cm2 * 1000 * 1e4      *1000 is kmol->mol, *1e4 is cm2->m2
PPall = refpro.gpart(:,gasID);   %% atm * 1.01e5 atm--> N/m2
Tall  = refpro.mtemp;            %% K
LenCell = convert_input_run8prof_to_gas_cell_length(Pall,PPall,Tall,Qall);

plot(LenCell/100/1000,'+-'); title('Cell length in Km')   %% /100 to m, /1000 to km

%}

MGC = 8.314674269981136;      %% J/mol/K

LenCell = (Qall * 1.0e7)./(PPall * 1.01e5) * MGC .* Tall; %% meters
LenCell = LenCell * 100;                                  %% cm
