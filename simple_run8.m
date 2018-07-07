function [outvectV,outvectD,outvectL] = ...
  simple_run8(wavenumberIN,mol_mass,line,qfcn,iWhich,profile)

%% function [outvectV,outvectD,outvectL] = ...
%%  simple_run8(wavenumberIN,mol_mass,line,qfcn,iWhich,profile)
%%
%% cut and paste heavily from run8 to compute spectra corresponding to ONE line
%% input 
%%   wavenumberIN = wavenumber array (in cm-1) over which spectra is computed
%%       mol_mass = molecular mass (in amu) eg 12+16 = 28 for CO
%%       qfcn     = partition function
%%       iWhich   = which of below lines to use
%%       line     = HITRAN spectral database info, should have following fields
%%                  line.linct  
%%                  line.igas   
%%                  line.iso    
%%                  line.wnum   
%%                  line.stren  
%%                  line.tprob  
%%                  line.abroad 
%%                  line.sbroad 
%%                  line.els    
%%                  line.abcoef 
%%                  line.tsp    
%%                  line.iusgq  
%%                  line.ilsgq  
%%                  line.gasid  
%%       profile  = [gasID p pp T q]   p,pp in atm; T in K; q in kilomoles/cm2
%%             
%% output
%%   voigt, doppler and lorentz spectral lineshapes computed using above info
%%   
%% Author : Sergio DeSouza-Machado, UMBC Nov 2011
%%
%%%%%%%%%%%%%%%%%%% IS THIS INTERACTIVE SESSION OR CLUNK THRU PROFILE %%%%%%%%%
useruser = -1;
if (useruser > 0)
  MGC=8.314674269981136  ;  
  gasID       = input('Enter HITRAN gasID ');
  press       = input('Enter total pressure (in atm) : ');
  partpress   = input('Enter gas partial pressure (in atm) : ');
  temperature = input('Enter temperature (in K) : ');
  GasAmt      = input('Enter path cell length (in cm) ');
  %change to kmoles cm-2 
  GasAmt = GasAmt*101325*partpress/1e9/MGC/temperature; %change to kmoles/cm2 
else
  gasID       = profile(1);
  press       = profile(2);
  partpress   = profile(3);
  temperature = profile(4);
  GasAmt      = profile(5);
  end

tempr     = temperature;
p         = press;
ps        = partpress;

pwr      = line.abcoef(iWhich);
for_brd  = line.abroad(iWhich);
self_brd = line.sbroad(iWhich);
freq     = line.wnum(iWhich) + press*line.tsp(iWhich);  %freq shift
energy   = line.els(iWhich);
s0       = line.stren(iWhich);
brd      = broad(p,ps,1.0,for_brd,self_brd,pwr,tempr,gasID);
strength = find_stren(qfcn(iWhich),freq,tempr,energy,s0,GasAmt);

outvectV = strength*voigt1_mat(wavenumberIN,freq,tempr,mol_mass,brd);
outvectL = strength*lorentz(wavenumberIN,freq,tempr,mol_mass,brd);
outvectD = strength*line_doppler(wavenumberIN,freq,tempr,mol_mass,brd);

outvectV = outvectV';

%semilogy(wavenumberIN,[outvectL; outvectD; outvectV]);
%hl = legend('lorentz','doppler','voigt');

%%%%%%%%%%%%%%%% you guys ignore this!
%% gas_cell_others_ppmv                                 
%%    Enter pressure units (1) atm (2) torr (3) mb : 1 
%%    Enter gasID 5 
%%    Enter [start stop] wavenumber : [2105 2205] 
%%    Enter total pressure  : 1 
%%    Enter gas ppmv  : 1000/1000      %% 1000 ppb = 1 ppmv 
%%    Enter temperature (in K) : 285 
%%    Enter path cell length (in cm) 4000 
%%
%% save qfcn and line in co_lines.mat
%% [wn,y] = run8(5,2105,2205,'IPFILES/gas_cell');  
%%
%% topts.LVG = 'G';
%% actually run it in 10 sec
%% [wn,y] = run8(5,2105,2205,'IPFILES/gas_cell',topts);  
%% load co_lines.mat    --- load in lines, qfc
%% save co_lines.mat    --- save EVERYTHING
%% profile = [5    1.00000000    0.00000100   285.000   1.71036e-10]
%%
%% [outvect] = simple_run8(wn,28,line,qfcn,184,profile);             
%%  plot(wn,y,wn,outvect)


