function lineX = adjust_for_uncertainty(line,adjust_str,gid,fr0,iDebug)

% there are 6 values that can be changed :
%     wavenumber, intensity, air- and self-
%     broadened halfwidths, temperature-dependence, and pressure shift 

% input
%   line       : HITRAN parameters straight from hitread
%   adjust_str : how to change the 6 parameters + for max, - for min, R for random, else (eg X) unchanged
%     wnum_unc_index   = str2num(unc_index(:,1));    line center
%     stren_unc_index  = str2num(unc_index(:,2));    line strength
%     abroad_unc_index = str2num(unc_index(:,3));    air broadening
%     sbroad_unc_index = str2num(unc_index(:,4));    self broadening
%     abcoef_unc_index = str2num(unc_index(:,5));    temp dependence of broadening
%     tsp_unc_index    = str2num(unc_index(:,6));    pressure shift
%          gid : GasID            ... needed for random number seed
%          fr0 : start wavenumber ... needed for random number seed
%       iDebug : optional argument
%
% output
%   lineX      : output

% example adjust_str = []       : nothing changed
%                    = '      ' : nothing changed
%                    = '000000' : nothing changed
%                    = 'XXXXXX' : nothing changed
%                    = '+XXXXX' : wavenumber changed to v+1dv, nothing else changed
%                    = '-XXXXX' : wavenumber changed to v-1dv, nothing else changed
%                    = 'RXXXXX' : wavenumber changed to v+Rdv, nothing else changed
%                                 where R is a random number between -1 and +1
%                    = 'RRRRRR' : wavenumber changed to v+Rdv, as is broadening, shift, line strength
%                                 where R is a random number between -1 and +1

%{
https://www.cfa.harvard.edu/hitran/uncertainty.html

Uncertainty Indices

The following table gives the definition of the uncertainty indices
used in the HITRAN database. These codes are found towards the end of
each line record (consult general format). The code for accuracy of
line position and the pressure shift of the line is fairly
straightforward: the code indicates the accuracy of the figures after
the decimal in wavenumber (hence there is no need to "write" the line
position with more digits than that). For the other parameters
(intensity, air-broadened half-width, self-broadened half-width, and
the coefficient for the temperature-dependence of the air-broadened
half-width) we have adopted a different but uniform coding.

Uncertainty Codes Used in HITRAN Database

Wavenumber and                               Intensity, Halfwidths, and
Pressure shift (cm-1)                        Temperature-dependence

Code            Uncertainty Range            Code        Uncertainty Range

0              ≥ 1. or Unreported              0         Unreported or Unavailable
1              ≥ 0.1 and < 1.                  1         Default or Constant
2              ≥ 0.01 and < 0.1                2         Average or Estimate
3              ≥ 0.001 and < 0.01              3         ≥ 20%
4              ≥ 0.0001 and < 0.001            4         ≥ 10% and < 20%
5              ≥ 0.00001 and < 0.0001          5         ≥ 5% and < 10%
6              ≥ 0.000001 and < 0.00001        6         ≥ 2% and < 5%
7              ≥ 0.0000001 and < 0.000001      7         ≥ 1% and < 2%
8              ≥ 0.00000001 and < 0.0000001    8         < 1%
9              Better than 0.00000001

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

also see poster by P_I_11_Jacquemart_D.pdf
Jacquemart and Gamache
http://imk-asf.webarchiv.kit.edu/sat/ClosedProjects/assfts/P_I_11_Jacquemart_D.pdf
(saved in HITRAN_PDFs)

[PDF]The Next HITRAN Edition: Description of new Parameters and Formats
imk-asf.webarchiv.kit.edu/sat/ClosedProjects/assfts/P_I_11_Jacquemart_D.pdf
by D Jacquemart - ‎Related articles

incumbent upon users to consult the references. Sources for these
parameters are provided by the reference indices (correspondence to
these references can be found in ref-table2003.pdf available at
ftp://cfa- ftp.harvard.edu/pub/HITRAN/Global_Data/). Links are
provided to many of the abstracts. <1%. 8. ≥1% and <2%. 7.

see page 40 of https://modis-images.gsfc.nasa.gov/JavaHAWKS/hawksman.pdf
FORTRAN Format (I2,I1,F12.6,1P2E10.3,0P2F5.4,F10.4,F4.2,F8.6,2A15,2A15,6I1,6I2,A1,2F7.1)
corresponding to:

6I1 : Uncertainty indices for wavenumber, intensity, air- and self-
broadened halfwidths, temperature-dependence, and pressure shift
(saved in HITRAN_PDFs)


JavaHAWKS Manual
36
Example of 100-character HITRAN line-transition format.
FORTRAN Format (I2,I1,F12.6,1P2E10.3,0P2F5.4,F10.4,F4.2,F8.6,2I3,2A9,3I1,3I2)
corresponding to:
Mol       I2      Molecule number
Iso       I1      Isotopologue number (1= most abundant, 2= second most abundant, etc.)
νij       F12.6   Wavenumber in cm-1
Sij       E10.3   Intensity in cm-1/(molecule x cm-2) @ 296K
Rij       E10.3   Weighted transition mo ment-squared in Debyes
γair      F5.4    Air-broadened halfwidth (HWHM) in cm-1/atm @ 296K
γself     F5.4    Self-broadened halfwidth (HWHM) in cm-1/atm @ 296K
E′′       F10.4   Lower state energy in cm-1
n         F4.2    Coefficient of temperature dependence of air-broadened halfwidth
δ         F8.6    Air-broadened pressure shift of line transition in cm-1/atm @ 296K
iv′, iv′′ 2I3     Upper-state global quanta index,lower-state global quanta indices
q′, q"    2A9     Upper-state local quanta, lower-state local quanta
ierr      3I1     Uncertainty indices for wavenumber, intensity, and air-broadened  halfwidth
iref      3I2     Indices for table of references corresponding to wavenumber, 
                  intensity, and halfwidth

Example of 160-character HITRAN line-transition format. 
FORTRAN Format (I2,I1,F12.6,1P2E10.3,0P2F5.4,F10.4,F4.2,F8.6,2A15,2A15,6I1,6I2,A1,2F7.1) 
corresponding to:
Mol        I2        Molecule  number        
Iso        I1        Isotopologue number (1= most abundant, 2= second most abundant, etc.) 
νij        F12.6     Wavenumber         in         cm-1
Sij        E10.3     Intensity in cm-1/(molecule x cm-2) @ 296K 
Aij        E10.3     Einstein-A         coefficient         
γair       F5.4      Air-broadened halfwidth (HWHM) in cm-1/atm @ 296K 
γself      F5.4      Self-broadened halfwidth (HWHM) in cm-1/atm @ 296K
E′′        F10.4     Lower state energy in cm-1
n          F4.2      Coefficient of temperature dependence of air-broadened halfwidth 
δ          F8.6      Air-broadened pressure shift of line transition in cm-1/atm @ 296K 
v′, v′′    2A15      Upper-state global quanta, lower-state global quanta  
q′, q′′    2A15      Upper-state local quanta, lower-state local quanta 
ierr       6I1       Uncertainty indices for wavenumber, intensity, air- and self-broadened halfwidths, temperature-dependence, and pressure shift 
iref       6I2       Indices for table of references corresponding to wavenumber, 
                     intensity, air- and self-broadened halfwidths, temeperature-
                     dependence, and pressure shift
Flag       A1        Flag (*) for lines supplied with line-coupling algorithm 
g′         F7.1      Upper-state statistical weight
g′′        F7.1      Lower-state statistical weight 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dec 6, 2017 at 1:25 PM,

Dear Sergio,

In the case of the shift it is pretty simple, in most of the cases
when the uncertainty is zero it usually means that we just do not have
a value and the default is zero. For most of the atmospheric gases at
terrestrial temperatures the shifts range from -.02 to +0.005
cm-1. Therefore zero values are not that far off and you can
conservatively use +/- 0.02 cm-1 as the uncertainty in that particular
case.

Things are slightly different for the wavenumber as in that case, the
values either come from pre HITRAN1986 sources or in the cases when
there is no information available at all (we have very few cases like
that). I do not see any of these values deviating from real numbers by
more than 1 cm-1 in the worst cases, moreover if we are just talking
pre-HITRAN1986 era then I would even say better than 0.01 cm-1. To
play it safe you can set this to 1 cm-1, but it will be an
overestimation of uncertainty on average. Good news is that we do not
have many cases like that at all for the dominant absorbers so this
should be Ok.

I hope this helps,
Regards,
Iouli

%%%%%%%%%%%%%%%%%%%%%%%%%

December 15, 2017

Hi Marco,

I just talked to Iouli Gordon. For both line center and pressure
shift, when the uncertainty index is "0" he suggested using a max of
0.01 cm-1 instead (he said the number "0" was inserted in the H2008
database, but since then there have been more reliable measurements
whose uncertainty has not yet made it into the H201X databases, so
that is a safe number)

As far as linestrength and broadenings, he said when the index is 0,1,2 just use 20%

It'll take me a while to re-run my stuff (I may stop including HNO3,
so many lines!) but hopefully early next year it should be done.

He also recommended trying the speed dependent lineshapes in the
future, plus he said the line mixing code they have provided is in
good shape and works well with the H2016 database.

The new CrIS instrument is coming online soon, so trying this may be
pushed down the order for a while.

Cheers

Sergio

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 4
  iDebug = +1;
  iDebug = -1;    %% turn off debug
end  

%% see  /home/sergio/SPECTRA/read_hitr06/read_hitran.c
lineX = line;

if length(adjust_str) == 0
  disp('no need to adjust line strength parameeters : adjust_str is null string')
  return
elseif length(adjust_str) < 6
  adjust_str(length(adjust_str)+1:6) = 'X';
end

unc_index = lineX.ai;
  wnum_unc_index   = str2num(unc_index(:,1));
  stren_unc_index  = str2num(unc_index(:,2));
  abroad_unc_index = str2num(unc_index(:,3));
  sbroad_unc_index = str2num(unc_index(:,4));
  abcoef_unc_index = str2num(unc_index(:,5));
  tsp_unc_index    = str2num(unc_index(:,6));
  
for ii = 1 : 6
  uncvalue = adjust_str(ii);
  if uncvalue == '+'
    uncvalue = +1;    %% set to Y + max(dY)
  elseif uncvalue == '-'
    uncvalue = -1;    %% set to Y - max(dY) 
  elseif uncvalue == 'R' | uncvalue == 'r'
    uncvalue = +10;   %% set to Y + random(dY)
  else
    uncvalue = -10;   %% do nothing, keep as Y
  end

  if ii == 1
    unctype = 1; str = 'wnum'; unique(wnum_unc_index)'
    x_with_unc = find_unc(unctype,uncvalue,lineX.wnum,wnum_unc_index,gid,fr0,iDebug);
    lineX.wnum = x_with_unc;
    figure(1); clf; plot(line.wnum); 
    figure(2); clf; plot(wnum_unc_index); 
    figure(3); clf; plot(lineX.wnum);
    figure(4); clf; plot(lineX.wnum-line.wnum);
    figure(5); clf; plot(lineX.wnum./line.wnum);
    figure(6); clf; hist(lineX.wnum-line.wnum);
  elseif ii == 2
    unctype = 2; str = 'stren'; unique(stren_unc_index)' 
    x_with_unc = find_unc(unctype,uncvalue,lineX.stren,stren_unc_index,gid,fr0,iDebug);
    lineX.stren = x_with_unc;
    figure(1); clf; plot(line.stren); 
    figure(2); clf; plot(stren_unc_index); 
    figure(3); clf; plot(lineX.stren);
    figure(4); clf; plot(lineX.stren-line.stren);
    figure(5); clf; plot(lineX.stren./line.stren);
    figure(6); clf; hist(lineX.stren-line.stren);
  elseif ii == 3
    unctype = 2; str = 'abroad'; unique(abroad_unc_index)' 
    x_with_unc = find_unc(unctype,uncvalue,lineX.abroad,abroad_unc_index,gid,fr0,iDebug);
    lineX.abroad = x_with_unc;
    figure(1); clf; plot(line.abroad); 
    figure(2); clf; plot(abroad_unc_index); 
    figure(3); clf; plot(lineX.abroad);
    figure(4); clf; plot(lineX.abroad-line.abroad);
    figure(5); clf; plot(lineX.abroad./line.abroad);
    figure(6); clf; hist(lineX.abroad-line.abroad);
  elseif ii == 4
    unctype = 2; str = 'sbroad'; unique(sbroad_unc_index)' 
    x_with_unc = find_unc(unctype,uncvalue,lineX.sbroad,sbroad_unc_index,gid,fr0,iDebug);
    lineX.sbroad = x_with_unc;
    figure(1); clf; plot(line.sbroad); 
    figure(2); clf; plot(sbroad_unc_index); 
    figure(3); clf; plot(lineX.sbroad);
    figure(4); clf; plot(lineX.sbroad-line.sbroad);
    figure(5); clf; plot(lineX.sbroad./line.sbroad);
    figure(6); clf; hist(lineX.sbroad-line.sbroad);
  elseif ii == 5
    unctype = 2; str = 'abcoef'; unique(abcoef_unc_index)' 
    x_with_unc = find_unc(unctype,uncvalue,lineX.abcoef,abcoef_unc_index,gid,fr0,iDebug);
    lineX.abcoef = x_with_unc;
    figure(1); clf; plot(line.abcoef); 
    figure(2); clf; plot(abcoef_unc_index); 
    figure(3); clf; plot(lineX.abcoef);
    figure(4); clf; plot(lineX.abcoef-line.abcoef);
    figure(5); clf; plot(lineX.abcoef./line.abcoef);
    figure(6); clf; hist(lineX.abcoef-line.abcoef);
  elseif ii == 6
    unctype = -1; str = 'tsp'; unique(tsp_unc_index)' 
    x_with_unc = find_unc(unctype,uncvalue,lineX.tsp,tsp_unc_index,gid,fr0,iDebug);
    lineX.tsp = x_with_unc;
    figure(1); clf; plot(line.tsp); 
    figure(2); clf; plot(tsp_unc_index); 
    figure(3); clf; plot(lineX.tsp);
    figure(4); clf; plot(lineX.tsp-line.tsp);
    figure(5); clf; plot(lineX.tsp./line.tsp);
    figure(6); clf; hist(lineX.tsp-line.tsp);
  end
  figure(1); title([str ' orig']);      figure(2); title([str ' uncertainty index']);
  figure(3); title([str ' new']);       figure(4); title([str ' new-orig']);
  figure(5); title([str ' new/orig']);  figure(6); title([str ' hist(new-orig)']); 
  if iDebug == +1
    disp('ret to continue'); pause
  else
    pause(1)
  end
end
