if exist('mass.dat')
  load mass.dat        %get the mass of the isotopes
else
  disp('mass.dat DNE')
  fprintf(1,'this should have een made by     hitread.m     for your hitran version %s \n',hitran_version)
  fprintf(1,'so for eg H2020, see /umbc/xfs3/strow/asl/rta/hitran/H2020/ISOTOPES/Readme \n')
  fprintf(1,'                 then run convert2mass20.m \n')
  fprintf(1,'                 then  cp /asl/data/hitran/H2020/ISOTOPES/mass20.dat /home/sergio/SPECTRA/MASS_ISOTOPES/mass20.dat \n')
  disp(' ')
  fprintf(1,'so for eg H2024, see /umbc/xfs3/strow/asl/rta/hitran/H2024/ISOTOPES/Readme \n')
  fprintf(1,'                 then run convert2mass2X.m \n')
  fprintf(1,'                 then  cp /asl/data/hitran/H2020/ISOTOPES/mass24.dat /home/sergio/SPECTRA/MASS_ISOTOPES/mass24.dat \n')
  disp(' ')  
  fprintf(1,' the copy is in HITRAN_Database_important_scripts/2020/ISOTOPES \n')
  fprintf(1,' the copy is in HITRAN_Database_important_scripts/2024/ISOTOPES \n')
  disp(' ' )
  error('please fix!!!')
end  
%%pause(0.1);

[x,y] = size(mass);

if y == 5
  disp('   ** you have 5 columns in your mass.dat file ==> H2K or H04 or H08 or H12 or H16 or H20 or H24')
  disp('   ** data was used to save the nuclear degeneracy factor(unused even')
  disp('   ** by QTIPS???) and QTIPS partition fcns at 296K');
  %%% [gasID   mol.mass  abundance   deGeneracy    Q(296)]
  end

ladeda = mass(:,3);
iZero  = find(ladeda < -0.1);
iZero  = max(iZero);   %% this should tell where we have ended putting ZEROS
                       %% ie  where <iGasID  NumIsotopes  -1.0> ends
                       %% and where the <iGasID  MassIsotopes Abundance> begins

% for H92 .. H98, 
% first iNumGasIDS are the numbers of the isotopes  <iGasID  NumIsotopes  0.0>
% next bunch are the masses of the isotopes  <iGasID  MassIsotopes Abundance>
% for H2k and H04, this is not necessarily true; 
% dunno about H08
% there were 38,39 gases respectively, but truncated make_mass_HITRAN2MATLAB.f
% at the first 32 gases

iNumGasIDS = 42;
iNumGasIDS = 32;

if (iZero < iNumGasIDS)
  disp('Whoops : fix make_mass_HITRAN2MATLAB.f so at least 32 gases saved');
  error('either that or make load_mass_dat.m smarter!!!')
else
  iNumGasIDS = iZero;
end

% first iNumGasIDS are the numbers of the isotopes  <iGasID  NumIsotopes  0.0>
% next bunch are the masses of the isotopes  <iGasID  MassIsotopes Abundance>

massblock = mass(1:iNumGasIDS,1);
hoho      = find(massblock == gasID);
liso      = mass(hoho,2);               %%number of isotopes

datablock = mass(iNumGasIDS+1:x,1);
hoho      = find(datablock == gasID);
datablock = mass(iNumGasIDS+1:x,:);
mass_iso  = datablock(hoho,2);

if y == 5 | y == 6
  mass_abn   = datablock(hoho,3);
  mass_dgn   = datablock(hoho,4);
  mass_QT296 = datablock(hoho,5);
else
  mass_abn   = ones(liso,1);
  mass_dgn   = ones(liso,1);
  mass_QT296 = ones(liso,1);
end

%mass_iso  = round(mass_iso);  %% original, before June 2018
mass_info  = [mass_abn mass_dgn mass_QT296];
