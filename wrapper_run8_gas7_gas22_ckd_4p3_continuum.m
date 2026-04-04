function [fr43,od43] = wrapper_run8_gas7_gas22_ckd_4p3_continuum(gasID,outwave,profname);

fr43 = outwave;
profile = load(profname,'-ascii');
[mm,nn] = size(profile);

od43 = zeros(mm,length(outwave));

Pall  = profile(:,2);            %% atm * 1.01e5 atm--> N/m2
Qall  = profile(:,5);            %% kmol/cm2 * 1000 * 1e4      *1000 is kmol->mol, *1e4 is cm2->m2
PPall = profile(:,3);            %% atm * 1.01e5 atm--> N/m2
Tall  = profile(:,4);            %% K
LenCell = convert_input_run8prof_to_gas_cell_length(Pall,PPall,Tall,Qall);

n2_o2_continuum_range_ckd3p2_ckd4p3

if gasID == 22
  gasYESrange = n2range;
else
  gasYESrange = o2range;
end

addpath /home/sergio/HITRAN2UMBCLBL/MAKEIR/H2024/MAKEIR_CO2_O3_N2O_CO_CH4_othergases_LBLRTM_v12.17_lnflv3.8.1
addpath /home/sergio/SPECTRA/CKDLINUX

fmin = outwave(1);
fmax = outwave(end);
JOBB = 0;
for ll = 1 : mm
  fip = ['temp_n2ORo2_input_' num2str(ll)];
  fip = mktempS(fip);  
  fid = fopen(fip,'w');
  fxout = dowrite_n2_or_o2_new(fid,gasID,ll,fmin,fmax,JOBB,Pall*1013.25,Tall,PPall*1013.25,LenCell);
  fclose(fid);

  fprintf(1,'for CKD4.3 O2/N2 cntnm_sergio_v4.3_linux_gnu_dbl fip,fxout = %s %s \n',fip,fxout);
  % catter = ['!cat ' fip]; eval(catter);
  str = ['!../cntnm_sergio_v4.3_linux_gnu_dbl < ' fip];
  eval(str)
  junk = read_MT_CKD4p3_sergio(fxout);
  rmer = ['!/bin/rm ' fxout ' ' fip]; eval(rmer)

  w1 =  junk.data(:,1);
  od1 = junk.data(:,2);
  odx = interp1(w1,od1,outwave,[],'extrap');
  bad = find(outwave < gasYESrange(1) | outwave > gasYESrange(2)); odx(bad) = 0.0;
  od43(ll,:) = odx;
end
