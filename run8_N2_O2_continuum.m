%now do continuum .. it should be smooth enough that no need to do at
%fine mesh ... just do at output mesh
%       SUBROUTINE CALCON23( IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P, 
%     $    PARTP, AMNT, CON, CKD, whichlayer) 
%                    has been changed to
%      con = CALCON23( IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P, 
%     $    PARTP, AMNT, CKD, whichlayer)

docontinuum = -1;

if (gasID == 1)         %water
  docontinuum = -1;
  fprintf(1,'\n this code WILL NOT DO water continuum ... \n');
  fprintf(1,'\n use run8water/run8watercontinuum instead ... \n');
  error('exiting out of run8');
end

if (gasID == 7)         %o2
  docontinuum = 1;
  fprintf(1,'\n doing O2 continuum ... \n');
end

if (gasID == 22)         %n2
  docontinuum = 1;
  fprintf(1,'\n doing N2/02 continuum ... \n');
end

if docontinuum < 0
  fprintf(1,'run8_N2_O2_continuum.m : not doing continuum for gasID %3i \n',gasID)
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (docontinuum == 1 & O2O3N2continuumVers < 0)
  fprintf(1,'GasID = %3i : needs continuum \n',gasID);
  disp('but you have set O2O3N2continuumVers = -1 so NO continuum added on!')

elseif (docontinuum == 1 & O2O3N2continuumVers > 0)
  
  disp('N2 or O2 : adding on continuum OD')
  [mm0,nn0] = size(out_array);
  out_array0 = out_array;

  iCheckOld = -1;
  if iCheckOld > 0
    plot(outwave,out_array)
    cd /home/sergio/SPECTRA
    run8_N2_O2_continuum_beforeApr2026
    od32cont = out_array - out_array0;
  end
  
  homedirx = pwd;
  cd /home/sergio//SPECTRA/CKDLINUX/MT_CKD-4.3/run_example
  [fr43,od43cont] = wrapper_run8_gas7_gas22_ckd_4p3_continuum(gasID,outwave,profname);
  if iCheckOld > 0  
    semilogy(outwave,od32cont,'b',outwave,od43cont,'r')
  end
  
  out_array = out_array0;
  out_array = out_array + od43cont;
  cder = ['cd ' homedirx];
  eval(cder);

elseif (docontinuum == 2)
  disp('N2 or O2 : replacing OD with continuum OD')
  error('why bother with this???????')
end

