%now do continuum .. it should be smooth enough that no need to do at
%fine mesh ... just do at output mesh
%       SUBROUTINE CALCON23( IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P, 
%     $    PARTP, AMNT, CON, CKD, whichlayer) 
%                    has been changed to
%      con = CALCON23( IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P, 
%     $    PARTP, AMNT, CKD, whichlayer)


docontinuum = -1;

if (gasID == 1)         %water
  if ((CKD == 0) | (CKD == 21) | (CKD == 23) | (CKD == 24))
    docontinuum = -1;
    fprintf(1,'\n this code WILL NOT DO water continuum ... \n');
    fprintf(1,'\n use run7water/ru76watercontinuum instead ... \n');
    error('exiting out of run7');
  end
end

if (gasID == 7)         %o2
  if (CKD == 1)
    docontinuum = 1;
    fprintf(1,'\n doing O2 continuum ... \n');
  end
end

if (gasID == 22)         %n2
  if (CKD == 1)
    docontinuum = 1;
    fprintf(1,'\n doing N2 continuum ... \n');
  end
end

%%%%%%%%%%%%%%%%% this is the orig code .... might blow up if input vectors
%%%%%%%%%%%%%%%%  are tooo large
%if (docontinuum > 0)
%    for jj=MinLayer:Step:MaxLayer %INNER LOOP 1..100 = bottom -> top
%      nn=(jj-MinLayer)/Step+1;
%      tempjunk=calcon(gasID,length(outwave),outwave,1,length(press),...
%                    temperature,press,partpress,GasAmt,CKD,jj);
%      out_array(nn,:)=out_array(nn,:)+tempjunk;
%    end
%end      

%******************************************************

%%%%from FORTRANFILES/max.inc
%c this is max length of arrays that can be used in the Mex Files  
%c this number came out of 
%c   200000 = max number of elements in mesh 
%c              eg (755-655)/0.0005 = 160000
%c        4 = number tacked on to arrays so boxint(y,5) can be done  
%      integer MaxLen
%      parameter(MaxLen=200010)
MaxLen = 200010;

%% docontinuum  =  docontinuum*2   %%% this way, only pass out CONTINUUM

dff = ffin*nbox;

%% to figurre out V1 V2 V3 V4 go to CKDLINUX/calcon.F
%% and uncomment the relevant iVers+includefile, compile with 
%% makewatermexLinux, then 
%%    mv  calcon.mexa64 calconVN.mexa64
%%

if (docontinuum == 1 & O2O3N2continuumVers < 0)
  fprintf(1,'GasID = %3i : needs continuum \n',gasID);
  disp('but you have set O2O3N2continuumVers = -1 so NO continuum added on!')

elseif (docontinuum == 1 & O2O3N2continuumVers > 0)
  disp('N2 or O2 : adding on continuum OD')
  [mm0,nn0] = size(out_array);
  if O2O3N2continuumVers == 5 & nn0 <= 10001
     homedirx = pwd;
     cd /home/sergio//SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/
     [fr5,od5] = wrapper_run8_gas7_gas22_continuum(gasID,fmin,fmax,0.0005,profname);
     out_array = out_array + od5;
     cder = ['cd ' homedirx];
     eval(cder);
  elseif O2O3N2continuumVers == 5 & nn0 > 10001
     error('ooops need to increase array sizes for wrapper_run8_gas7_gas22_continuum MTCKD3.2')
     [fr5,od5] = wrapper_run8_gas7_gas22_continuum(gasID,fmin,fmax,0.0005,profname);
     out_array = out_array + od5;
  elseif (O2O3N2continuumVers ~= 5) & (length(outwave) <= MaxLen)          %can do it in one big chunk!
    for jj=MinLayer:Step:MaxLayer %INNER LOOP 1..100 = bottom -> top
      nn=(jj-MinLayer)/Step+1;
      switch O2O3N2continuumVers
        case 1
          tempjunk = calconV1(gasID,length(outwave),outwave,1,length(press),...
                  temperature,press,partpress,GasAmt,CKD,jj);
        case 3	  
          tempjunk = calconV3(gasID,length(outwave),outwave,1,length(press),...
                  temperature,press,partpress,GasAmt,CKD,jj);
        case 4
          tempjunk = calconV4(gasID,length(outwave),outwave,1,length(press),...
                  temperature,press,partpress,GasAmt,CKD,jj);
        case 5
          tempjunk = calconV_MTCKD3p2(gasID,length(outwave),outwave,1,length(press),...
                  temperature,press,partpress,GasAmt,CKD,jj);
        otherwise
          error('only have o2/n2/o3 continuum Vers 1,3,4,5');
      end
      out_array(nn,:) = out_array(nn,:)+tempjunk;
    end     %% for loop

  elseif O2O3N2continuumVers ~= 5 & (length(outwave) > MaxLen) %have to break it into frequency intervals
    for jj = MinLayer:Step:MaxLayer %LOOP 1..100 = bottom -> top
      nn = (jj-MinLayer)/Step+1;
      index0 = 1:10000;
      number = floor(length(outwave)/10000);

      for kk = 1:number             %LOOP OVER FREQUENCY INTERVALS
        index = index0+(kk-1)*10000;
        switch O2O3N2continuumVers
          case 1
            tempjunk = calconV1(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 3
            tempjunk = calconV3(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 4
            tempjunk = calconV4(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 5
            tempjunk = calconV_MTCKD3p2(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          otherwise
            error('only have o2/n2/o3 continuum Vers 1,3,4,5');
      end

        out_array(nn,index) = out_array(nn,index)+tempjunk;
      end

      %%see if anything left over
      if (index(length(index))+1 <= length(outwave)) 
        index = index(length(index))+1:length(outwave);
        switch O2O3N2continuumVers
          case 1
            tempjunk = calconV1(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 3
            tempjunk = calconV3(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 4
            tempjunk = calconV4(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          case 5
            tempjunk = calconV_MTCKD3p2(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
          otherwise
            error('only have o2/n2/o3 continuum Vers 1,3,4,5');
      end
        out_array(nn,index) = out_array(nn,index)+tempjunk;
      end

    end
  end            %%%%if then else loop

elseif (docontinuum == 2)
  disp('N2 or O2 : replacing OD with continuum OD')
  error('why bother with this???????')
end

