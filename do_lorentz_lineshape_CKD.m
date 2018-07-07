function [out_array] = do_lorentz_lineshape_CKD(outwave,out_array,AVOG,c2,...
                        temperature,press,partpress, GasAmt,...
                        CKD,CKD_0,selfmult,formult,profname,...
                        local,divide,ffin,nbox,MaxLen,MinLayer,Step,MaxLayer);

%% used by run8watercontinuum.m

if (ismember(CKD,[1 2 3 4 5 6 24 25 50 51 52 53 55 56 60]))
  error('v 1,2,3,4,5,6,24,25,50-56,60 impossible with lorentz lineshape!!');
end

if (ismember(CKD,[0 21 23]))
  fprintf(1,'\n CKD continuum %3i Lorentz lineshape = %2i \n',CKD,local);
  if (length(outwave) <= MaxLen)          %can do it in one big chunk!
    for jj=MinLayer:Step:MaxLayer %INNER LOOP 1..100 = bottom -> top
      nn=(jj-MinLayer)/Step+1;
      scum=calconwater(1,length(outwave),outwave,dff,length(press),...
          temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
      if (divide  == -1)
        tempfreq=ones(size(scum));
      else
        tempfreq = AVOG*GasAmt(jj)*outwave .* ...
                     tanh(c2*outwave/2/temperature(jj))* ...
                     (296/temperature(jj));
        if ((selfmult >= 0.999999) & (formult <= 0.00000001))
          tempfreq=tempfreq * partpress(jj);
        elseif ((formult >= 0.999999) & (selfmult <= 0.00000001))
          tempfreq=tempfreq * (press(jj)-partpress(jj));
        end
      end
      if max(tempfreq) <= 1e-20
        scum=scum./(tempfreq+1e-20);   %% else we get divide by 0 = NAN
      else
        scum=scum./(tempfreq);
      end
      out_array(nn,:)=out_array(nn,:)+scum;
    end

  else %have to break it into frequency intervals
    fprintf(1,'breaking into smaller intervals, length (MaxLen-10) \n');
    for jj=MinLayer:Step:MaxLayer %LOOP 1..100 = bottom -> top
      nn=(jj-MinLayer)/Step+1;
      index0=1:(MaxLen-10);
      number=floor(length(outwave)/(MaxLen-10));
      for kk=1:number             %LOOP OVER FREQUENCY INTERVALS
        index=index0+(kk-1)*(MaxLen-10);
        scum=calconwater(1,length(index),outwave(index),dff,length(press),...
            temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
        if (divide  == -1)
          tempfreq=ones(size(scum));
        else
          tempfreq = AVOG*GasAmt(jj)*outwave(index) .* ...
                       tanh(c2*outwave(index)/2/temperature(jj))* ...
                       (296/temperature(jj));
          if ((selfmult >= 0.999999) & (formult <= 0.00000001))
            tempfreq=tempfreq * partpress(jj);
          elseif ((formult >= 0.999999) & (selfmult <= 0.00000001))
            tempfreq=tempfreq * (press(jj)-partpress(jj));
          end
        end
        if max(tempfreq) <= 1e-20
          scum=scum./(tempfreq+1e-20);   %% else we get divide by 0 = NAN
        else
          scum=scum./(tempfreq);
        end
        out_array(nn,index)=out_array(nn,index)+scum;
      end
       %%see if anything left over
      if (index(length(index))+1 <= length(outwave)) 
        index=index(length(index))+1:length(outwave);
        scum=calconwater(1,length(index),outwave(index),dff,length(press),...
                temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
        if (divide  == -1)
          tempfreq=ones(size(scum));
        else
          tempfreq = AVOG*GasAmt(jj)*outwave(index) .* ...
                       tanh(c2*outwave(index)/2/temperature(jj))* ...
                       (296/temperature(jj));
          if ((selfmult >= 0.999999) & (formult <= 0.00000001))
            tempfreq=tempfreq * partpress(jj);
          elseif ((formult >= 0.999999) & (selfmult <= 0.00000001))
            tempfreq=tempfreq * (press(jj)-partpress(jj));
          end
        end
        if max(tempfreq) <= 1e-20
          scum=scum./(tempfreq+1e-20);   %% else we get divide by 0 = NAN
        else
          scum=scum./(tempfreq);
        end
        out_array(nn,index)=out_array(nn,index)+scum;
      end

    end
  end            %%%%if then else loop
end
