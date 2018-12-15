function [out_array] = do_local_lineshape_CKD(outwave,out_array,AVOG,c2,...
                        temperature,press,partpress, GasAmt,...
                        CKD,CKD_0,selfmult,formult,profname,...
                        local,divide,ffin,nbox,MaxLen,MinLayer,Step,MaxLayer);

%% used by run8watercontinuum.m

MaxLenX = 50000;
if (ismember(CKD,[0 21 23 24   1 2 3 4 5 6   25 27 32]))
  fprintf(1,'\n doing CKD continuum %3i lineshape = %2i \n',CKD,local);
  haha = [outwave(1) outwave(end) mean(diff(outwave))];
  fprintf(1,'  outwave(1) outwave(end) diff(outwave) = %20.12f %20.12f %20.12f \n',haha);
  fprintf(1,'  length(outwave)  MaxLen = %8i %8i \n',length(outwave),MaxLen)  
  dff=ffin*nbox;

  if (length(outwave) <= MaxLenX)          %can do it in one big chunk!
    for jj=MinLayer:Step:MaxLayer %INNER LOOP 1..100 = bottom -> top
      nn=(jj-MinLayer)/Step+1;
      if (CKD_0 == 50)
        scum = mst50(outwave,temperature,press,partpress,...
                   GasAmt,CKD_0,selfmult,formult,jj,profname);
      elseif (CKD_0 == 51)
        %%%scum = mst51(outwave,temperature,press,partpress,...
        scum = mst51_tobin(outwave,temperature,press,partpress,...
                   GasAmt,CKD_0,selfmult,formult,jj,profname);
      elseif (CKD_0 == 52)
        scum = mst52(outwave,temperature,press,partpress,...
                   GasAmt,CKD_0,selfmult,formult,jj,profname);
      elseif (CKD_0 == 53)
        scum = mst53(outwave,temperature,press,partpress,...
                   GasAmt,CKD_0,selfmult,formult,jj,profname);
      elseif (CKD_0 == 55)
        scum = mst55(outwave,temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
      elseif (CKD_0 == 56)
        scum = mst56(outwave,temperature,press,partpress,...
                   GasAmt,CKD_0,selfmult,formult,jj,profname);
      elseif (CKD_0 == 60)
        scum = mst60(outwave,temperature,press,partpress,...
                   GasAmt,CKD_0,selfmult,formult,jj,profname);
      elseif (CKD_0 == 5)
          scum = tobin5(outwave,temperature,press,partpress,...
                   GasAmt,CKD_0,selfmult,formult,jj,profname);
      elseif (ismember(CKD_0,[0 21 22 23 24  1 2 3 4 6]))
        ponk = [jj press(jj) partpress(jj) temperature(jj) CKD];
        str = 'using calconwater_loc in one gulp ... ii [p T pp] CKD = ';
        fprintf(1,' %s %3i %8.6e %8.6e %8.6f %3i\n',str,ponk);
        scum = calconwater_loc(1,length(outwave),outwave,dff,...
                               length(press),temperature,press,partpress,...
                               GasAmt,CKD,selfmult,formult,jj);
      elseif (ismember(CKD_0,[25 27 32]))
        ponk = [jj press(jj) partpress(jj) temperature(jj) CKD];
        str  = ...
           'using calconwater_loc_ckd2p5,2p7,3p2 in one gulp ... ii [p T pp] CKD = ';
        fprintf(1,'%s %3i %8.6e %8.6e %8.6f %3i\n',str,ponk)
	if (ismember(CKD_0,25))
          scum = calconwater_loc_ckd2p5(1,length(outwave),outwave,dff,...
                               length(press),temperature,press,partpress,...
                               GasAmt,CKD,selfmult,formult,jj);
	elseif (ismember(CKD_0,27))
	  error('huh? CKD 27?')
          scum = calconwater_loc_ckd2p7(1,length(outwave),outwave,dff,...
                               length(press),temperature,press,partpress,...
                               GasAmt,CKD,selfmult,formult,jj);
	elseif (ismember(CKD_0,32))
          scum = calconwater_loc_ckd3p2(1,length(outwave),outwave,dff,...
                               length(press),temperature,press,partpress,...
                               GasAmt,CKD,selfmult,formult,jj);
        end			   
      end

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


  else % have to break it into frequency intervals
    disp('using calconwater_loc in mini gulps...')
    for jj=MinLayer:Step:MaxLayer %LOOP 1..100 = bottom -> top
      nn=(jj-MinLayer)/Step+1;
      index0=1:(MaxLenX-10);
      number=floor(length(outwave)/(MaxLenX-10));

      for kk=1:number             %LOOP OVER FREQUENCY INTERVALS
        index=index0+(kk-1)*(MaxLenX-10);
        if (CKD_0 == 50)
          scum = mst50(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 51)
          %%%scum = mst51(outwave(index),temperature,press,partpress,...
          scum = mst51_tobin(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 52)
          scum = mst52(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 53)
          scum = mst53(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 55)
          scum = mst55(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 56)
          scum = mst56(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 60)
          scum = mst60(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 5)
          scum = tobin5(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (ismember(CKD_0,[0 21 22 23 24  1 2 3 4 6]))
          scum = calconwater_loc(1,length(index),outwave(index),dff,...
                             length(press),temperature,press,partpress,...
                             GasAmt,CKD,selfmult,formult,jj);
        elseif (ismember(CKD_0,[25]))
          scum =calconwater_loc_ckd2p5(1,length(index),outwave(index),dff,...
                             length(press),temperature,press,partpress,...
                             GasAmt,CKD,selfmult,formult,jj);
        elseif (ismember(CKD_0,[27]))
	  error('huh? CKD 27?')
          scum =calconwater_loc_ckd2p7(1,length(index),outwave(index),dff,...
                             length(press),temperature,press,partpress,...
                             GasAmt,CKD,selfmult,formult,jj);
        elseif (ismember(CKD_0,[32]))
          scum =calconwater_loc_ckd3p2(1,length(index),outwave(index),dff,...
                             length(press),temperature,press,partpress,...
                             GasAmt,CKD,selfmult,formult,jj);
        end
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
        if (CKD_0 == 50)
          scum = mst50(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 51)
          %%%scum = mst51(outwave(index),temperature,press,partpress,...
          scum = mst51_tobin(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 52)
          scum = mst52(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 53)
          scum = mst53(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 55)
          scum = mst55(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 56)
          scum = mst56(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 60)
          scum = mst60(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (CKD_0 == 5)
          scum = tobin5(outwave(index),temperature,press,partpress,...
                     GasAmt,CKD_0,selfmult,formult,jj,profname);
        elseif (ismember(CKD_0,[0 21 22 23 24  1 2 3 4 6]))
          scum = calconwater_loc(1,length(index),outwave(index),dff,...
                 length(press),...
                 temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
        elseif (ismember(CKD_0,[25]))
          scum =calconwater_loc_ckd2p5(1,length(index),outwave(index),dff,...
                   length(press),...
                 temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
        elseif (ismember(CKD_0,[27]))
	  error('huh? CKD 27?')
          scum =calconwater_loc_ckd2p7(1,length(index),outwave(index),dff,...
                   length(press),...
                 temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
        elseif (ismember(CKD_0,[32]))
          scum =calconwater_loc_ckd3p2(1,length(index),outwave(index),dff,...
                   length(press),...
                 temperature,press,partpress,GasAmt,CKD,selfmult,formult,jj);
        end
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
