function [hartmann_bands,lineOUT] = co2lines_jmhartmann(...
                          wn1,wn2,lineIN,which_isotope,exchangelinecenters);

%%% this loads in the data from JM Hartmann IR and sees which bands he uses
%%% for linemixing; 
%%% input : 
%%%   wn1,wn2 are the start/stop of the chunk you want to compute optdepths
%%%   lineIN   is from hitread
%%%   isotopes are the isotope list
%%% output
%%%   hartmann_bands = see if we can identify them with usual run8co2 IDs
%%%                    (if you find NaN, then it is not identified)
%%%   lineOUT  are the remaining HITRAN lineparams
%%% remove these from the lineIN structure and then use use regular run8.m 
%%% code for the rest of the lines

fhartdir = '/home/sergio/SPECTRA/JMHARTMANN_2007/LM_PQR_CO2_2.0/Data/';
fhart    = [fhartdir 'BandInfo.dat'];
hartmann = load(fhart);

[m,n] = size(hartmann);
usethis = [];
if which_isotope == 0
  %% all isotopes
  for ii = 1 : m
    if (hartmann(ii,6) <= wn2 & hartmann(ii,7) >= wn1)
      usethis = [usethis ii];
      end
    end
else
  for jj = 1 : length(which_isotope)
    jx = which_isotope(jj);
    for ii = 1 : m
      if (hartmann(ii,6) <= wn2 & hartmann(ii,7) >= wn1 & hartmann(ii,1) == jx)
        usethis = [usethis ii];
        end
      end
    end
  end

hartmann_bands = usethis;

ggall   = [];
bandall = [];

if exchangelinecenters < 0
  %% then we are removing the lines, called from run8co2_hartmann
  if length(usethis) > 0
    str='  #  Iso    US    LS  len(Hart)  len(Hitran) len(common)  Sum(SoFar)';
    disp(str);
    disp('-------------------------------------------------------------------')
    for ii = 1 : length(usethis)
      iii = usethis(ii);
      xx = hartmann(iii,:);

      %%% these are ALL the lines we can find in current HITRAN database
      [gg,lenHitr,lenHart] = find_hartmann_hitran(lineIN,fhartdir,xx);
      ggall = [ggall gg];
      figure(1); 
      semilogy(lineIN.wnum,lineIN.stren,'b.',...
               lineIN.wnum(ggall),lineIN.stren(ggall),'g.',...
               lineIN.wnum(gg),lineIN.stren(gg),'r.'); 
      title('looking for JM Hartmann bands (in red)')
      data = [ii xx(1) xx(2) xx(3) lenHart lenHitr length(gg) length(ggall)];
      fprintf(1,'%3i  %3i   %3i    %3i    %3i     %4i     %7i     %7i\n',data')
      pause(0.1);
      end
    end

  allin = 1 : length(lineIN.wnum);
  fprintf(1,'--> original number of lines = %7i \n',length(lineIN.wnum));
  fprintf(1,'--> number of lines in JMHartmanns bands = %7i \n',length(ggall));

  if length(usethis) > 0
    allout = setdiff(allin,ggall);
    lineOUT.linct  = length(allout);
    lineOUT.iso    = lineIN.iso(allout);  
    lineOUT.wnum   = lineIN.wnum(allout); 
    lineOUT.stren  = lineIN.stren(allout); 
    lineOUT.tprob  = lineIN.tprob(allout); 
    lineOUT.abroad = lineIN.abroad(allout); 
    lineOUT.sbroad = lineIN.sbroad(allout); 
    lineOUT.els    = lineIN.els(allout); 
    lineOUT.abcoef = lineIN.abcoef(allout); 
    lineOUT.tsp    = lineIN.tsp(allout); 
    lineOUT.iusgq  = lineIN.iusgq(allout); 
    lineOUT.ilsgq  = lineIN.ilsgq(allout);
    lineOUT.uslq   = lineIN.uslq(allout,:); 
    lineOUT.bslq   = lineIN.bslq(allout,:);
    lineOUT.ai     = lineIN.ai(allout,:); 
    lineOUT.ref    = lineIN.ref(allout,:);
    lineOUT.gasid  = lineIN.gasid(allout);  
    lineOUT.igas   = lineIN.igas;
  else
    lineOUT = lineIN;
    end

elseif exchangelinecenters > 0
  %% then we are exchanging linecenters, called from run8co2_linemixUMBC
  if length(usethis) > 0
    str='  #  Iso    US    LS  len(Hart)  len(Hitran) len(common)  Sum(SoFar)';
    disp(str);
    disp('-------------------------------------------------------------------')
    for ii = 1 : length(usethis)
      iii = usethis(ii);
      xx = hartmann(iii,:);

      %%% exchange ALL the lines we can find in current HITRAN database
      [lineOUT,gg,lenHitr,lenHart] = ...
          exchange_hartmann_hitran(lineIN,fhartdir,xx);
      lineIN = lineOUT;
      ggall = [ggall gg];
      figure(1); 
      semilogy(lineIN.wnum,lineIN.stren,'b.',...
               lineIN.wnum(ggall),lineIN.stren(ggall),'g.',...
               lineIN.wnum(gg),lineIN.stren(gg),'r.'); 
      title('looking for JM Hartmann bands (in red)')
      data = [ii xx(1) xx(2) xx(3) lenHart lenHitr length(gg) length(ggall)];
      fprintf(1,'%3i  %3i   %3i    %3i    %3i     %4i     %7i     %7i\n',data')
      pause(0.1);
      end
    end

  allin = 1 : length(lineIN.wnum);
  fprintf(1,'--> original number of lines = %7i \n',length(lineIN.wnum));
  fprintf(1,'--> number of lines in JMHartmanns bands = %7i \n',length(ggall));
  end
