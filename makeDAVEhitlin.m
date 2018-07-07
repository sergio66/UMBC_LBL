function [PQRline] = makeDAVEhitlin(band,vers,strengthM,homepath,...
                                  hitlin_fname,xfar,exchangelinecenters)

% fcn [PQRline]=makeDAVEhitlin(band,vers,strengthM,homepath,hitlin_fname,xfar)
% makes a dave tobin style hit_c02 file ... does not overstamp the file if it
% previously exits with correct parameters

global quiet stretch p2311_21_jmax

if ~exist('p2311_21_jmax')
  p2311_21_jmax = 600;
end
if length(p2311_21_jmax) == 0
  p2311_21_jmax = 600;
end

gasID=2;

%assume isotope ======= 1 unless stated
isotope = 1;

%assume we have to make a hit* file
do_makefile = +1;

%find all PQR lines for isotope 1 (or 2 or 3) : these are PQR_sigpi
if (band == 618)
  v_l = 2; v_u = 3;
elseif (band == 648)
  v_l = 1; v_u = 2;        isotope = 2;
elseif (band == 662)
  v_l = 1; v_u = 2;        isotope = 3;
elseif (band == 667)
  v_l = 1;  v_u = 2; 
elseif (band == 720)
  v_l = 2;  v_u = 5; 
elseif (band == 791)
  v_l = 3; v_u = 8;
elseif (band == 1932)
  v_l = 1; v_u = 6;
elseif (band == 2080)
  v_l = 1; v_u = 8;
elseif (band == 2129)
  v_l = 2; v_u = 15;

%find all PQR lines for isotope 1 : these are PQR_deltpi
elseif (band == 668)
  v_l = 2; v_u = 4;
elseif (band == 740)
  v_l = 4; v_u = 8;
elseif (band == 2093)
  v_l = 2; v_u = 14;

%find all PQR lines for isotope 1 : these are PQR_sigsig
elseif (band == 2350)
  v_l = 1; v_u = 9;
elseif (band == 2351)
  v_l = 1; v_u = 9;        isotope = 2;
elseif (band == 2352)
  v_l = 1; v_u = 9;        isotope = 3;
elseif (band == 2353)
  v_l = 3; v_u = 23; 
elseif (band == 2354)
  v_l = 5; v_u = 25; 
%elseif (band == 2355)  isotope = 4;
%  v_l = 1; v_u = 9; 

%find all PQR lines for isotope 1 : these are PQR_pipi
elseif (band == 2320) 
  v_l = 2; v_u = 16; 
elseif (band == 2321) 
  v_l = 2; v_u = 16;       isotope = 2;
elseif (band == 2322) 
  v_l = 2; v_u = 16;       isotope = 3;

%find all PQR lines for isotope 1 : these are PQR_deltdelt
elseif (band == 2310) 
  v_l = 4; v_u = 24;
elseif (band == 2311) 
  v_l = 4; v_u = 24;       isotope = 2;

else
  fprintf('no need to generate hitBLAH.mat file!!');
  do_makefile = -1;
end

if (do_makefile > 0)
  %%%%% some of these CO2 bands go quite far!!!!!! eg P2350 is from 2160-2350
  %call hittomat to get the lines
  start = band-10*xfar;  
  if (start < 0) 
    start = 0.0;
  end
  stop = band+10*xfar;

  if ~exist('stretch')
    start
    stop
    stretch = 500;
  end
  if length(stretch) == 0
    start
    stop
    stretch = 500;
  end
  
  start = band-stretch;  
  stop = band+stretch;

  %now see if the .mat file exists with these parameters
  doexist = alreadyexist(gasID,homepath,band,start,stop,vers,strengthM,1);
  
  if (doexist  < 0)  %.mat file DNE .... create it
    strQ = ['save ' homepath 'CO2_MATFILES/hit'];

    %line = hittomat2(start,stop,vers,strengthM,gasID);
    [line,hitran_version] = hitread(start,stop,strengthM,gasID,hitlin_fname);
    
    line = should_I_translate2oldHITparams(line,gasID,hitran_version);

    %get the branch info P Q or R
    pqr = line.bslq(:,5);

    pband_iso = [2311 2321 2351];
    inside    = length(intersect(band,pband_iso));
    if (inside <= 0)
      %%%%%no big deal, just find all P Q R lines
      ind = find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
        (line.iso == isotope) & ...
         ((pqr' == 'P') | (pqr' == 'Q') | (pqr' == 'R'))  );
    elseif ((inside > 0) & (p2311_21_jmax < 0))
      %%%%%no big deal, just find all P Q R lines as above
      ind = find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
        (line.iso == isotope) & ...
        ((pqr' == 'P') | (pqr' == 'Q') | (pqr' == 'R'))  );
    elseif ((inside > 0) & (p2311_21_jmax > 0))
      %%%%%%only want the lower P values, all R values
      find_lowest = p2311_21_jmax;
      indR = find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
          (line.iso == isotope) & ((pqr' == 'Q') | (pqr' == 'R'))  );
      indP = find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
        (line.iso == isotope) & (pqr' == 'P')  );
      jval = str2num(line.bslq(indP,6:8));

      [yyy,iii] = sort(jval);           %%%%%resort from lowest to highest  
      jval = jval(iii);
      indP = indP(iii);

      index = find(jval <= find_lowest);
      indP = indP(index);

      ind = [indP indR];               %%combine P and R
    end

    ii = length(ind);
    fprintf(1,'number of PQR lines selected for band %5i = %5i \n',band,ii);

    accuracy    = line.ai(ind,1:3);
    dipole      = line.tprob(ind)';    %transition probablility  
    elower      = line.els(ind)';
    freq        = line.wnum(ind)';
    gas_id      = gasID;
    iso         = line.iso(ind)';
    j_lower     = line.bslq(ind,1:9);
    j_upper     = line.uslq(ind,1:9);
    line_status = ones(length(ind),1)*vers;
    p_shift     = line.tsp(ind)';
    reference   = line.ref(ind,1:6);
    stren       = line.stren(ind)';
    v_lower     = line.ilsgq(ind)';
    v_upper     = line.iusgq(ind)';
    w           = line.abroad(ind)';  %air broadenend widths
    w_s         = line.sbroad(ind)';  %self broadened widths
    w_temp      = line.abcoef(ind)';  %tempr correct for air broadened widths

    if exchangelinecenters > 0
      disp('exchanging the HITRAN linecenters with HARTMANN');
      freqold = freq;
      [freq,str] = exchange_f1_f2(line,ind,iso(1),v_upper(1),v_lower(1));
      plot(freq,freq-freqold); title(str); pause(0.1);
    end     
  
    str1 = ' accuracy dipole elower freq gas_id iso j_lower j_upper ';
    str2 = 'line_status p_shift reference stren v_lower v_upper w w_s w_temp';

    eval([strQ num2str(band) '.mat' str1 str2])

    %save hit720.mat accuracy dipole elower freq gas_id iso j_lower j_upper ...
    %     line_status p_shift reference stren v_lower v_upper w w_s w_temp

    PQRline.ai    = line.ai(ind,1:3);
    PQRline.tprob = line.tprob(ind)';    %transition probablility
    PQRline.els   = line.els(ind)';
    PQRline.wnum  = line.wnum(ind)';
    PQRline.gasID  = gasID;
    PQRline.iso   = line.iso(ind)';
    PQRline.bslq  = line.bslq(ind,1:9);
    PQRline.uslq  = line.uslq(ind,1:9);
    PQRline.tsp   = line.tsp(ind)';
    PQRline.ref   = line.ref(ind,1:6);
    PQRline.stren = line.stren(ind)';
    PQRline.ilsgq = line.ilsgq(ind)';
    PQRline.iusgq = line.iusgq(ind)';
    PQRline.abroad= line.abroad(ind)';  %air broadenend widths
    PQRline.sbroad= line.sbroad(ind)';  %self broadened widths
    PQRline.abcoef= line.abcoef(ind)';  %tempr correct for air broadened widths

    semilogy(freq,stren);

    doexist = alreadyexist(gasID,homepath,band,start,stop,vers,strengthM,-1);

  end

  if (doexist  > 0)      %.mat file EXISTS .... read it
    strQ = ['load ' homepath 'CO2_MATFILES/hit'];
    eval([strQ num2str(band) '.mat'])

    PQRline.ai     = accuracy;
    PQRline.tprob  = dipole;    %transition probablility
    PQRline.els    = elower;
    PQRline.wnum   = freq;
    PQRline.gasID  = gas_id;
    PQRline.iso    = iso;
    PQRline.bslq   = j_lower;
    PQRline.uslq   = j_upper;
    PQRline.tsp    = p_shift;
    PQRline.ref    = reference;
    PQRline.stren  = stren;
    PQRline.ilsgq  = v_lower;
    PQRline.iusgq  = v_upper;
    PQRline.abroad = w;  %air broadenend widths
    PQRline.sbroad = w_s;  %self broadened widths
    PQRline.abcoef = w_temp;  %tempr correct for air broadened widths

    semilogy(freq,stren); pause(0.1);
  end

elseif (do_makefile < 0)
  %%call hittomat to get the lines;   %this is really useless
  disp('nothing to do in makeDAVEhitlin')
  iDO = -1; 
  if iDO > 0
    start = band-stretch;
    stop = band+stretch;
    %%PQRline = hittomat2(start,stop,vers,strengthM,gasID);
    [PQRline,hitran_version]=hitread(start,stop,strengthM,gasID,hitlin_fname);
    PQRline = should_I_translate2oldHITparams(PQRline,gasID,hitran_version);
  end
end
