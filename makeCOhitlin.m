function [PQRline] = makeCOhitlin(band,vers,strengthM,homepath,hitlin_fname,xfar);

% fcn [PQRline]=makeDAVEhitlin(band,vers,strengthM,homepath,hitlin_fname,xfar)
% makes a dave tobin style hit_c02 file ... does not overstamp the file if it
% previously exits with correct parameters

gasID=5;

%assume isotope ======= 1 unless stated
isotope = 1;

%assume we have to make a hit* file
do_makefile = +1;

%find all PQR lines for isotope 1 (or 2 or 3) : these are PQR_sigpi
if (band == 2150)
  v_l = 1; v_u = 2;
else
  fprintf('no need to generate hitBLAH.mat file!!');
  do_makefile = -1;
end

if (do_makefile > 0)
  %%%%% some of these CO bands go quite far!!!!!! eg P2350 is from 2160-2350
  %call hittomat to get the lines
  start = band-10*xfar;  
  if (start < 0) 
    start = 0.0;
  end
  stop = band+10*xfar;

  if ~exist('stretch')
    [start stop]
    stretch = 500;
  end
  if length(stretch) == 0
    [start stop]
    stretch = 500;
  end
  
  start = band-stretch;  
  stop = band+stretch;

  %now see if the .mat file exists with these parameters
  doexist = alreadyexist(gasID,homepath,band,start,stop,vers,strengthM,1);
  
  if (doexist  < 0)  %.mat file DNE .... create it
    strQ = ['save ' homepath 'CO_MATFILES/hit'];

    %line = hittomat2(start,stop,vers,strengthM,gasID);
    [line,hitran_version] = hitread(start,stop,strengthM,gasID,hitlin_fname);
    
    line = should_I_translate2oldHITparams(line,gasID,hitran_version);

    %get the branch info P Q or R
    pqr = line.bslq(:,5);

    ind = find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
      (line.iso == isotope) & ...
       ((pqr' == 'P') | (pqr' == 'Q') | (pqr' == 'R'))  );

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
    strQ = ['load ' homepath 'CO_MATFILES/hit'];
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
