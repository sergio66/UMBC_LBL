function [hh,ll] = band_end(gasID,thebands,vers,strengthM,fname,hitran_version);
%this function finds the start/stop points of the bands, so that they can
%be included as bands in the computation

global stretch  p2311_21_jmax

hh = -100.0e5*ones(size(thebands));
ll = +100.0e5*ones(size(thebands));   
%set up dummy numbers as initial values

%gasID = 2;

if gasID == 2

  for iiii = 1:length(thebands)
    band = thebands(iiii);
  
    %assume isotope ======= 1 unless stated
    isotope = 1;

    %assume we have to make a hit* file
    do_makefile = +1;

    %from makeDAVEhitlin.m
    %find all PQR lines for isotope 1 (or 2 or 3) : these are PQR_sigpi
    if (band == 618)
      v_l = 2; v_u = 3;
    elseif (band == 648)
      v_l = 1; v_u = 2; isotope = 2;
    elseif (band == 662)
      v_l = 1; v_u = 2; isotope = 3;
    elseif (band == 667)
      v_l = 1; v_u = 2; 
    elseif (band == 720)
      v_l = 2; v_u = 5; 
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
      v_l = 1; v_u = 9; isotope = 2;
    elseif (band == 2352)
      v_l = 1; v_u = 9; isotope = 3;
    elseif (band == 2353)
      v_l = 3; v_u = 23; 
    elseif (band == 2354)
      v_l = 5; v_u = 25; 
    %elseif (band == 2355)    isotope = 4;
    %    v_l = 1; v_u = 9; 

    %find all PQR lines for isotope 1 : these are PQR_pipi
    elseif (band == 2320) 
      v_l = 2; v_u = 16; 
    elseif (band == 2321) 
      v_l = 2; v_u = 16; isotope = 2;
    elseif (band == 2322) 
      v_l = 2; v_u = 16; isotope = 3;

    %find all PQR lines for isotope 1 : these are PQR_deltdelt
    elseif (band == 2310) 
      v_l = 4; v_u = 24;
    elseif (band == 2311) 
      v_l = 4; v_u = 24; isotope = 2;

    else
      do_makefile =    -1;
      error('check that band_end bands same as those in makeDAVE');
    end

    %call hittomat to get the lines
    start = band-stretch;
    stop = band+stretch;

    %line    =    hittomat2(start,stop,vers,strengthM,gasID);
    line    =    hitread(start,stop,strengthM,gasID,fname);
    line = should_I_translate2oldHITparams(line,2,hitran_version);

    %get the branch info P Q or R
    pqr = line.bslq(:,5);

    pband_iso = [2311 2321 2351];
    inside = length(intersect(band,pband_iso));
    if (inside <= 0)
      %%%%%no big deal, just find all P Q R lines
      ind = find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
                  (line.iso == isotope) & ...
                  ((pqr' == 'P') | (pqr' == 'Q') | (pqr' == 'R'))    );
      ii = length(ind);
    elseif ((inside > 0) & (p2311_21_jmax < 0))
      %%%%%no big deal, just find all P Q R lines as above
      ind = find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
                 (line.iso == isotope) & ...
                ((pqr' == 'P') | (pqr' == 'Q') | (pqr' == 'R'))    );
      ii = length(ind);
    elseif ((inside > 0) & (p2311_21_jmax > 0))
      %%%%%%only want the lower P values, all R values
      find_lowest = p2311_21_jmax;
      indR = find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
                (line.iso == isotope) & ((pqr' == 'Q') | (pqr' == 'R'))    );
      indP = find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
                (line.iso == isotope) & (pqr' == 'P')    );
      jval = str2num(line.bslq(indP,6:8));

      [yyy,iii] = sort(jval);                     %%%%%resort from lowest to highest    
      jval = jval(iii);
      indP = indP(iii);

      index = find(jval <= find_lowest);
      indP = indP(index);

      ind = [indP indR];                             %%combine P and R
    end

    freq         = line.wnum(ind)';
    hh(iiii) = max(freq);
    ll(iiii) = min(freq);
  end
elseif gasID == 5
  for iiii = 1:length(thebands)
    band = thebands(iiii);
  
    %assume isotope ======= 1 unless stated
    isotope = 1;

    %assume we have to make a hit* file
    do_makefile = +1;

    %from makeDAVEhitlin.m
    %find all PQR lines for isotope 1 (or 2 or 3) : these are PQR_sigpi
    if (band == 2150)
      v_l = 1; v_u = 2;
    else
      do_makefile =    -1;
      error('check that band_end bands same as those in makeDAVE');
    end

    %call hittomat to get the lines
    start = band-stretch;
    stop = band+stretch;

    %line    =    hittomat2(start,stop,vers,strengthM,gasID);
    line    =    hitread(start,stop,strengthM,gasID,fname);
    line = should_I_translate2oldHITparams(line,2,hitran_version);

    %get the branch info P Q or R
    pqr = line.bslq(:,5);

    %%%%%no big deal, just find all P Q R lines
    ind = find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
                  (line.iso == isotope) & ...
                  ((pqr' == 'P') | (pqr' == 'Q') | (pqr' == 'R'))    );
    ii = length(ind);

    freq         = line.wnum(ind)';
    hh(iiii) = max(freq);
    ll(iiii) = min(freq);
  end

end

