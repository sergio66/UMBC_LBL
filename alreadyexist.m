function [doexist] = alreadyexist(gasID,homepath,band,start,stop,vers,strengthM,rw);
% function [doexist] = alreadyexist(gasID,homepath,band,start,stop,vers,strengthM,rw);
% this function sees if correct mat file exists for the current band

doexist = -1;    %assume it DNE

if gasID == 2
  filename=[homepath 'CO2_MATFILES/hit_info' num2str(band) '.mat'];
elseif gasID == 5
  filename=[homepath 'CO_MATFILES/hit_info' num2str(band) '.mat'];
end

if rw > 0                     %open the file for reading info
  check=1;              

  fid=fopen(filename,'r');

  if (fid == -1)
     check = -1;        %file DNE
  elseif (fid  ~= -1)
    ab = fread(fid,inf,'float64');
    fclose(fid);

    start_from_file     = ab(1);
    stop_from_file      = ab(2);
    vers_from_file      = ab(3);
    strengthM_from_file = ab(4);

    if (abs(start-start_from_file) > 0.1)
      check=-1;
    end

    if (abs(stop-stop_from_file) > 0.1)
      check=-1;
    end

    if (abs(vers-vers_from_file) > 0.1)
      check=-1;
    end

    ratio=0.0;
    if (strengthM_from_file >= 1e-100)
      ratio=strengthM/strengthM_from_file;
    elseif ((strengthM_from_file <= 1e-100) & (strengthM <= 1e-100))
      ratio=1.0;
    end  

    if ( (ratio < 0.99999) | (ratio > 1.00001))
      check=-1;
    end

  end
  doexist=check;
end

if rw <  0                     %open the file for writing info
  check=1;              
  fid=fopen(filename,'w');

  ab(1)=start;
  ab(2)=stop;
  ab(3)=vers*1.0;
  ab(4)=strengthM;

  fwrite(fid,ab,'real*8');
  fclose(fid);

  doexist=check;
end
  

