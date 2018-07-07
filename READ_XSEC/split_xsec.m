function split_xsec(gf)

% split xsec file gf into component gasses
%

[fin, msg] = fopen(gf, 'r');
if fin == -1
  error(msg);
end

% read first header
rhdr = fgetl(fin); 

fout = -1;
rgstr = '';
rgstr_prev = '';

while isstr(rhdr) & length(rhdr) >= 10

  rgstr_prev = rgstr;
  rgstr = rhdr(1:10);

  if ~strcmp(rgstr, rgstr_prev)
     if fout ~= -1
       fclose(fout);
     end
     [fout, msg] = fopen([strtok(rgstr),'.xsc'], 'w');
     if fout == -1
        error(msg);
     end
  end

  % get header fields
  rgstr = rhdr(1:10);
  rv1   = str2num(rhdr(11:20));
  rv2   = str2num(rhdr(21:30));
  rnpts = str2num(rhdr(31:40));
  rtemp = str2num(rhdr(41:50));
  rpres = str2num(rhdr(51:60));
  rmaxi = str2num(rhdr(61:70));
  rjunk = rhdr(71:length(rhdr));

  % sanity check
  if rv1 < 10 | 50000 < rv2 | rtemp < 100 | 400 < rtemp 
    error('bad header record');
  end
    
  % write header
  fprintf(fout, '%s\n', rhdr);

  % write abs data
  nabsrow = floor((rnpts-1)/10) + 1;
  for j = 1:nabsrow
    absrow = fgetl(fin);
    fprintf(fout, '%s\n', absrow);
  end

  % read next header
  rhdr = fgetl(fin); 

end % end of xsec record read loop

fclose(fin);
fclose(fout);

