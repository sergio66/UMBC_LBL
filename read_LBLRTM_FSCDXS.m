function [iYes,lines] = read_LBRTM_FSCDXS(f1,f2,gid,iVers)

%% reads the LBLRTM LNFL data

if nargin == 3
  iVers = 12.4;
  iVers = 12.8;  
end

if iVers == 12.2
  fname0 = ['/home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LNFL2.6/aer_v_3.2/line_files_By_Molecule/'];
  error('not worried about this vers')
  if gid < 51 | gid > 63
    error('can only do gases 51-63');
  end
elseif iVers == 12.4
  if gid < 51 | gid > 83
    error('can only do gases 51-83');
  end
  fname0 = ['/home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.4/LBLRTM/lblrtm/run_examples/xs_files_v3.4/FSCDXS'];
elseif iVers == 12.8
  if gid < 51 | gid > 83
    error('can only do gases 51-83');
  end
  fname0 = ['/home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.8/LBLRTM/lblrtm/run_examples/xs_files/FSCDXS'];
else
  error('can only do LBLRTM 12.2 or 12.4 or 12.8')
end

iYes = -1; %% pretend no lines

iFound  = dir(fname0);
if length(iFound) == 0
  disp('could not find LBLRTM XCS files!')
  return
else
  fprintf(1,'looking for XSEC data in %s \n',fname0);
end

fid = fopen(fname0);
tline = fgetl(fid);
tline = fgetl(fid);
%fprintf(1,'%s\n',tline);
iCnt = 0;
while ischar(tline)
  tline = fgetl(fid);
  if ~strcmp(tline(1),' ') & ischar(tline)
    %fprintf(1,'%s\n',tline);
    name = deblank(tline(1:10));
   
    if iCnt == 0
      iCnt = iCnt + 1;
      lblrtm_xsc{iCnt} = name;
      iaFound(iCnt) = 1;
      bandStart(iCnt,1) = str2num(tline(11:20));
      bandStop(iCnt,1)  = str2num(tline(21:30));      
    else
      iFound = -1;
      ii = 1;
      while ii <= iCnt & iFound < 0
        oldname = char(lblrtm_xsc{ii});
        if strcmp(oldname,name)
          iFound = 1;
          iaFound(ii) = iaFound(ii) +1;
          bandStart(ii,iaFound(ii)) = str2num(tline(11:20));
          bandStop(ii,iaFound(ii))  = str2num(tline(21:30));      
        else
          ii = ii + 1;
        end
      end
      if iFound < 0
        iCnt = iCnt + 1;
        lblrtm_xsc{iCnt} = name;
        iaFound(iCnt) = 1;
        bandStart(iCnt,1) = str2num(tline(11:20));
        bandStop(iCnt,1)  = str2num(tline(21:30));      	
      end
    end   %% if iCnt
  end     %% if ~strcmp
end       %% while
fclose(fid);

%{
%% this is the summary of what was just read in
for ii = 1 : iCnt
  fprintf(1,'%2i) %s has %2i bands : \n',ii,char(lblrtm_xsc{ii}),iaFound(ii))
  for jj = 1 : iaFound(ii)
    fprintf(1,'    %10.6f   %10.6f \n',bandStart(ii,jj),bandStop(ii,jj))
  end
end
%}

umbc_xscname = gid_to_gname(gid);
umbc_xscname = strtrim(umbc_xscname(2:end));
iYes = -1;
lines = struct;
for ii = 1 : iCnt
  lblrtm_name = strtrim(char(lblrtm_xsc{ii}));
  if strcmp(umbc_xscname,lblrtm_name)
    for jj = 1 : iaFound(ii)
      %fprintf(1,'    %10.6f   %10.6f \n',bandStart(ii,jj),bandStop(ii,jj))
      bS = bandStart(ii,jj);
      bE = bandStop(ii,jj);
      if f1 < bE & f2 > bS
        iYes = 1;
	lines.bStart = bS;
	lines.bSop   = bE;
	fprintf(1,'%2i --> %s  %10.6f   %10.6f \n',gid,lblrtm_name,bandStart(ii,jj),bandStop(ii,jj))
      end	
    end
  end
end

