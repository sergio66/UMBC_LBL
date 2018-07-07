function [iYes,lines] = read_LBLRTM_LNFL(f1,f2,gid,iVers)

%% reads the LBLRTM LNFL data

lines = struct;

if nargin == 3
  iVers = 12.4;
  iVers = 12.8;  
end

if iVers == 12.2
  dir0 = ['/home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LNFL2.6/aer_v_3.2/line_files_By_Molecule/'];  
  if gid > 39
    error('can only do gases 1-39');
  end
elseif iVers == 12.4
  dir0 = ['/home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.4/LNFL/run_examples/run_example_infrared/TAPE7_aer_v_3.4_ex'];
  %% aslo see binary file                                      LBLRTM/lblrtm/run_examples/TAPE3_files/TAPE3_aer_v_3.4_ex_little_endian
  %% where the header shows you names of 47 gases (which correspond to those in /home/sergio/KCARTA/DOC/gasids_H2012
  
  dir0 = ['/home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.4/LNFL/aer_v_3.4_sergio/line_files_By_Molecule/'];
  if gid > 47
    error('can only do gases 1-47');
  end
elseif iVers == 12.8
  dir0 = ['/home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.8/LNFL/lnfl/run_examples/run_example_infrared_sergio/TAPE7_aer_v_3.6'];
  %% aslo see binary file                                      LBLRTM/lblrtm/run_examples/TAPE3_files/TAPE3_aer_v_3.6_ex_little_endian
  %% where the header shows you names of 47 gases (which correspond to those in /home/sergio/KCARTA/DOC/gasids_H2012
  
  dir0 = ['/home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.8/LNFL/aer_v_3.4_sergio/line_files_By_Molecule/'];
  dir0 = ['/home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.8/LINEDATAFILE/aer_v_3.6/line_files_By_Molecule/'];  
  if gid > 47
    error('can only do gases 1-47');
  end
else
  error('can only do LBLRTM 12.2 or 12.4 or 12.8')
end

iYes = -1; %% pretend no lines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
see make_breakapart.m in
  dir0 = ['/home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.4/LNFL/aer_v_3.4_sergio/line_files_By_Molecule/'];  
  [sergio@maya-usr1 SPECTRA]$ ls /home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.4/LNFL/aer_v_3.4_sergio/
  junk.lines  junk.txt  line_files_By_Molecule  make_breakapart.m
%}

fprintf(1,'looking for MOLGAS data in %s \n',dir0);

thedir = dir([dir0 '/*_*']);
iFound = -1;
gg = 1;
while iFound < 0 & gg <= 47
  zname = thedir(gg).name;
  zgid = str2num(zname(1:2));
  if zgid == gid
    iFound = +1;
  else
    gg = gg + 1;
  end
end

if iFound < 0
  disp('could not find LNFL files!')
  return
else
  fname = [dir0 '/' zname '/' zname];
  fprintf(1,'will be looking at %s \n',fname)
end

fnamex  = ['/home/sergio/xlnfl_junk_' num2str(gid,'%02d')];
fnamey  = ['/home/sergio/ylnfl_junk_' num2str(gid,'%02d')];

i3or4 = 3; %% old
i3or4 = 4; %% new

%% i3or4 = input('enter i3or4 : ');

if i3or4 == 3
  %%  ORIG
  sedder = ['!sed  -e ''/>/d'' -e ''/%/d''  ' fname ' > ' fnamex];
  eval(sedder);
  awker = ['!awk ''{print $1 " "  $2 " "  $3}'' ' fnamex ' > ' fnamey];
  eval(awker)
  %% 3 columns are [gid+iso   wnum stren]

  a = load(fnamey);
  woo = find(a(:,2) >= f1 & a(:,2) <= f2);
  if length(woo) > 0
    iYes = 1;
    lines.wnum = a(woo,2);
    lines.stren = a(woo,3);
    semilogy(a(woo,2),a(woo,3),'.')
    clf; scatter(a(woo,2),log10(a(woo,3)),10,ones(size(woo)),'filled'); colorbar; colormap jet    
  end

else

  %% NEW
  fnamex1  = ['/home/sergio/x1lnfl_junk_' num2str(gid,'%02d')];
  sedder = ['!sed  -e ''s/\(^ \{0,1\}[0-9]\{1,2\}\)\([0-9]\{1\} \)/\1 \2/p''  ' fname ' > ' fnamex1];
  eval(sedder);

  sedder = ['!sed  -e ''/>/d'' -e ''/%/d''  ' fnamex1 ' > ' fnamex];
  eval(sedder);
  awker = ['!awk ''{print $1 " "  $2 " "  $3 " "  $4}'' ' fnamex ' > ' fnamey];
  eval(awker)
  %% 4 columns are [gid  iso   wnum stren]

  rmer = ['!/bin/rm ' fnamex1];
  eval(rmer);
  
  a = load(fnamey);
  woo = find(a(:,3) >= f1 & a(:,3) <= f2);
  if length(woo) > 0
    iYes = 1;
    lines.iso = a(woo,2);    
    lines.wnum = a(woo,3);
    lines.stren = a(woo,4);
    clf; scatter(a(woo,3),log10(a(woo,4)),10,a(woo,2),'filled'); colorbar; colormap jet
  end

end

ax = axis;
axis([min(lines.wnum) max(lines.wnum) ax(3) ax(4)]); grid

rmer = ['!/bin/rm ' fnamex ' ' fnamey];
eval(rmer);

%{
disp('tic')
tic
fid = fopen(fname);
tline = fgetl(fid);
while ischar(tline)
   tline = fgetl(fid);
end
fclose(fid);
toc
%}
