%% first find the LONG names of the dis of the XSEC gases (51-63)
disp('final files will be in IR-XSect/Uncompressed-files/')
disp('final files will be in IR-XSect/Uncompressed-files/')
disp('final files will be in IR-XSect/Uncompressed-files/')
disp(' ')


disp('there are 1 + 31 places to change to 2024 : lser, and ALL the *_IR2024.xsc')
disp('there are 1 + 31 places to change to 2024 : lser, and ALL the *_IR2024.xsc')
disp('there are 1 + 31 places to change to 2024 : lser, and ALL the *_IR2024.xsc')


%lser = ['!ls /asl/data/hitran/H2016/XSEC/*_g[56]*// | grep -i asl >& ugh'];
%lser = ['!ls /asl/data/hitran/H2016/XSEC/*_g[5678]*// | grep -i asl >& ugh'];
%% asl/data/hitran/H2016/ --> /umbc/xfs3/strow/asl/rta/hitran/H2024/
lser = ['!ls /umbc/xfs3/strow/asl/rta/hitran/H2024/XSEC/*_g[56]*// | grep -i asl >& ugh'];
lser = ['!ls /umbc/xfs3/strow/asl/rta/hitran/H2024/XSEC/*_g[5678]*// | grep -i asl >& ugh'];
eval(lser)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iStart = 01;  iStop = 13;   %% gids 51-63
iStart = 01;  iStop = 31;   %% gids 51-81

%% have to redo iCnt = 10 (gid 72)
%% iStart = 11;  iStop = 31;   %% gids 51-81

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iCnt = 0;
fid = fopen('ugh');
tline = fgetl(fid);
while ischar(tline)
  iCnt = iCnt + 1;
  fprintf(1,'%d:%s\n',iCnt,tline);
  str(iCnt).tline = tline(1:end-2);
  tline = fgetl(fid);
end
fclose(fid);
rmer = ['!/bin/rm ugh']; eval(rmer);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iSkip = 0;
%% then find the output file names, and also read in the multiple files, concat them together
for iCnt = iStart : iStop
  tline = str(iCnt).tline;
  gid = str2num(tline(end-1:end));
  fprintf(1,'iCnt = %2i  tline = %s gid = %2i \n',iCnt,tline,gid);
  
  if gid == 51
    outname = 'CFC-11_IR2024.xsc';
  elseif gid == 52
    outname = 'CFC-12_IR2024.xsc';
  elseif gid == 53
    outname = 'CFC-13_IR2024.xsc';
  elseif gid == 54
    outname = 'CFC-14_IR2024.xsc';
  elseif gid == 55
    outname = 'HCFC-21_IR2024.xsc';
  elseif gid == 56
    outname = 'HCFC-22_IR2024.xsc';
  elseif gid == 57
    outname = 'CFC-113_IR2024.xsc';
  elseif gid == 58
    outname = 'CFC-114_IR2024.xsc';
  elseif gid == 59
    outname = 'CFC-115_IR2024.xsc';
  elseif gid == 60
    outname = 'CCl4_IR2024.xsc';
  elseif gid == 61
    outname = 'ClONO2_IR2024.xsc';
  elseif gid == 62
    outname = 'N2O5_IR2024.xsc';
  elseif gid == 63
    outname = 'HNO4_IR2024.xsc';
  elseif gid == 64
    outname = 'C2F6_IR2024.xsc';
  elseif gid == 65
    outname = 'HCFC-123_IR2024.xsc';
  elseif gid == 66
    outname = 'HCFC-124_IR2024.xsc';
  elseif gid == 67
    outname = 'HCFC-141b_IR2024.xsc';
  elseif gid == 68
    outname = 'HCFC-142b_IR2024.xsc';  
  elseif gid == 69
    outname = 'HCFC-225ca_IR2024.xsc';  
  elseif gid == 70
    outname = 'HCFC-225cb_IR2024.xsc';  
  elseif gid == 71
    outname = 'HFC-32_IR2024.xsc';
  elseif gid == 72
    outname = 'HFC-134a_IR2024.xsc';
  elseif gid == 73
    outname = 'HFC-143a_IR2024.xsc';
  elseif gid == 74
    outname = 'HFC-152a_IR2024.xsc';
  elseif gid == 75
    outname = 'C6H6_IR2024.xsc';
  elseif gid == 76
    outname = 'HFC-125_IR2024.xsc';
  elseif gid == 77
    outname = 'HFC-134_IR2024.xsc';
  elseif gid == 78
    outname = 'SF5CF3_IR2024.xsc';    
  elseif gid == 79
    outname = 'PAN_IR2024.xsc';
  elseif gid == 80
    outname = 'CH3CN_IR2024.xsc';    
  elseif gid == 81
    outname = 'SF6_IR2024.xsc';    
  end

  iBahSkip = -1;
  thedir = dir([tline '/*.xsc']);
  thedirTorr = dir([tline '/*Torr*.xsc']);
  fprintf(1,'%2i %2i %s has %3i files for output of which %3i have Torr in their name %s \n',iCnt,gid,tline,length(thedir),length(thedirTorr),outname);
  if length(thedir) == length(thedirTorr)
    disp('oh no : all files have Torr in the name!!! Skipping!!');
    iSkip = iSkip + 1;
    iaSkip(iSkip) = iCnt;
    iaGasIDSkip(iSkip) = gid;
    iBahSkip = +1;
  end

  if iBahSkip < 0
    clear okname yes
    okname = ones(1,length(thedir));
    if length(thedirTorr) > 0
      for ii=1:length(thedir)
        fname = thedir(ii).name;
        if length(strfind(fname,'Torr')) == 0
          okname(ii) = 1;
        else
          okname(ii) = -1;      
        end
      end
    end
  
    yes = find(okname == 1);
    
    %% typical name : CCl3F_232.7_336.0_810.0-880.0_00.xsc
    clear T P startf stopf array barray index
    for ii = 1 : length(yes)
      fname = [thedir(yes(ii)).name];
      boo1 = findstr(fname,'_');
      boo2 = findstr(fname,'-');    
      T(ii) = str2num(fname(boo1(1)+1:boo1(2)-1));
      P(ii) = str2num(fname(boo1(2)+1:boo1(3)-1));
      startf(ii) = str2num(fname(boo1(3)+1:boo2(1)-1));
      stopf(ii) = str2num(fname(boo2(1)+1:boo1(4)-1));    
    end
    array = [startf; stopf; T; P]';
    %[barray,index] = sortrows(array,[1 3 4],{'ascend','descend','descend'});
    [barray,index] = sortrows(array,[+1 -3]);   
    [barray,index] = sortrows(array,[+1 -3 -4]);
    %[barray,index] = sortrows(array,[+1 -4 -3]);     
  
    if length(intersect(index,1:length(yes)))-length(yes) ~= 0    
      error('oh oh while sorting, something went wrong')
    end
  
    outnameX = ['IR-XSect/Uncompressed-files/' outname];
    fname = [tline '/' thedir(yes(index(1))).name];
    if gid == 79
      fname = strrep(fname,'(','\(');
      fname = strrep(fname,')','\)');    
    end    
    catter = ['!cat ' fname ' > ' outnameX]; eval(catter);
    for ii = 2 : length(yes)
      fname = [tline '/' thedir(yes(index(ii))).name];
      if gid == 79
        fname = strrep(fname,'(','\(');
        fname = strrep(fname,')','\)');    
      end        
      catter = ['!cat ' fname ' >> ' outnameX]; eval(catter);  
    end
    fprintf(1,'outnameX = %s \n',outnameX);
  end   %% iBahSkip < 0
  %disp('ret'); pause
  disp(' ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iSkip > 0
  disp('had to skip the following since all the filenames had Torr (iCnt & iGasID)')
  [iaSkip; iaGasIDSkip]'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to continue making plots ...')

plot_convert_xsec_multiplefiles2onefile
