function x = read_MT_CKD4p3_sergio(fin)

%% /home/sergio/SPECTRA/CKDLINUX/MT_CKD-4.3/cntnm_sergio_v4.3_linux_gnu_dbl < n2_input

iCnt = 0;
iData = 0;
fid = fopen(fin,'r');

tline = fgets(fid);
iCnt = iCnt + 1; 
while ischar(tline)
  tline = fgets(fid);
  iCnt = iCnt + 1;
  
  if iCnt == 3
    x.P_T = tline;
    fprintf(1,'%s \n',tline);
  elseif iCnt == 8
    x.Q= tline;
    fprintf(1,'%s \n',tline);    
  end
  
  if iCnt >= 12
    if ischar(tline)
      iData = iData + 1;
      x.data(iData,1:2) = sscanf(tline,'%f',2);
    end
  end
end
fclose(fid);
