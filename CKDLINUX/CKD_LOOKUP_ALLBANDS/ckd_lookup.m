%script file to produce fortran look up table at 1.0 cm-1

%%%%%%%%%%%%%%do the stuff at 1 cm-1 output spacing; no need to avg
%output spacing = ffin*nbox == 1.0 cm-1
ffin=0.0025; nbox=1/ffin;

%%%note we go from 100 to 3001 at 1cm-1 spacing; remember code does not do last freq point, so it will stop
%%%at 3000 cm-1

ckd_version

% CKD=24;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do self
[fr,ks]=run6watercontinuum(1,100,3001,ffin,0.1,0.5,1,1,2,25,nbox,0,0,1,CKD,1,0,1,1,'../SPECTRA/watercont');

fid=fopen('CKDself24.dat','w');

fprintf(fid,'c code made from /salsify/users/sergio/KCARTA/SPECTRA/CKD/ckd_lookup.m \n');
fprintf(fid,'      data CKDvers/ %3i / \n',CKD);
fprintf(fid,'      \n');
fprintf(fid,'      data rfr1,rfr2,npts / %10.4f,  %10.4f, %5i/ \n',fr(1),fr(length(fr)),length(fr));

numblocks=(length(fr)-1)/5;          %5 points per line in file
i=1;  
for ii=150:10:370
  ind=1:5;
  fprintf(fid,'c this is for T = %4i \n',ii);
  fprintf(fid,'      data ((ckdself(i,j),j=1,2901),i=%4i,%4i)/ \n',i,i);
  for jj=1:numblocks
    data=ks(i,ind);
    fprintf(fid,'     + %12.10e, %12.10e, %12.10e, %12.10e, %12.10e, \n',data);
    ind=5+ind;
    end
  data=ks(i,length(fr));
  fprintf(fid,'     + %12.10e / \n',data);
  i=i+1;  
  end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do for
[fr,kf]=run6watercontinuum(1,100,3001,ffin,0.1,0.5,1,1,2,25,nbox,0,0,1,CKD,0,1,1,1,'../SPECTRA/watercont');

fid=fopen('CKDfor24.dat','w');

fprintf(fid,'c code made from /salsify/users/sergio/KCARTA/SPECTRA/CKD/ckd_lookup.m \n');
fprintf(fid,'      data CKDvers/ %3i / \n',CKD);
fprintf(fid,'      \n');
fprintf(fid,'      data rfr1,rfr2,npts / %10.4f,  %10.4f, %5i/ \n',fr(1),fr(length(fr)),length(fr));

numblocks=(length(fr)-1)/5;          %5 points per line in file
i=1;  
for ii=150:10:370
  ind=1:5;
  fprintf(fid,'c this is for T = %4i \n',ii);
  fprintf(fid,'      data ((ckdfor(i,j),j=1,2901),i=%4i,%4i)/ \n',i,i);
  for jj=1:numblocks
    data=kf(i,ind);
    fprintf(fid,'     + %12.10e, %12.10e, %12.10e, %12.10e, %12.10e, \n',data);
    ind=5+ind;
    end
  data=kf(i,length(fr));
  fprintf(fid,'     + %12.10e / \n',data);
  i=i+1;  
  end

fclose(fid);
