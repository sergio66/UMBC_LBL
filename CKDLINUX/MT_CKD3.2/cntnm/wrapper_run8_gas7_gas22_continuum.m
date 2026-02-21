function [fr5,od5] = wrapper_run8_gas7_gas22_continuum(gid,fr1,fr2,df,profile)

%% see /home/sergio/SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/wrapper_run8_gas7_gas22_continuum.m

%{
profile = '/home/sergio/SPECTRA/IPFILES/std_gx7x_6';
fr1 = 1555; fr2 = fr1 + 25;
%fr1 = fr1 - 2*0.0005; fr2 = fr2 - 0.0025 + 2*0.0005 + 0.0005;
[fr,od] = wrapper_run8_gas7_gas22_continuum(7,fr1,fr2,0.0005,profile);
oldrun8 = load('/asl/s1/sergio/H2016_RUN8_NIRDATABASE/IR_605_2830/g7.dat/std1555_7_6.mat');
plot(fr,sum(od),'b.',oldrun8.w,sum(oldrun8.d),'r')

profile = '/home/sergio/SPECTRA/IPFILES/std_gx22x_6';
fr1 = 2355;
fr1 = 2380; fr2 = fr1 + 25;
%fr1 = fr1 - 2*0.0005; fr2 = fr2 - 0.0025 + 2*0.0005 + 0.0005;
[fr,od] = wrapper_run8_gas7_gas22_continuum(22,fr1,fr2,0.0005,profile);
oldrun8 = load('/asl/s1/sergio/H2016_RUN8_NIRDATABASE/IR_605_2830/g22.dat/std2380_22_6.mat');
plot(fr,sum(od),'b.',oldrun8.w,sum(oldrun8.d),'r')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
o2 = [1338 1850];
n2 = [1996 2902];

%% this is going to assume a 25 cm-1 chunk

if (fr2-fr1) > 26
  error('only doing 25 cm-1 chunks')
end

fr1 = fr1 - 2*0.0005; fr2 = fr2 - 0.0025 + 2*0.0005 + 0.0005;
dalen0 = (fr2-fr1)*400;
if gid ~= 7 & gid ~= 22
  error('only doing gid 7,22')
end

prof = load(profile,'-ascii');
numlay = length(prof(:,4));
for ii = 1 : length(prof(:,4))
   fprintf(1,'doing layer %3i of %3i \n',ii,numlay);
  
   irand = floor(rand(1)*1e8);   
   while irand < 1e7;   
     irand = floor(rand(1)*1e8);
   end

   fin   = ['input.txt_' num2str(gid) '_' num2str(round(fr1))     '_' num2str(irand)];
   fout  = ['CNTNM.OPTDPT_' num2str(gid) '_' num2str(round(fr1))  '_' num2str(irand)];
   xfout = ['xCNTNM.OPTDPT_' num2str(gid) '_' num2str(round(fr1)) '_' num2str(irand)];

   fid = fopen(fin,'w');

   junk = prof(ii,:);
   junk(1) = gid;   
   fprintf(fid,'%3i %8.6e %8.6e %8.6f %8.6e \n',junk);
   
   junk = [fr1 fr2 0.0005];
   fprintf(fid,'%8.6e %8.6e %8.6e \n',junk);
   fprintf(fid,'%7i \n',irand);   
   fclose(fid);

   %% SEE ~/SPECTRA/Readme_make_new_qtipsH24
   %% see CKDLINUX/calcon_3p2.F:5:c see /home/sergio/SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/build/make_cntnm_sergio_run8
   %% first line says
   %%   Also try this                     gmake -f make_cntnm linuxINTELdbl
   %%   From cntnm/build directory, type: gmake -f make_cntnm_sergio_run8 <TARGET>
   %% <<<<< these steps are all done by make_sergio_run8.sc    these steps are all done by make_sergio_run8.sc  >>>>>
   
   runner = ['!time cntnm_sergio_run8_v3.2_linux_gnu_dbl < ' fin];
   
   eval(runner);
   sedder = ['!sed ''1,73d'' ' fout ' > ' xfout];
   eval(sedder);
   boo = load(xfout);
   boo = boo(1:dalen0*5,:);
   fr = boo(:,1);
   od(:,ii) = boo(:,2);

   rmer = ['!/bin/rm ' fin ' ' fout ' ' xfout];
   eval(rmer)
end

addpath /home/sergio/SPECTRA
od5 = boxint_many(od',5);
fr5 = boxint_many(fr,5);

% min(fr)
% min(fr5)
% max(fr)
% max(fr5)
