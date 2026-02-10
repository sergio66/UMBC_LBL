%% this is default
[w3,d3_24,hitranfile24] = run8(3,980,1080,'IPFILES_EXAMPLE/test_one.txt');
[w1,d1_24,hitranfile24] = run8water(1,1280,1580,'IPFILES_EXAMPLE/test_one.txt');

%% now artificially set topts.HITRAN to use H2024
topts.HITRAN = [hitranpath '/h24.by.gas/'];
[w3,xd3_24,xhitranfile24] = run8(3,980,1080,'IPFILES_EXAMPLE/test_one.txt',topts);
[w1,xd1_24,xhitranfile24] = run8water(1,1280,1580,'IPFILES_EXAMPLE/test_one.txt',topts);

%% now artificially set topts.HITRAN to use H2020
topts.HITRAN = [hitranpath '/h20.by.gas/'];
[w3,xd3_20,xhitranfile20] = run8(3,980,1080,'IPFILES_EXAMPLE/test_one.txt',topts);
[w1,xd1_20,xhitranfile20] = run8water(1,1280,1580,'IPFILES_EXAMPLE/test_one.txt',topts);

figure(1);semilogy(w1,xd1_24,'b.-',w1,xd1_20,'g.-',w1,d1_24,'r');
  title('WV'); legend('x24','x20','default 24','location','best');
figure(2); plot(w1,xd1_24./d1_24,'b.-',w1,xd1_20./d1_24,'g');
  title('WV'); legend('x24/24','x20/24','location','best');
disp('ret to continue'); pause

figure(1);semilogy(w3,xd3_24,'b.-',w3,xd3_20,'g.-',w3,d3_24,'r');
  title('O3'); legend('x24','x20','default 24','location','best');
figure(2); plot(w3,xd3_24./d3_24,'b.-',w3,xd3_20./d3_24,'g');
  title('O3'); legend('x24/24','x20/24','location','best');
disp('ret to continue'); pause

fprintf(1,' hitranfile24 = %s \n', hitranfile24)
fprintf(1,'xhitranfile24 = %s \n',xhitranfile24)
fprintf(1,'xhitranfile20 = %s \n',xhitranfile20)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath READ_XSEC
%% [iYes,gf] = findxsec_plot(fmin,fmax,gasid,20XY);
[iYes,gfx20] = findxsec_plot(780,1080,51,2020);
[iYes,gfx24] = findxsec_plot(780,1080,51,2024);
fprintf(1,'gfx20 = %s \n',gfx20)
fprintf(1,'gfx24 = %s \n',gfx24)

%% see /home/sergio/git/HITRAN2UMBCLBL/MAKEIR/H2020/MAKEIR_ALL_H20/clust_runXtopts_savexsecN_file.m
%% [d,w] = calc_xsec(gasid,fmin,fmax-dvv,dvv,tp,pL,figno,2020);
boo = load('IPFILES_EXAMPLE/test_one.txt'); tp = boo(4); pL = boo(2);
[xgid51_20,w51,xgfout20] = calc_xsec(51,780,1080-0.0025,0.0025,tp,pL,1,2020);
[xgid51_24,w51,xgfout24] = calc_xsec(51,780,1080-0.0025,0.0025,tp,pL,2,2024);
[gid51_24,w51,gfout24] = calc_xsec(51,780,1080-0.0025,0.0025,tp,pL,2);

fprintf(1,' xsecfile24 = %s \n',gfout24)
fprintf(1,'x_secfile24 = %s \n',xgfout24)
fprintf(1,'x_secfile20 = %s \n',xgfout20)

figure(1);semilogy(w51,xgid51_24,'b.-',w51,xgid51_20,'g.-',w51,gid51_24,'r');
  title('CFC11'); legend('x24','x20','default 24','location','best');
figure(2); plot(w51,xgid51_24./gid51_24,'b.-',w51,xgid51_20./gid51_24,'g');
  title('CFC11'); legend('x24/24','x20/24','location','best');
disp('ret to continue'); pause
