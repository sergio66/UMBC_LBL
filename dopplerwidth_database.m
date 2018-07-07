% function [y] = dopplerwidth_database(v0,T,m)

figure(1); clf;

%% eg this should show 0.0025 cm-1 spacing for 0800-2800 cm-1
%% eg this should show 0.0015 cm-1 spacing for 0500-0800 cm-1

v0 = 0.1:0.1:4000;    %% original only upto NIR
v0 = 0.1:0.1:20000;   %% new upto VIS
v0 = 0.1:0.1:40000;   %% new upto VIS/UV
T  = 200;
m  = 48;

% optimized for dopplerwidth_database(0.1:0.1:4000,150,50);
% optimized for dopplerwidth_database(0.1:0.1:4000,100,50); is also good!
% EVRYTHING HERE DONE at input mono spacing ie BEFORE 5 point boxcar
%this function computes the doppler width for spieces mass m, temp T, 
% at center freq v0
%line center  = v0 
%atmomic mass = m (amu) 
%temperature  = T  (k)

k       = 1.380658e-23; 
c_light = 2.99792458e8;         %ms-1 
amu     = 1.6605402e-27;        %nucleon mass/kg 
mass    = m*amu;

y =  v0.*sqrt(2*log(2)*k*T/mass/c_light/c_light);
plot(v0,y); xlabel('v0 cm-1'); ylabel('\delta \nu (doppler)');  pause

%%%%%%first find where the MONOCHROMATIC (before 5 pt boxcar) = 0.0005 cm-1
ii    = find(abs(y - 5e-4) <= 1e-6);
v0avg = sum(v0(ii))/length(ii);
v0avg = 1.1*v0avg;
jj0   = ii(1);            %%%%%this is where we start searching

currentchunks = 605:25:2805;
differ = abs(v0avg - currentchunks);
mn = min(differ);
iii = find(differ == mn);

datastart = currentchunks(iii);
fprintf(1,'monochromatic width = 5e-3 at v0 = %12.6f \n',v0avg);
fprintf(1,'closest chunk beginning to v0 = %12.6f \n',datastart);

%%%% start going to the right of jj0, in 10000 pts %%%%%%%%%%%%%%%%%%%%%%%%
%%%% at datastart, monochromatic = 0.0005 ==> 5 pt avg = 0.0025

%first restart jj0
ii    = find(abs(datastart - v0) <= 1e-6); 
jj0   = ii(1);            %%%%%this is where we start searching 

startstep = 5e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii=1; stepsize(ii) = startstep; chunkstart(ii) = datastart;  
fprintf(1,'index = %3i  stepsize = %12.6f chunkstart = %12.6f \n',ii,...
	startstep*1000*10000,datastart);

disp('going ++++ ...')
jj = jj0;
while (jj < length(v0))
  startstep = startstep*2;
  while ((jj <= length(v0)) & (y(jj) <= startstep))
    jj = jj+1;   %%recall y(jj) is the doppler width
    end

  if (jj < length(v0))
    currentchunks = chunkstart(ii):10000*stepsize(ii):5000;
    currentchunks = chunkstart(ii):10000*stepsize(ii):max(v0);
    v0avg = v0(jj);
    differ = abs(v0avg - currentchunks);
    mn = min(differ);
    iii = find(differ == mn);

    datastart = currentchunks(iii);

    ii = ii+1; stepsize(ii) = startstep; chunkstart(ii) = datastart;  
    fprintf(1,'index = %3i  stepsize = %12.6f chunkstart = %12.6f \n',ii,...
	startstep*1000*10000,datastart);

    end

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% start going to the left of jj0, in 10000 pts 
%%%% at datastart, monochromatic = 0.0005 ==> 5 pt avg = 0.0025
disp('going ---- ...')

startstep = 5e-4;
jj = jj0;
while ((jj >= 1) & (chunkstart > 40.0))
  startstep = startstep/2;
  while ((jj >= 1) & (y(jj) >= startstep))
    jj = jj-1;
    end

  if (jj >= 1)
    currentchunks = 0.0:10000*stepsize(ii):chunkstart(ii);
    v0avg = v0(jj);
    differ = abs(v0avg - currentchunks);
    mn = min(differ);
    iii = find(differ == mn);

    datastart = currentchunks(iii);

    ii = ii+1; stepsize(ii) = startstep; chunkstart(ii) = datastart;  
    fprintf(1,'index = %3i  stepsize = %12.6f chunkstart = %12.6f \n',ii,...
	startstep*1000*10000,datastart);

    end

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,' \n');
fprintf(1,'   i   chunkstart  (monostep)*1e-7   (5*monostep)*1e-7  chunksize   chunkstop    numchunks\n');
fprintf(1,'----------------------------------------------------------------------------------------\n');

barr=zeros(size(v0));
[yy,ii] = sort(chunkstart); 
  chunkstart = chunkstart(ii); stepsize = stepsize(ii);
jj_start = find(v0 <= chunkstart(1));  jj_start=jj_start(length(jj_start));
for ii=1:length(chunkstart)
  if (ii < length(chunkstart))
    numkc(ii) = (chunkstart(ii+1) - chunkstart(ii))/(5*stepsize(ii)*10000);
  else  
    numkc(ii) = (v0(length(v0)) - chunkstart(ii))/(5*stepsize(ii)*10000);
    end
  numkc(ii) = ceil(numkc(ii));
  chunkend(ii)   = chunkstart(ii) + 10000*5*stepsize(ii)*numkc(ii);
  chunksize(ii) = 10000*5*stepsize(ii); %%5 pt avg, 10000 pts/chunk

  fprintf(1,' %3i   %11.4f    %11.4f    %11.4f    %11.4f    %11.4f    %3i\n',ii,...
  chunkstart(ii),stepsize(ii)*1000*10000,5*stepsize(ii)*1000*10000,chunkend(ii),chunksize(ii),numkc(ii));

  %plot the doppler line
  v00=chunkstart(ii); dv=stepsize(ii); v=v00-5*dv:dv:v00+5*dv;
  y00=line_doppler(v,v00,T,m,1); plot(v,y00,'+',v,y00); grid; pause

  %do the bar graph info
  if (ii < length(chunkstart))
    jj_end = find(v0 <= chunkstart(ii+1)); jj_end=jj_end(length(jj_end));
    barr(jj_start:jj_end) = dv;
    jj_start = jj_end+1;
  else
    jj_end = length(v0); 
    barr(jj_start:jj_end) = dv;
    jj_start = jj_end+1;
    end
  end
monospace = stepsize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(' ')
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
disp('Scott/Sergio preference on May 2007-8')
disp('0.00250 cm^-1 (0.00200 x 5): 100.0 cm^-1 chunks: 3500 to 5500 cm^-1 ')
disp('0.00250 cm^-1 (0.00050 x 5):  25.0 cm^-1 chunks:  800 to 3500 cm^-1 ')
disp('0.00150 cm^-1 (0.00030 x 5):  15.0 cm^-1 chunks:  500 to 800 cm^-1 ');
disp('0.00100 cm^-1 (0.00020 x 5):  10.0 cm^-1 chunks:  300 to 500 cm^-1 ');
disp('0.00050 cm^-1 (0.00010 x 5):   5.0 cm^-1 chunks:  140 to 300 cm^-1 ');
disp('0.00025 cm^-1 (0.00005 x 5):   2.5 cm^-1 chunks:   80 to 140 cm^-1 ');
disp('0.00015 cm^-1 (0.00003 x 5):   1.5 cm^-1 chunks:   50 to  80 cm^-1 ');
disp('0.00010 cm^-1 (0.00002 x 5):   1.0 cm^-1 chunks:   30 to  50 cm^-1 ');
disp('0.00005 cm^-1 (0.00001 x 5):   0.5 cm^-1 chunks:   14 to  30 cm^-1 ');
disp(' ')
disp('and so on (same pattern with different factor of 10) down into ');
disp('the microwave region.  It would be 3x more work to extend      ');
disp('the database down to 1 cm^-1 than to stop around 140 cm^-1.    ');

scott_dv = ...
  [0.00200 0.00050 0.00030 0.00020 0.00010 0.00005 0.00003 0.00002 0.00001];
scott_v0 = ...
  [3500    800     500     300     140     80      50      30      14];
sarr=zeros(size(v0));
[yy,ii] = sort(scott_v0); scott_v0=scott_v0(ii); scott_dv=scott_dv(ii);
jj_start = find(v0 <= scott_v0(1));  jj_start=jj_start(length(jj_start));
for ii=1:length(scott_v0)
  if (ii < length(scott_v0))
    numkc(ii) = (scott_v0(ii+1) - scott_v0(ii))/(5*scott_dv(ii)*10000);
  else  
    numkc(ii) = (v0(length(v0)) - scott_v0(ii))/(5*scott_dv(ii)*10000);
    end
  numkc(ii) = ceil(numkc(ii));
  fprintf(1,' %3i   %12.6f      %12.6f       %3i\n',ii,...
       scott_v0(ii),scott_dv(ii)*5,numkc(ii));

  %plot the doppler line
  v00=scott_v0(ii); dv=scott_dv(ii); v=v00-5*dv:dv:v00+5*dv;
  y00=line_doppler(v,v00,T,m,1); plot(v,y00,'+',v,y00); grid; pause

  %do the bar graph info
  if (ii < length(scott_v0))
    jj_end = find(v0 <= scott_v0(ii+1)); jj_end=jj_end(length(jj_end));
    sarr(jj_start:jj_end) = dv;
    jj_start = jj_end+1;
  else
    jj_end = length(v0); 
    sarr(jj_start:jj_end) = dv;
    jj_start = jj_end+1;
    end
  end

semilogx(v0,y,v0,5e-4*ones(size(v0)),'b.-',chunkstart,monospace,'ro',...
         v0,barr,'r',scott_v0,scott_dv,'ko-',v0,sarr,'k'); 
plot(v0,y,v0,5e-4*ones(size(v0)),'b.-',chunkstart,monospace,'ro',...
         v0,barr,'r',scott_v0,scott_dv,'ko-',v0,sarr,'k'); 
grid; title('Monochromatic Doppler Widths'); 
xlabel('v0 cm-1'); ylabel('mono width (before boxcar) cm-1');

figure(2); rad = ttorad(v0,300);
           radSun = ttorad(v0,6000)*6.785087652174316e-5;
  h1 = subplot(211); plot(v0,rad,v0,radSun); ylabel('ttorad(v,300)')
  h2 = subplot(212); plot(v0,cumsum(rad)/max(cumsum(rad)));
                                   ylabel('cumsum(rad)')
  adjust21(h1,h2,'even')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
