function x_with_unc = find_unc(unctype,uncvalue,X,dX_index,gid,fr0,iDebug);

%% input
%%   unctype  = 1   for wavenumber or pressure shift in cm-1
%%              2   for intensity, or self/foreign halfwidths or temp dependence
%%   uncvalue = +1  for use MAX uncertainty
%%              -1  for use MIN uncertainty = -1*MAX uncertainty
%%              +10 for use random value between MIN and MAX
%%              -10 for zero uncertainty ie keep things unchanged
%%   X        = line parameter to be adjusted
%%   dX_index = line.ai = uncertainty index associated with X
%%        gid = GasID            ... needed for random number seed
%%        fr0 = start wavenumber ... needed for random number seed
%% output
%%  x_with_unc = X + uncertainty from above

if nargin == 7
  iDebug = +1;
  iDebug = -1;    %% turn off debug
end  
iDebugSum = 0;

dX_index = dX_index + 1;   %% matlab arrays lowest index is 1 not 0, so need to offset

if uncvalue == 1
  delta = ones(size(X));
elseif uncvalue == -1
  delta = ones(size(X));
elseif uncvalue == -10
  delta = zeros(size(X));
elseif uncvalue == +10
  %rng('shuffle');
  
  %seed = [num2str(gid,'%02d') num2str(fr0,'%05d')];
  %seed = str2num(seed);
  %rng(seed,'twister');
  
  delta  = min(max(0.5*randn(size(X)),0),1);
  delta  = 2*(rand(size(X))-0.5);
  
  delta2 = randn(size(X));
  hist(delta2)

end

if abs(unctype) == 1
  %% wavenumber or pressure shift in cm-1
  themin(1) = 1;        themax(1) = 2;        %% code 0
  themin(2) = 1.0e-1;   themax(2) = 1.0e+0;   %% code 1
  themin(3) = 1.0e-2;   themax(3) = 1.0e-1;   %% code 2
  themin(4) = 1.0e-3;   themax(4) = 1.0e-2;   %% code 3      
  themin(5) = 1.0e-4;   themax(5) = 1.0e-3;   %% code 4
  themin(6) = 1.0e-5;   themax(6) = 1.0e-4;   %% code 5
  themin(7) = 1.0e-6;   themax(7) = 1.0e-5;   %% code 6      
  themin(8) = 1.0e-7;   themax(8) = 1.0e-6;   %% code 7
  themin(9) = 1.0e-8;   themax(9) = 1.0e-7;   %% code 8
  themin(10) = 1.0e-99; themax(10) = 1.0e-8;  %% code 9

  if unctype == -1
    themin(1) = 0;        themax(1) = 0.02;        %% code 0, tsp         Dec 05, 2017 old, keep as is
    themin(2) = 0;        themax(2) = 0.02;        %% code 1, tsp         May 15, 2018 new
    themin(3) = 0;        themax(3) = 0.02;        %% code 2, tsp         May 15, 2018 new    
  elseif unctype == +1
    %themin(1) = 0;        themax(1) = 1.0;         %% code 0, wavenumber  Dec 05, 2017 old
    themin(1) = 0;        themax(1) = 0.01;        %% code 0, wavenumber  Dec 15, 2017 new
    themin(2) = 0;        themax(2) = 0.01;        %% code 1, wavenumber  May 15, 2018 new
    themin(3) = 0;        themax(3) = 0.01;        %% code 2, wavenumber  May 15, 2018 new    
  end
  
  logthemin = log10(themin);
    bwah = find(themin < eps);
    logthemin(bwah) = log10(themax(bwah)/10);
  logthemax = log10(themax);
  
  for ii = 1 : length(themax)
    oo = find(dX_index == ii);
    if uncvalue == 1
      dX = +themax(ii)*delta(oo);
    elseif uncvalue == -1
      dX = -themax(ii)*delta(oo);
    elseif uncvalue == -10
      dX = 0*delta(oo);
    elseif uncvalue == +10
      %% old
      %gamma = (themax(ii)-themin(ii))/2;
      %dX = themin(ii) + (themax(ii)-themin(ii))/2 + delta(oo)*gamma;

      %% new
      gamma = (logthemax(ii)-logthemin(ii))/2;
      dX = logthemin(ii) + (logthemax(ii)-logthemin(ii))/2 + delta(oo)*gamma;
      dX = exp10(dX);      
      boo = find(delta2(oo) < 0);
      dX(boo) = -dX(boo);
    end
    if iDebug > 0
      iDebugSum = iDebugSum + length(oo);
      array = [unctype ii length(oo) length(X) iDebugSum];
      fprintf(1,'unctype ii length(oo) length(X) sumSOfar = %2i %2i %6i %6i %6i\n',array);
      figure(2); clf; plot(dX); 
      figure(1); clf; plot(X(oo)); pause
    end
    x_with_unc(oo) = X(oo) + dX;
  end
  
elseif unctype == 2
  %% intensity, or self/foreign halfwidths or temp dependence
  themin(1) = 0;        themax(1) = 0;        %% code 0, unknown/unavailable
  themin(2) = 0;        themax(2) = 0;        %% code 1, default or constant
  themin(3) = 0;        themax(3) = 0;        %% code 2, average or estimate
  themin(4) = 20;       themax(4) = 25;       %% code 3
  themin(5) = 10;       themax(5) = 20;       %% code 4
  themin(6) = 05;       themax(6) = 10;       %% code 5
  themin(7) = 02;       themax(7) = 05;       %% code 6
  themin(8) = 01;       themax(8) = 02;       %% code 7
  themin(9) = 00;       themax(9) = 01;       %% code 9

  themax(1:3) = 20;   %% codes 0,1,2  Dec 15, 2017 new
  if gid == 3 & abs(fr0-1000) < 250
    themax(1:3) = 3;   %% codes 0,1,2  from Brian Drouin email and JSQRT 2017 paper
  end

  themin = themin/100;
  themax = themax/100;

  for ii = 1 : length(themax)
    oo = find(dX_index == ii);
    if uncvalue == 1
      dX = +themax(ii)*delta(oo);
    elseif uncvalue == -1
      dX = -themax(ii)*delta(oo);
    elseif uncvalue == -10
      dX = 0*delta(oo);
    elseif uncvalue == +10
      gamma = (themax(ii)-themin(ii))/2;    
      dX = themin(ii) + (themax(ii)-themin(ii))/2 + delta(oo)*gamma;
      boo = find(delta2(oo) < 0);
      dX(boo) = -dX(boo);
    end
    if iDebug > 0
      iDebugSum = iDebugSum + length(oo);
      array = [unctype ii length(oo) length(X) iDebugSum];
      fprintf(1,'unctype ii length(oo) length(X) sumSOfar = %2i %2i %6i %6i %6i\n',array);
      figure(2); clf; plot(dX); 
      figure(1); clf; plot(X(oo)); pause
    end
    x_with_unc(oo) = X(oo).*(1 + dX);
  end

else
  unctype
  error('unknown unctype')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = exp10(x);

y = 10.^x;
