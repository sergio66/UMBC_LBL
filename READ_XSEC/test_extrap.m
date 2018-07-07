
% Basic test of temperature extrapolation 
%

% select gas ID
gf = 51;

% extrapolation parameter, Ke = Km*(Tm/Te)^C1, where 
%    Km is measured absorption at temperature Tm, and 
%    Ke is estimated absorption at temperature Te.
C1 = 0.75;

xs = read_xsec(gf);

[nrec, nband] = size(xs);

for b = 1:nband

  fprintf(1, 'band %d\n', b);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % scatter plot of (t,p) sample points
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(1)
  clf
  hold on
  for r = 1:nrec
    if isempty(xs(r,b)), break, end
    plot(xs(r,b).temp, xs(r,b).pres, '+');
  end
  xlabel('temperature')
  ylabel('pressure')
  title([xs(r,b).gstr, ' sample points'])
  hold off

  % get points of interest
  fprintf(1, 'select two points\n', b);
  [x,y] = ginput(2);
  
  % search for xsec records closest to ginputs
  q = zeros(length(x),1);
  for j = 1:length(x)
    min = 1000;
    for i = 1:nrec
      d = sqrt((xs(i,b).temp - x(j))^2 + (xs(i,b).pres - y(j))^2);
      if d < min
        min = d;
        q(j) = i;
      end
    end
  end

  hold on
  for j = 1:length(x)
    plot(xs(q(j),b).temp, xs(q(j),b).pres, 'o')
    text(xs(q(j),b).temp+1, xs(q(j),b).pres+4, num2str(j))
  end
  hold off
  
  %%%%%%%%%%%%%%%%%%%%%
  % extrapolation test
  %%%%%%%%%%%%%%%%%%%%%
  figure(2)
  clf
  clist = 'ymcrgbwkymcrgbwk';
  hold on
  j2 = length(x);
  j1 = j2-1;
  
  % extrapolate from record k1 to record k2:
  k1 = q(j1);
  k2 = q(j2);
  ax = xs(k1,b).absc * (xs(k1,b).temp/xs(k2,b).temp)^C1;
  dvp = (xs(k1,b).v2 - xs(k1,b).v1) / (xs(k1,b).npts - 1);
  f = xs(k1,b).v1:dvp:xs(k1,b).v2;
  plot(f, ax, clist(j2+1));
  
  % measured value
  dvp = (xs(k2,b).v2 - xs(k2,b).v1) / (xs(k2, b).npts - 1);
  f = xs(k2,b).v1:dvp:xs(k2,b).v2;
  plot(f, xs(k2,b).absc, clist(j2));
  title('Extrapolation and Measurement Compared');
  ylabel('kmoles/cm^2')
  xlabel('wavenumber')
  legend('extrapolated', 'measured');
  hold off
  
  fprintf(1, '<return> for next band\n');
  pause

end

return

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % absorption spectra for selected (t,p) points
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(3)
  clf
  clist = 'rgbwkymcrgbwkymc';
  hold on
  a = [];
  for j = 1:length(x)
    k = q(j);
    dvp = (xs(k,b).v2 - xs(k,b).v1) / (xs(k,b).npts - 1);
    f = xs(k,b).v1:dvp:xs(k,b).v2;
    plot(f, xs(k,b).absc, clist(j));
    a = [a, '''', num2str(j), ''', '];
  end
  eval(['legend(', a, '0)']); 
  title('Selected Absorption Spectra');
  ylabel('molecules/cm^2')
  xlabel('wavenumber')
  hold off
