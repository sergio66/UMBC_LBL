function browse_xsec(gid);

% Select and examine absorption spectra at selected 
% tabulated and interpolated or extrapolated values.
%
% This is the main tool for checking the sensibility
% of interpolated and extrapolated spectra
%

% gid = input('xsec gas ID > ');

xs = read_xsec(gid);

[nrec, nband] = size(xs);

for b = 1:nband

  fprintf(1, 'band %d\n', b);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % scatter plot of (t,p) pairs for this band
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(1); clf
  plot([xs(:,b).temp], [xs(:,b).pres], 'o');
  axis([200,310,0,1000])
  xlabel('temperature')
  ylabel('pressure')
  title(sprintf('%s (gid %d)  tabulation points, band %d', ...
		 xs(1,b).gstr, gid, b));
  grid

  % get points of interest
  fprintf(1, 'select tabulated point(s) to plot; <return> to continue\n');
  [st, sp] = ginput;
  
  % search xsec records for values closest to (st,sp)
  sr = zeros(length(st),1);
  for j = 1:length(st)
    dmin = 1000;
    for i = 1:nrec
      d = sqrt((xs(i,b).temp - st(j))^2 + (xs(i,b).pres - sp(j))^2);
      if d < dmin
        dmin = d;
        sr(j) = i;
      end
    end
  end

  hold on
  plot([xs(sr(:),b).temp], [xs(sr(:),b).pres], '*')
  text([xs(sr(:),b).temp]+1, [xs(sr(:),b).pres]+4, num2str((1:length(st))'))
  hold off

  fprintf(1, 'select extrapolated point(s) to plot; <return> to continue\n');
  [et, ep] = ginput;

  hold on
  plot(et, ep, '+')
  text(et+1, ep+4, num2str((length(st)+1:length(st)+length(et))'))
  hold off

  %%%%%%%%%%%%%%%%%%%%%%%%
  % plot selected spectra 
  %%%%%%%%%%%%%%%%%%%%%%%%
  figure(2); clf
  clist = {'y','r','g','b','k','y--','r--','g--','b--','k--'};

  % plot selected tabulated spectra
  for j = 1:length(st)
    v1 = xs(sr(j), b).v1;
    v2 = xs(sr(j), b).v2;
    npts = xs(sr(j), b).npts;
    dvp = (v2 - v1) / (npts - 1);
    plot(v1:dvp:v2, xs(sr(j), b).absc, clist{j});
    hold on
  end

  % plot selected extrapolated spectra
  ev1 = xs(1, b).v1;
  ev2 = xs(1, b).v2;
  edv = 0.01;
  nevpts = 1 + round((ev2-ev1)/edv);  % points in output grid
  evgrid = ev1 + (0:nevpts-1) * edv;  % output wavenumber grid
  if length(et) > 0
    fprintf(1, 'calling calc_xsec()...\n');
    absc = calc_xsec(gid, ev1, ev2, edv, et, ep, 3);
  end
  figure(2)
  for j = 1:length(et)
    plot(evgrid, absc(:,j), clist{j+length(st)});
    hold on
  end

  legend(num2str((1:length(st)+length(et))'))
  title('Selected Absorption Spectra');
  title(sprintf('%s (gid %d)  selected absorption spectra, band %d', ...
		 xs(1,b).gstr, gid, b));

  ylabel('kmoles/cm^2')
  xlabel('wavenumber')
  grid
  hold off

  if b < nband
    fprintf(1, '<return> for next band\n');
    pause
  end
end

