
% step thru all xsec gas IDS and bands, and display
% all tabulated absorptions (in yellow) together with
% selected extrapolations, as set below.  Mainly a
% quick sanity check of extrapolation.
%
% for easier-to-use inspection of small extrapolation
% sets, use browse_xsec

% sample test profile
et = [300, 280, 260, 280, 290];	  % test temperatures
ep = [1,  300, 600, 850, 1000];   % test pressures
% et = [260, 300];
% ep = [10, 850];

% loop on gas IDs
for gf = 51:81
  disp(' ')
  disp(' ')
  disp(' ')

  xs = read_xsec(gf);
  
  [nrec, nband] = size(xs);
  
  % loop on bands
  for b = 1:nband

    fprintf(1, 'band %d\n', b);
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % scatter plot of (t,p) pairs for this band
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
    clf
    plot([xs(:,b).temp], [xs(:,b).pres], 'o');
    axis([200,310,0,1000])
    xlabel('temperature')
    ylabel('pressure')
    title(sprintf('%s sample points, band %d', xs(1,b).gstr, b))
    grid

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot tabulated and selected absorption spectra 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    clf
    % clist = 'ymcrgbkymcrgbkymcrgbkymcrgbkymcrgbk';
    clist = {'r','g','b','k','r--','g--','b--','k--'};

    % plot test profile spectra for this band
    ev1 = xs(1, b).v1;
    ev2 = xs(1, b).v2;
    edv = 0.05;
    nevpts = 1 + round((ev2-ev1)/edv);  % points in output grid
    evgrid = ev1 + (0:nevpts-1) * edv;  % output wavenumber grid
    if length(et) > 0
      fprintf(1, '  calling calc_xsec() . . .\n');
      absc = calc_xsec(gf, ev1, ev2, edv, et, ep);
    end
    figure(2)
    for j = 1:length(et)
      plot(evgrid, absc(:,j), clist{j});
      hold on
    end

    % plot all tabulated spectra for this band
    for j = 1:nrec
      v1 = xs(j, b).v1;
      v2 = xs(j, b).v2;
      npts = xs(j, b).npts;
      dvp = (v2 - v1) / (npts - 1);
      plot(v1:dvp:v2, xs(j, b).absc, 'y');
      hold on
    end

    legend(num2str([ep', et']));

    title(sprintf('%s absorption coefficients, band %d', xs(1,b).gstr, b));
    ylabel('kmoles/cm^2')
    xlabel('wavenumber')
    grid
    hold off

    fprintf(1,'  <enter> for next band\n');
    pause

  end % band loop

end % gas loop

