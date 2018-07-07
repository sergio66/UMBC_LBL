
% xsec summary plot -- step thru all xsec gas IDS and display
% selected extrapolation values across a broad spectral range
%
% this is mainly to get an overview of where the various xsec
% regions are within a broader IR spectra

% bands
v1 = 500;
v2 = 1800;
dv = .05;

npts = (v2-v1)/dv + 1;

et = 290;  % sample temperature
ep = 850;  % sample pressure

h1 = figure(1); 
clf
clist = {'y',   'm',   'c',   'r',   'g',   'b',   'k',   ...
         'y--', 'm--', 'c--', 'r--', 'g--', 'b--', 'k--', ...
         'y-.', 'm-.', 'c-.', 'r-.', 'g-.', 'b-.', 'k-.', ...
         'y-s', 'm-s', 'c-s', 'r-s', 'g-s', 'b-s', 'k-s', ...
         'y-d', 'm-d', 'c-d', 'r-d', 'g-d', 'b-d', 'k-d', ...
         };

% loop on gas IDs
for gf = 51:81

  absbuf = zeros(npts,1);

  xs = read_xsec(gf);

  [nrec, nband] = size(xs);

  % loop on bands
  for b = 1:nband

    fprintf(1, 'calling calc_xsec, gid %d, band %d...\n', gf, b);
    absc = calc_xsec(gf, v1, v2, dv, et, ep);
    absbuf = absbuf + absc;

  end % band loop

  plot(v1:dv:v2, absbuf, clist{gf-50});
  axis([v1, v2, 0, 4e10])
  hold on

end % gas loop

hl = legend(num2str((51:81)')); set(hl,'fontsize',10);
title('IR Cross Section Overview');
ylabel('molecules/cm^2'); xlabel('wavenumber')
grid
hold off

orient landscape
% print -dpsc summary.ps
saveas(h1, 'summary.fig');

