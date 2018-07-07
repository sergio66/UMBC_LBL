function xsectab25(vfirst, vlast, gdir, cdir, refpro);

% function xsectab25(vfirst, vlast, gdir, cdir, refpro);
%
% abstab25 -- simple driver for tabulation of xsec absorptions
%
% inputs
%   vfirst  - wavenumber of first chunk 
%   vlast   - wavenumber of last chunk
%   gdir    - temp. dir. for uncompressed absorptions
%   cdir    - output dir. for compressed absorptions
%   refpro  - optional reference profile
%
% outputs - absorption data in 25 1/cm chunks, in .mat files of
% the form g<gid>v<v1>.mat, where <v1> is first wavenumber for
% that chunk

% set defaults
if nargin < 1,  vfirst = 605;  end
if nargin < 2,  vlast = 2805;  end
if nargin < 3,  gdir = '/home/motteler/absdat/abs.xsec';  end
if nargin < 4,  cdir = '/home/motteler/absdat/kcomp.xsec';  end
if nargin < 5,  refpro = '/home/motteler/abscmp/refpro';  end

% minimum significant absorption value
abseps = 1e-8;

% load reference profile
eval(sprintf('load %s', refpro))

% use only xsec gasses from reference profile
xgind = find(refpro.glist > 50);
refpro.glist = refpro.glist(xgind);
refpro.gamnt = refpro.gamnt(:,xgind);
refpro.gpart = refpro.gpart(:,xgind);

% tabulation frequency spacing
dvk = 0.0025;

% temperature offsets
toffset = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50];

% space for absorption chunks
k = zeros(1e4, 100, length(toffset));

% loop on profile xsec gasses
for gind = 1 : length(refpro.glist)

  gid = refpro.glist(gind);

  % use reference profile values
  tref = refpro.mtemp;		 % reference temperature profile
  pref = refpro.mpres;		 % reference profile combined pressures
  gref = refpro.gamnt(:,gind);	 % reference profile gas amount
  gref2d = ones(1e4,1) * gref';  % gas amount (duplicate rows)

  % loop on 25 1/cm chunks
  for v1 = vfirst : 25 : vlast
    v2 = v1 + 25 - dvk;

    fout = sprintf('g%dv%d.mat', gid, v1);

    % check this interval for any absorption
    abstest = calc_xsec(gid, v1, v2, 0.2, 250, 850);
    if sum(abstest) > 0 
  
      fprintf(1, 'gas %d, v1 = %g...\n', gid, v1);

      % loop on temperature offsets
      for tind = 1:length(toffset);
  
        tref2 = tref + toffset(tind);
  
        % extrapolat absorptions from tabulated values
	k1 = calc_xsec(gid, v1, v2, dvk, tref2, pref);
  
        % scale absorptions by gas amounts
        k(:, :, tind) = k1 .* gref2d;
  
      end % temp offset loop
  
     if max(max(max(k))) > abseps

        % save this completed chunk
        fr = v1 + (0:9999) * dvk; % output wavenumber grid
        eval(sprintf('save %s/%s fr k gid', gdir, fout));
        clear k
  
        % compress the data we just saved
        B = absbasis(gid, gdir, v1);
        absbcmp(gid, gdir, cdir, v1, B);

        % clean up
        % delete([gdir,'/',fout]);

      end % abseps test
    end % gas action test
  end % chunk loop
end % gass ID loop

