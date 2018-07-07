function [absc, vgrid] = calc_xsec(gf, v1, v2, dv, tp, pL, db);

% function [absc, vgrid] = calc_xsec(gf, v1, v2, dv, tp, pL, db);
%
% calc_xsec returns IR absorptions calculated by
% interpolation or extrapolation of HITRAN tabulated
% IR cross-section data
%
% inputs
%   gf  - gas id or file
%   v1  - lower wavenumber
%   v2  - upper wavenumber
%   dv  - wavenumber increment
%   tp  - temperature profile
%   pL  - pressure levels for tp
%   db  - plot figure number (default is no plot)
%
% implicit inputs
%   nvpts = (v2-v1)/dv + 1, number of wavenumber points
%   nplev = length(pL), number of pressure levels
%
% outputs
%   absc  - nvpts x nplev array of layer absorptions, kmoles/cm^2
%   vgrid - nvpts vector of output frequencies
%

% Original Version:
% - H. Motteler, 15 Dec 98
% Revisions:
% - interpolation bug fix; fancier plotting; removed xsec data 
%   directory param; renamed "calc_xsec"; HM, 4 Apr 00

% physical extrapolation parameter, Ke = Km*(Tm/Te)^C1, where 
%   - Km is measured absorption at temperature Tm, and 
%   - Ke is estimated absorption at temperature Te.
%
C1 = 0.75;

% input defaults
if nargin == 6, 
  db = 0; 
end
if (nargin < 6 | 7 < nargin)
  error('wrong number of arguments')
end

% read_xsec() returns a structure with xsec data for this gas
xs = read_xsec(gf);
[nrec,nband] = size(xs);

% initializations
nvpts = 1 + round((v2-v1)/dv);  % points in output grid
vgrid = v1 + (0:nvpts-1) * dv;  % output wavenumber grid
v2 = vgrid(nvpts);              % make v2 a multiple of dv
nplev = length(pL); 	        % number of pressure levels
wc = zeros(nplev, nrec);        % "combination" weights
we = zeros(nplev, nrec);        % physical extrapolation weights
absc = zeros(nvpts, nplev);     % return zeros by default

for b = 1:nband

  % find intersection of xsec and requested wavenumber ranges
  % (xsec band edges may vary somewhat from record to record)
  v1C = max([v1, xs(:,b).v1]);
  v1C = ceil(v1C / dv) * dv;
  v2C = min([v2, xs(:,b).v2]);
  v2C = floor(v2C / dv) * dv;

  % see if the band is non-empty
  if v1C < v2C
  
    % wavenumber indices in output array
    v1Cind = interp1(vgrid, 1:nvpts, v1C, '*nearest');
    v2Cind = interp1(vgrid, 1:nvpts, v2C, '*nearest');
  
    % initalize intermediate absorption arrays
    % because xsec data may be tabulated at varying points even within 
    % a nominal band, we must first interpolate xsec data to the output 
    % grid spacing before taking weighted combinations
    absbuf = zeros(1+v2Cind-v1Cind, nrec);
   
    % find ratio of temperature and pressure ranges.  If the pressure
    % fields are 0, C2, the ratio of temperature and pressure ranges,
    % is set to 0.  (C2 is used to scale pressure, for extrapolation.)
    tmin = min([xs(:,b).temp]);
    tmax = max([xs(:,b).temp]);
    pmin = min([xs(:,b).pres]);
    pmax = max([xs(:,b).pres]);
    if pmax - pmin == 0
      C2 = 0;
    else
      C2 =  (tmax - tmin) / (pmax - pmin) ;
    end
  
    % loop on records and build interpolation weights
    for r = 1:nrec
  
      if isempty(xs(r,b).v1), break, end
  
      % xsec record wavenumber grid
      xdv = (xs(r,b).v2 - xs(r,b).v1) / (xs(r,b).npts - 1);
      xvgrid = xs(r,b).v1 + (0:xs(r,b).npts - 1) * xdv;
  
      % interpolate to output array spacing
      absbuf(:,r) = interp1(xvgrid, xs(r,b).absc, ...
			    vgrid(v1Cind:v2Cind)', '*linear');
  
      % compute weights for this record, for each pressure level
      for L = 1:nplev
  
        dt = (xs(r,b).temp - tp(L))^2;	  % temperature difference
        dp = (xs(r,b).pres - pL(L))^2;	  % pressure difference
        dp = C2^2 * dp;			  % balance t & p distance
  
        wc(L,r) = 1 / (dt + dp)^2;	  % combination weight
  
        we(L,r) = (xs(r,b).temp/tp(L))^C1;  % extrapolation weight
  
      end % pressure loop
  
    end % record loop
  
    % normalize combination weights, for each level
    for L = 1:nplev  
       wc(L,:) = wc(L,:) / sum(wc(L,:));
    end
  
    % combine absorptions, for each xsec record 
    for r = 1:nrec
      if isempty(xs(r,b).v1), break, end
      for L = 1:nplev  
  
        absc(v1Cind:v2Cind, L) = absc(v1Cind:v2Cind, L) + ...
  			     absbuf(:,r) * we(L,r) * wc(L,r);
      end
    end
  
    if db
      % option to visualize extrapolation weights
      figure(db)
      for L = 1:nplev
        clf
        stem3([xs(:,b).temp], [xs(:,b).pres], [wc(L,:)])
        hold on
        stem3(tp(L), pL(L), 0, 'r');
        xlabel('temperature');
        ylabel('pressure');
        if isnumeric(gf)
          spstr = '%s (gid %d)  extrapolation weights for %.0f K, %.0f mb';
          title(sprintf(spstr, xs(1,b).gstr, gf, tp, pL));
        else
          spstr = '%s extrapolation weights for %.0f K, %.0f mb';
          title(sprintf(spstr, xs(1,b).gstr, tp, pL));
        end
        hold off
        if L < nplev
          fprintf(1, '<return> for next extrapolation point\n');
          pause
        end
      end
    end % if db

  end % if v1C < v2C
  
end % band loop

