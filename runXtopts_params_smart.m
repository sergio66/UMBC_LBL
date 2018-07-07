function topts = runXtopts_params_smart(fmin0,iAlways_1_2_25)

%% based on fmin0, this subroutine suggests 
%%   ffin, fmed, fcor, fstep,  xnear, xmed, xfar based on "dopplerwidth_database"
%% typically nbox = 5, pointsPerChunk = 10000)
%% stuff for fmin >= 800 cm-1 is taken from run7* defaults
%% iAlways_1_2_25 = -1 (default) so that if current box edges are [f1 f2] then
%%                                       xnear ~ box of                     [-ffin+f1     f2+ffin]
%%                                       xmed  ~ boxes of [-2*ffin+f1    -1*ffin+f1] and [f2+1*ffin f2+2*ffin]
%%                                       xfar  ~ boxes of [-25-2*ffin+f1 -2*ffin+f1] and [f2+2*ffin f2+2*ffin+25] 
%%
%%                  +1 (force overkill) so that if current box edges are [f1 f2] then
%%                                       xnear ~ box of           [-1+f1     f2+1]
%%                                       xmed  ~ boxes of [-2+f1    -1+f1] and [f2+1 f2+2]
%%                                       xfar  ~ boxes of [-25-2+f1 -2+f1] and [f2+2 f2+2+25] 

%% note Scott had suggested for v0 = 500.0,50.0 to use ffin = 0.0003,0.00003
%% but this seems to catch "checkpar" at ---> if (isint(xmed/fmed) ~= 1) <----

%% default behaviour is for 605-2830 cm-1 to have 
%%     ffin: 5.0000e-04
%%     fmed: 0.1000
%%     fcor: 0.5000
%%    fstep: 1
%%    xnear: 1
%%     xmed: 2
%%     xfar: 25

%% finally goes through  "checkpar.m"

if nargin == 1
  iAlways_1_2_25 = -1;
end

topts.ffin = 5.0000e-04;
topts.fmed = 0.1000;
topts.fcor = 0.5000;
topts.fstep = 1.0;
topts.xnear = 1.0;
topts.xmed  = 2.0;
topts.xfar  = 25.0;

%if (fmin0 >= 605 & fmin0 < 2830)
if (fmin0 >= 605 & fmin0 < 800)
  fmin = 900;
  fprintf(1,'for original kCARTA database, using 0.0005 cm-1 x 5 %6i\n',fmin0)
elseif (fmin0 >= 1400 & fmin0 < 2830)
  fmin = 900;
  fprintf(1,'for original kCARTA database, using 0.0005 cm-1 x 5 %6i\n',fmin0)
else
  fmin = fmin0;
end

nbox = 5;
pointsPerChunk = 10000;

%% Scotts B
v00    = [800.0   500.0   300.0   140.0  ];
ffin0  = [5.00e-4 3.00e-4 2.00e-4 1.00e-4];

v0 = [];
ffin = [];

%%go from 0.1 um to 10000 um
%xx = 2 : -1 : 0;
xx = 4 : -1 : -2;
xx = 10.^xx;
for ii = 1 : length(xx)
  mult = xx(ii)/100;
  v00*mult;
  v0   = [v0   v00*mult];
  ffin = [ffin ffin0*mult];
end

iPrint = -1;
if iPrint > 0
  disp('fcenter (cm-1)       10000 pt chunksize (cm-1)      fcenter (um)');
  disp('-----------------------------------------------------------------');
  aa = [v0; ffin*nbox*pointsPerChunk; 10000./v0];
  aa = fliplr(aa);
  boo = find(aa(1,:) <= 50); 
  fprintf(1,'  %10.4f          %10.4f              %10.4f\n',aa(:,boo));
  disp('-----------------------------------------------------------------');
  boo = find(aa(1,:) > 50 & aa(1,:) <= 15000); 
  fprintf(1,'  %10.4f          %10.4f              %10.4f\n',aa(:,boo));
  disp('-----------------------------------------------------------------');
  boo = find(aa(1,:) > 15000); 
  fprintf(1,'  %10.4f          %10.4f              %10.4f\n',aa(:,boo));
end

fcor  = ffin*1000;
fmed  = fcor/5;
%%% --->>>>>>>>>>>>>>>>>>>>>>>>>

fstep = ffin*nbox*pointsPerChunk/25;
xnear = fstep;
xmed  = fstep*2;
%%xfar  = ffin*nbox*pointsPerChunk;
xfar  = ones(size(fstep))*25.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v0 = fliplr(v0);
ffin = fliplr(ffin);
fmed = fliplr(fmed);
fcor = fliplr(fcor);
fstep = fliplr(fstep);
xnear = fliplr(xnear);
xmed  = fliplr(xmed);
xfar  = fliplr(xfar);

if fmin < min(v0)
  error('oops fmin < min(v0)!!!!');
end

ok = find(fmin >= v0);
ok = ok(length(ok));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
topts.ffin = ffin(ok);     %% fine   point spacing, before the nbox(=5) avg
topts.fmed = fmed(ok);     %% med    point spacing, before the nbox(=5) avg
topts.fcor = fcor(ok);     %% coarse point spacing, before the nbox(=5) avg
topts.fstep = fstep(ok);   %% divide 10000 pts chunk into boxes of this width

fprintf(1,'In runXtopts_params_smart.m, you sent in a "wavenumber" = %8.2f cm-1 \n',fmin0);
fprintf(1,'  \n')
fprintf(1,'    ffin = fine point spacing nbox(=5) avg = %10.6f cm-1 \n',topts.ffin)
fprintf(1,'    fmed = med  point spacing nbox(=5) avg = %10.6f cm-1 \n',topts.fmed)
fprintf(1,'    fcor = cor  point spacing nbox(=5) avg = %10.6f cm-1 \n',topts.fcor)
fprintf(1,' so we would typically divide [fmin,fmax] into boxes of width = %10.5f cm-1 \n',topts.fstep);
fprintf(1,' while 10000 pts boxcar=5 avg would be %10.5f cm-1 \n',topts.ffin*5*10000);
fprintf(1,'  \n')

%%%%%%%%%%%%%%%%%%%%%%%%%

%%% this is BEFORE JAN 2014 code
topts.xnear = xnear(ok);        %% near wing distance
topts.xmed  = xmed(ok);         %% med wing distance
topts.xfar  = max(xmed(ok),25); %% far wing distance   %%% <<<<--- should have been xfar

fprintf(1,' based on this, we have before JAN 2014 (A) \n')
fprintf(1,'    near wing distance = %10.6f cm-1 \n',topts.xnear);
fprintf(1,'    med wing distance  = %10.6f cm-1 \n',topts.xmed);
fprintf(1,'    far wing distance  = %10.6f cm-1 \n',topts.xfar);

%%%%%%%%%%%%%%%%%%%%%%%%%

if iAlways_1_2_25 == +1
  %%% this is AFTER JAN 2014 code
  topts.xnear = min(max(xnear(ok),1),1);        %% near wing distance
  topts.xmed  = min(max(xmed(ok),2),2);         %% med wing distance
  topts.xfar  = min(max(xfar(ok),25),25);       %% far wing distance   

  fprintf(1,' AFTER Jan 2014, FORCING the near,med,far boxes to the following adjustments (B) \n')
  if topts.xnear < topts.fstep
    fprintf(1,'  need xnear >= fstep, so adjust xnear from %10.6f to %10.6f cm-1 \n',topts.xnear,topts.fstep)
    topts.xnear = topts.fstep;
  end

  if ~isint(topts.xnear/topts.ffin)
    boo = ceil(topts.xnear/topts.ffin)*topts.ffin;
    fprintf(1,'  to make number of points an integer, adjusting xnear from %10.6f to %10.6f cm-1 \n',topts.xnear,boo)
    topts.xnear = boo;
  end

  if topts.xmed < topts.fstep
    fprintf(1,'  need xmed >= fstep, so adjust xnear from %10.6f to %10.6f cm-1 \n',topts.xmed,2*topts.fstep)
    topts.xmed = 2*topts.fstep;
  end

  if ~isint(topts.xmed/topts.fmed)
    boo = ceil(topts.xmed/topts.fmed)*topts.fmed;
    fprintf(1,'  to make number of points an integer, adjusting xmed  from %10.6f to %10.6f cm-1 \n',topts.xmed,boo)
    topts.xmed = boo;
  end

  fprintf(1,'    near wing distance = %10.6f cm-1 \n',topts.xnear);
  fprintf(1,'    med wing distance  = %10.6f cm-1 \n',topts.xmed);
  fprintf(1,'    far wing distance  = %10.6f cm-1 \n',topts.xfar);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fmax = fmin + topts.ffin*nbox*pointsPerChunk;   %%assume this

z = checkpar(topts.fstep,topts.ffin,topts.fmed,nbox,topts.fcor,...
             fmax,fmin,...
             topts.xnear,topts.xmed,topts.xfar);

