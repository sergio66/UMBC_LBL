function topts = runXtopts_params(fmin)

%% based on fmin, this subroutine suggests 
%% ffin, fmed, fcor, fstep,  xnear, xmed, xfar based on "dopplerwidth_database"
%% typically nbox = 5, pointsPerChunk = 10000)
%% stuff for fmin >= 800 cm-1 is taken from run7* defaults

%% note Scott had suggested for v0 = 500.0,50.0 to use ffin = 0.0003,0.00003
%% but this seems to catch "checkpar" at ---> if (isint(xmed/fmed) ~= 1) <----

%% finally goes through  "checkpar.m"

nbox = 5;
pointsPerChunk = 10000;

%%% --->>>>>>>>>>>>>>>>>>>>>>>>>
%% Sergio
v0    = [800.0   500.0   300.0   140.0   80.0    50.0    30.0    14.0   ...
                                         8.000   5.000   3.000   1.400  ];
ffin  = [5.00e-4 2.50e-4 2.00e-4 1.00e-4 5.00e-5 2.50e-5 2.00e-5 1.00e-5...
                                         5.00e-6 2.50e-6 2.00e-6 1.00e-6];
fmed  = [1.00e-1 1.00e-1 5.00e-2 2.50e-2 1.00e-2 1.00e-2 5.00e-3 2.50e-3 ...
                                         1.00e-3 1.00e-3 5.00e-4 2.50e-4];
fcor  = [5.00e-1 2.50e-1 2.00e-1 1.00e-1 5.00e-2 2.50e-2 2.00e-2 1.00e-2 ...
                                         5.00e-3 2.50e-3 2.00e-3 1.00e-3];

%% Scotts A
v0    = [800.0   500.0   300.0   140.0   80.0    50.0    30.0    14.0   ...
                                         8.000   5.000   3.000   1.400  ];
ffin  = [5.00e-4 3.00e-4 2.00e-4 1.00e-4 5.00e-5 3.00e-5 2.00e-5 1.00e-5...
                                         5.00e-6 3.00e-6 2.00e-6 1.00e-6];
fmed  = [1.00e-1 1.00e-1 5.00e-2 2.50e-2 1.00e-2 1.00e-2 5.00e-3 2.50e-3 ...
                                         1.00e-3 1.00e-3 5.00e-4 2.50e-4];
fcor  = [5.00e-1 3.00e-1 2.00e-1 1.00e-1 5.00e-2 3.00e-2 2.00e-2 1.00e-2 ...
                                         5.00e-3 3.00e-3 2.00e-3 1.00e-3];

%% Scotts B
v0    = [800.0   500.0   300.0   140.0   80.0    50.0    30.0    14.0   ...
                                         8.000   5.000   3.000   1.400  ];
ffin  = [5.00e-4 3.00e-4 2.00e-4 1.00e-4 5.00e-5 3.00e-5 2.00e-5 1.00e-5...
                                         5.00e-6 3.00e-6 2.00e-6 1.00e-6];

%% Scotts B .. out 605.0 instead of 800 cm-1 till we get JM Hartmann linemixing
v0    = [605.0   500.0   300.0   140.0   80.0    50.0    30.0    14.0   ...
                                         8.000   5.000   3.000   1.400  ];
ffin  = [5.00e-4 3.00e-4 2.00e-4 1.00e-4 5.00e-5 3.00e-5 2.00e-5 1.00e-5...
                                         5.00e-6 3.00e-6 2.00e-6 1.00e-6];

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

topts.ffin = ffin(ok);     %% fine   point spacing, before the nbox(=5) avg
topts.fmed = fmed(ok);     %% med    point spacing, before the nbox(=5) avg
topts.fcor = fcor(ok);     %% coarse point spacing, before the nbox(=5) avg
topts.fstep = fstep(ok);   %% divide 10000 pts chunk into boxes of this width
topts.xnear = xnear(ok);   %% near wing far distance
topts.xmed  = xmed(ok);    %% med wing far distance
topts.xfar  = xfar(ok);    %% far wing distance

fmax = fmin + topts.ffin*nbox*pointsPerChunk;   %%assume this

%[fmin fmax]

z = checkpar(topts.fstep,topts.ffin,topts.fmed,nbox,topts.fcor,...
             fmax,fmin,...
             topts.xnear,topts.xmed,topts.xfar);