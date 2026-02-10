function x = geisapath;

%% this is called by do_GEISA_vers.m and sets variable GEISA
%% and other subroutines

%% run8, run8water, calc_xsec use this as topts to override the default HITRAN
%% string    HITRAN         path to HITRAN database   /asl/data/hitran/h08.by.gas
%%                          path to GEISA  database   /asl/data/geisa/g15.by.gas

%% also see hitranpath.m
  
str = 'boo';
str = '/asl/data/geisa/';
str = '/asl/rta/geisa/';
str = '/umbc/xfs3/strow/asl/rta/geisa/';

x = str;
