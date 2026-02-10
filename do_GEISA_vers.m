%%% function x = do_GEISA_vers()

%% geisapath.m is a function that returns the geisa path eg /asl/data/geisa/ or /umbc/xfs3/strow/asl/rta/geisa/

%% please update this for latest GEISA version (G15 in Feb 2019)

%% run8, run8water, calc_xsec use this as topts to override the default GEISA
%% string    HITRAN         path to HITRAN database   /asl/data/hitran/h08.by.gas
%%           GEISA          path to GEISA  database   /asl/data/geisa/g15.by.gas

GEISApathNyear = '/asl/data/geisa/g15.by.gas/';
GEISApathNyear = '/umbc/xfs3/strow/asl/rta/geisa/g15.by.gas/';
GEISApathNyear = [geisapath '/g15.by.gas/'];

%%% x = GEISApathNyear;
