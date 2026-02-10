%% hitranpath.m is a function that returns the hitran path eg /asl/data/hitran/ or /umbc/xfs3/strow/asl/rta/hitran/

%% please update this for latest HITRAN version (HH24 in Feb 2026)

%% run8, run8water, calc_xsec use this as topts to override the default HITRAN
%% string    HITRAN         path to HITRAN database   /asl/data/hitran/h08.by.gas
%%                          path to GEISA  database   /asl/data/geisa/g15.by.gas

HITRAN        = '/asl/data/hitran/h2k.by.gas/';
HITRAN        = '/asl/data/hitran/h08.by.gas/';
HITRAN        = '/asl/data/hitran/h12.by.gas/';
HITRAN        = '/asl/data/hitran/h16.by.gas/';
HITRAN        = '/asl/rta/hitran/h16.by.gas/';
HITRAN        = '/asl/rta/hitran/h20.by.gas/';
HITRAN        = '/umbc/xfs3/strow/asl/rta/hitran/h20.by.gas/';
HITRAN        = '/umbc/xfs3/strow/asl/rta/hitran/h24.by.gas/';
HITRAN        = [hitranpath '/h12.by.gas/'];
HITRAN        = [hitranpath '/h16.by.gas/'];
HITRAN        = [hitranpath '/h20.by.gas/'];
HITRAN        = [hitranpath '/h24.by.gas/'];
