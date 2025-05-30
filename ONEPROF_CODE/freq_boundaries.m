dirout = ['/asl/s1/sergio/H2012_RUN8_NIRDATABASE/ONE_PROFILE_IR_605_2830/g' num2str(gid) '.dat'];
dirout = ['/asl/s1/sergio/H2020_RUN8_NIRDATABASE/ONE_PROFILE_IR_605_2830/g' num2str(gid) '.dat'];

topts = runXtopts_params_smart(2000); 
dv = topts.ffin*nbox*pointsPerChunk;

wn1 = 605;
wn2 = 2855-dv;   %% when checking against Howards results

wn1 = 1405;
wn2 = 1430-dv;   %% when checking against Howards results

wn1 = 1205;
wn2 = 1405-dv;   %% when checking against Howards results

wn1 = 1405;
wn2 = 1605-dv;   %% when checking against Howards results

wn1 = 1205;
wn2 = 1605-dv;   %% when checking against Howards results

%% these may be overwritten by the code that actually calls this
%% subroutine
fmin = wn1; 
fmax = wn2;
