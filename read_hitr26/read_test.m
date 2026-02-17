% quick tests of Matlab vs C versions of read_hitran

test2012 = '/umbc/xfs3/strow/asl/rta/hitran/h12.by.gas/';
test2024 = '/umbc/xfs3/strow/asl/rta/hitran/h24.by.gas/';

%% function s = read_hitran2(v1, v2, str, gid, hsrc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1 = read_hitran(1000, 1002, 0, 2, test2012); 
s2 = read_hitran2(1000, 1002, 0, 2, test2012);
if isequal(s1, s2); disp('s1,s2 are equal'); else disp('OH OH s1,s2 are different'); end

s1 = read_hitran(1000, 1002, 0, 3, test2012);
s2 = read_hitran2(1000, 1002, 0, 3, test2012);
if isequal(s1, s2); disp('s1,s2 are equal'); else disp('OH OH s1,s2 are different'); end

s1 = read_hitran(1010, 1020, 0, 3, test2012);
s2 = read_hitran2(1010, 1020, 0, 3, test2012);
if isequal(s1, s2); disp('s1,s2 are equal'); else disp('OH OH s1,s2 are different'); end

s1 = read_hitran(1210, 1220, 0, 6, test2012);
s2 = read_hitran2(1210, 1220, 0, 6, test2012);
if isequal(s1, s2); disp('s1,s2 are equal'); else disp('OH OH s1,s2 are different'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic; s1 = read_hitran(0605,2830,  0, 6, test2024); toc
tic; s2 = read_hitran2(0605,2830, 0, 6, test2024); toc
if isequal(s1, s2); disp('s1,s2 are equal'); else disp('OH OH s1,s2 are different'); end

tic; s1 = read_hitran(0605, 2830,  0, 2, test2024); toc
tic; s2 = read_hitran2(0605, 2830, 0, 2, test2024); toc
if isequal(s1, s2); disp('s1,s2 are equal'); else disp('OH OH s1,s2 are different'); end

tic; s1 = read_hitran(0605, 2830,  0, 1, test2024); toc
tic; s2 = read_hitran2(0605, 2830, 0, 1, test2024); toc
if isequal(s1, s2); disp('s1,s2 are equal'); else disp('OH OH s1,s2 are different'); end

tic; s1 = read_hitran(0605, 2830,  0, 12, test2024); toc
tic; s2 = read_hitran2(0605, 2830, 0, 12, test2024); toc
if isequal(s1, s2); disp('s1,s2 are equal'); else disp('OH OH s1,s2 are different'); end

