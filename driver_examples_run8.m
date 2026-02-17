%% checks water without basement mex
[d1,w1] = run8water(1,1505,1550,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt');

%% checks general voigt/vhh mex
[d3,w3] = run8(3,1050,1080,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt');

%% checks ckd mex for O2
[d7,w7] = run8(7,1500,1525,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt');

