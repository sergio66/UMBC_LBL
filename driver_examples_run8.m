%% checks water without basement mex
[d1,w1] = run8water(1,1505,1550,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt');

%% checks general voigt/vhh mex
[d3,w3] = run8(3,1050,1080,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt');

%% checks ckd mex for O2
[d7,w7] = run8(7,1500,1525,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt');


%% check CKD 2.5 vs 3.2 vs 4.3
topts.CKD = 25; [w,d2p5] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',topts);
topts.CKD = 32; [w,d3p2] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',topts);
topts.CKD = 43; [w,d4p3] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',topts);
