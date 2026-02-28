%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^
disp('see also driver_examples_run8.m')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% check CKD 2.5 vs 3.2 vs 4.3
% 
% %% these are defaults for wmin = 600 cm-1; so output would be at 3e-4 * 5 = 0.0015 cm-1 res
% topts.ffin = 3.0000e-04;
% topts.fmed = 0.06;
% topts.fcor = 0.3;
% %% these are defaults for wmin = 2600 cm-1; so output would be at 5e-4 * 5 = 0.0025 cm-1 res
% topts.ffin = 5.0000e-04;
% topts.fmed = 0.1;
% topts.fcor = 0.5;
% %% so do this to speed things to get 1 cm-1 output res
% topts.ffin = 5.0000e-04 * 400;
% topts.fmed = 0.1 * 400;
% topts.fcor = 0.5 * 400;
% 
% topts = rmfield(topts,'fmed');
% topts = rmfield(topts,'fcor');
% 
% topts.CKD = 25; [wx,d2p5] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',topts);
% topts.CKD = 32; [wx,d3p2] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',topts);
% topts.CKD = 43; [wx,d4p3] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',topts);
% 
% semilogy(wx,d2p5,'gx-',wx,d3p2,'b',wx,d4p3,'r.-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% so do this to speed things to get 1 cm-1 output res
topts.ffin = 5.0000e-04 * 400;
topts.CKD = 25; [wx,d2p5] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',topts);
topts.CKD = 32; [wx,d3p2] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',topts);
topts.CKD = 43; [wx,d4p3] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',topts);

%%%%%%%%%%%%%%%%%%%%%%%%%

%% See (3) Initial test, test 3
addpath /home/sergio/git/matlabcode
[x,attrx]   = read_netcdf_lls('/home/sergio/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/run_example/mt_ckd_h2o_output.nc');
[x0,attrx0] = read_netcdf_lls('/home/sergio/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/run_example/absco-ref_wv-mt-ckd.nc');
ii_x0 = find(x0.wavenumbers >= 600 & x0.wavenumbers <= 2800);
semilogy(x.wavenumbers,x.self_absorption,'k',x0.wavenumbers(ii_x0),x0.self_absco_ref(ii_x0),'m');
semilogy(x.wavenumbers,x.frgn_absorption,'k',x0.wavenumbers(ii_x0),x0.for_absco_ref(ii_x0),'m');
figure(4); plot(x0.wavenumbers,(296/300).^x0.self_texp); xlim([500 3000])  %% I have used 300 K in /home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt
figure(4); semilogy(x0.wavenumbers,x0.self_absco_ref,x0.wavenumbers,x0.for_absco_ref); xlim([500 3000]); legend('Self','Forn')
  ylabel('CS,CF  1/(molecules/cm2)'); xlabel('wavenumber cm-1');
  
tsopts.ffin = 5.0000e-04 * 400;
tsopts.divide = +1;
tsopts.selfmult = 1.0;
tsopts.formult  = 0.0;
tsopts.CKD = 25; [wx,sd2p5] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',tsopts);
tsopts.CKD = 32; [wx,sd3p2] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',tsopts);
tsopts.CKD = 43; [wx,sd4p3] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',tsopts);

tfopts.ffin = 5.0000e-04 * 400;
tfopts.divide = +1;
tfopts.selfmult = 0.0;
tfopts.formult  = 1.0;
tfopts.CKD = 25; [wx,fd2p5] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',tfopts);
tfopts.CKD = 32; [wx,fd3p2] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',tfopts);
tfopts.CKD = 43; [wx,fd4p3] = run8watercontinuum(1,600,2800,'/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt',tfopts);

mult = 0.0;
mult = 1.0;
mult_aer = 0.00990098; %% see ps/p in /home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt
figure(1); clf; semilogy(wx,sd2p5,'gx-',wx,sd3p2,'b.-',wx,mult*sd4p3,'rx-',x.wavenumbers,x.self_absorption/mult_aer,'k',x0.wavenumbers(ii_x0),x0.self_absco_ref(ii_x0),'m.-');
  legend('25','32','43 sergio','43 aer is adjusted for VMW','43 aer ref','location','best'); title('CSelf')
  ylabel('CS  1/(molecules/cm2)'); xlabel('wavenumber cm-1');  
figure(2); semilogy(wx,fd2p5,'gx-',wx,fd3p2,'b.-',wx,mult*fd4p3,'rx-',x.wavenumbers,x.frgn_absorption/(1-mult_aer),'k',x0.wavenumbers(ii_x0),x0.for_absco_ref(ii_x0),'m.-');
  legend('25','32','43 sergio','43 aer is adjusted for (1-vmr)','43 aer ref','location','best'); title('CForn')
  ylabel('CF  1/(molecules/cm2)'); xlabel('wavenumber cm-1');

figure(3); semilogy(wx,d2p5,'gx-',wx,d3p2,'b',wx,d4p3,'r.-'); title('OD'); legend('CKD2/5','CKD3.2','CKD4.3','location','best')
  ylabel('OD'); xlabel('wavenumber cm-1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iDo = -1;
if iDo > 0
  %% see do_local_lineshape_CKD.m
  junk = load('/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt');
  TT = junk(4); QQ = junk(5);
  c2 = 1.4387863
  AVOG = 6.022045E+26
  tempfreq = AVOG * QQ * wx.*tanh(c2*wx/2/TT)*(296/TT);
  frgratio = fd4p3./interp1(x0.wavenumbers(ii_x0),x0.for_absco_ref(ii_x0),wx,[],'extrap');
  slfratio = sd4p3./interp1(x0.wavenumbers(ii_x0),x0.self_absco_ref(ii_x0),wx,[],'extrap');
  figure(3)
  semilogy(wx,1./(tempfreq/AVOG),'rx-',wx,frgratio,'b',wx,slfratio,'c.-');
    legend('tempfreq in do\_local\_lineshape','CF 4.3 Sergio/CF 4.3 ref','CS 4.3 Sergio/CS 4.3 ref','location','best')  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
