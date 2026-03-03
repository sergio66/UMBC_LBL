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

iTest = 1;
iTest = 2;

if iTest == 1
  prof_file = '/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt'; %% p0 = 1.0 atm. VMR ~ 0.01, T = 300K
  aer_file  = 'mt_ckd_h2o_output_FROM_sergio_test.config1.nc';
  %% in /home/sergio/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/run_example
  %%   ../mt_ckd_h2o_4.3_linux_gnu_dbl < sergio_test1.config
  %%   mv mt_ckd_h2o_output.nc mt_ckd_h2o_output_FROM_sergio_test.config1.nc

elseif iTest == 2
  prof_file = '/home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_two.txt'; %% p0 = 0.8 atm. VMR ~ 0.10, T = 280K
  aer_file  = 'mt_ckd_h2o_output_FROM_sergio_test.config2.nc';
  %% in /home/sergio/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/run_example
  %%   ../mt_ckd_h2o_4.3_linux_gnu_dbl < sergio_test2.config
  %%   mv mt_ckd_h2o_output.nc mt_ckd_h2o_output_FROM_sergio_test.config2.nc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_c1_c2_avog_pref_tref

%% so do this to speed things to get 1 cm-1 output res
topts.ffin = 5.0000e-04 * 400;
topts.CKD = 25; [wx,d2p5] = run8watercontinuum(1,600,2800,prof_file,topts);
topts.CKD = 32; [wx,d3p2] = run8watercontinuum(1,600,2800,prof_file,topts);
topts.CKD = 43; [wx,d4p3] = run8watercontinuum(1,600,2800,prof_file,topts);

%%%%%%%%%%%%%%%%%%%%%%%%%

%% See (3) Initial test, test 3
addpath /home/sergio/git/matlabcode
[x,attrx]   = read_netcdf_lls(['/home/sergio/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/run_example/' aer_file]);
[x0,attrx0] = read_netcdf_lls('/home/sergio/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/run_example/absco-ref_wv-mt-ckd.nc');
ii_x0 = find(x0.wavenumbers >= 600 & x0.wavenumbers <= 2800);
semilogy(x.wavenumbers,x.self_absorption,'k',x0.wavenumbers(ii_x0),x0.self_absco_ref(ii_x0),'m');
semilogy(x.wavenumbers,x.frgn_absorption,'k',x0.wavenumbers(ii_x0),x0.for_absco_ref(ii_x0),'m');
figure(4); plot(x0.wavenumbers,(296/300).^x0.self_texp); xlim([500 3000])  %% I have used 300 K in /home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt
figure(4); semilogy(x0.wavenumbers,x0.self_absco_ref,x0.wavenumbers,x0.for_absco_ref); xlim([500 3000]); legend('Self','Forn')
  ylabel('CS,CF  1/(molecules/cm2)'); xlabel('wavenumber cm-1');

%% /home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt says Q = 5.0e-3 kmol/cm2   and we have <<<< radflag = .TRUE. >>>>>
odx = 5.0e-3 * AVOG * (x.self_absorption + x.frgn_absorption);
modx = interp1(x.wavenumbers,odx,wx);
figure(5); semilogy(wx,d2p5,'gx-',wx,d3p2,'b',wx,d4p3,'r.-',x.wavenumbers,odx,'k'); title('OD');
  legend('CKD2/5','CKD3.2','CKD4.3','AER','location','best')

figure(6); plot(wx,d4p3./d3p2,'b',wx,modx./d3p2,'g',wx,modx./d4p3,'r.-'); legend('d4p3./d3p2','AER OD./d3p2','AER OD ./d4p3','location','best');

disp('main part done!! ret to continue'); pause
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('computing self coeffs')
tsopts.ffin = 5.0000e-04 * 400;
tsopts.divide = +1;
tsopts.selfmult = 1.0;
tsopts.formult  = 0.0;
tsopts.CKD = 25; [wx,sd2p5] = run8watercontinuum(1,600,2800,prof_file,tsopts);
tsopts.CKD = 32; [wx,sd3p2] = run8watercontinuum(1,600,2800,prof_file,tsopts);
tsopts.CKD = 43; [wx,sd4p3] = run8watercontinuum(1,600,2800,prof_file,tsopts);

disp('computing forn coeffs')
tfopts.ffin = 5.0000e-04 * 400;
tfopts.divide = +1;
tfopts.selfmult = 0.0;
tfopts.formult  = 1.0;
tfopts.CKD = 25; [wx,fd2p5] = run8watercontinuum(1,600,2800,prof_file,tfopts);
tfopts.CKD = 32; [wx,fd3p2] = run8watercontinuum(1,600,2800,prof_file,tfopts);
tfopts.CKD = 43; [wx,fd4p3] = run8watercontinuum(1,600,2800,prof_file,tfopts);

mult_aer = 0.00990098; %% see ps/p in /home/sergio/git/UMBC_LBL/IPFILES_EXAMPLE/test_one.txt
junk = load(prof_file);
%% rho_rat = (p_atm/dat%ref_press)*(dat%ref_temp/t_atm)
TT = junk(4); QQ = junk(5); vmr = junk(3)/junk(2); rho_rat = junk(2)/1 * 296/junk(4); rad = x.wavenumbers .* tanh(c2*x.wavenumbers/TT);

mult = 0.0;
mult = 1.0;
figure(1); clf; semilogy(wx,sd2p5,'gx-',wx,sd3p2,'b.-',wx,mult*sd4p3,'rx-',x.wavenumbers,x.self_absorption/(vmr * rho_rat)./rad,'k',x0.wavenumbers(ii_x0),x0.self_absco_ref(ii_x0),'m.-');
  legend('25','32','43 sergio','43 aer is adjusted for VMW','43 aer ref 296 K','location','best'); title('CSelf')
  ylabel('CS  1/(molecules/cm2)'); xlabel('wavenumber cm-1');  
figure(2); semilogy(wx,fd2p5,'gx-',wx,fd3p2,'b.-',wx,mult*fd4p3,'rx-',x.wavenumbers,x.frgn_absorption/(rho_rat * (1-vmr))./rad,'k',x0.wavenumbers(ii_x0),x0.for_absco_ref(ii_x0),'m.-');
  legend('25','32','43 sergio','43 aer is adjusted for (1-vmr)','43 aer ref 296 K','location','best'); title('CForn')
  ylabel('CF  1/(molecules/cm2)'); xlabel('wavenumber cm-1');

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5); semilogy(wx,d2p5,'gx-',wx,d3p2,'b',wx,d4p3,'r.-'); title('OD'); legend('CKD2/5','CKD3.2','CKD4.3','location','best')
  ylabel('OD'); xlabel('wavenumber cm-1');
figure(6); plot(wx,d4p3./d3p2,'b',wx,modx./d3p2,'g',wx,modx./d4p3,'r.-'); legend('d4p3./d3p2','AER OD./d3p2','AER OD ./d4p3','location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iDo = -1;
if iDo > 0
  %% see do_local_lineshape_CKD.m
  tempfreq = AVOG * QQ * wx.*tanh(c2*wx/2/TT)*(296/TT);
  frgratio = fd4p3./interp1(x0.wavenumbers(ii_x0),x0.for_absco_ref(ii_x0),wx,[],'extrap');
  slfratio = sd4p3./interp1(x0.wavenumbers(ii_x0),x0.self_absco_ref(ii_x0),wx,[],'extrap');
  figure(3)
  semilogy(wx,1./(tempfreq/AVOG),'rx-',wx,frgratio,'b',wx,slfratio,'c.-');
    legend('tempfreq in do\_local\_lineshape','CF 4.3 Sergio/CF 4.3 ref','CS 4.3 Sergio/CS 4.3 ref','location','best')  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
