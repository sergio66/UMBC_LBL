% compile_fortran_mex_with_modules.m
% Compile a Fortran MEX-file when module sources (.f90) or .mod files are in another folder
% Sergio-friendly version — Feb 2025 style

addpath /home/sergio/git/matlabcode

clear mex   % just in case something is already loaded

% ────────────────────────────────────────────────
%               USER SETTINGS – change these
% ────────────────────────────────────────────────

mexF77 = '/usr/local/matlab/bin/mex';
mexF77 = '/usr/cluster/matlab/r2013a/bin/mex';
mexF77 = '/usr/cluster/matlab/r2016b/bin/mex';  %% this worked way back
mexF77 = '/asl/opt/matlab/R2009b/bin/mex';
mexF77 = '/usr/ebuild/software/MATLAB/2023b/bin/mex';
mexF77 = '/usr/ebuild/software/MATLAB/2021b/bin/mex';
mexF77 = '/usr/ebuild/installs/software/MATLAB/2023b/bin/mex';

main_f90      = 'calconwater_loc_ckd4p3.F90';                       % your mexFunction file (must be in current folder)
module_folder = fullfile('/home/sergio/git/UMBC_LBL/CKDLINUX/MT_CKD_H2O-4.3/build/');    % ←←← most important: where your modules live
module_o      = fullfile('/home/sergio/git/UMBC_LBL/CKDLINUX/MT_CKD_H2O-4.3/build/mt_ckd_h2o_4.3_linux_gnu_dbl.obj/');  %% where the .o lives

output_name   = '../calconwater_loc_ckd4p3';              % name of the final .mexw64 / .mexa64 / .mexmaci64 file

% List module source files (only needed if you want this script to compile them)
% Leave empty {} if the .mod files are already compiled and up-to-date
module_sources = {};
o_sources = {};

% Extra source files in the SAME folder as main_f90 (if any)
extra_sources_here = {'calcon_loc_mtckd_43_wrap.F90','MT_CKD_H2O-4.3/src/cntnm_progr_sergio.f90'};   % example: {'utils.f90', 'helpers.f90'}
extra_sources_here = {'calcon_loc_mtckd_43_wrap.F90'}; % MT_CKD_H2O-4.3/src/cntnm_progr_sergio.f90 is already called by calcon_loc_mtckd_43_wrap.F90

% ────────────────────────────────────────────────
%           Usually no need to change below
% ────────────────────────────────────────────────

% Build include path for .mod files
include_flag = ['-I' module_folder];

% 1. Compile module sources → produce .o and .mod files (if requested)
obj_files = {};

if ~isempty(module_sources)
    fprintf('Compiling modules ...\n');
    for i = 1:length(module_sources)
        f = fullfile(module_folder, module_sources{i});
        fprintf('  %s\n', module_sources{i});
        mex('-c', f, include_flag);   % puts .o and .mod in module_folder
        [~, name] = fileparts(f);
        obj_files{end+1} = fullfile(module_folder, [name '.o']);
    end
    fprintf('\n');
end

% if ~isempty(o_sources)
%     fprintf('Compilingo bjs ...\n');
%     for i = 1:length(o_sources)
%         f = fullfile(module_folder, o_sources{i});
%         fprintf('  %s\n', o_sources{i});
%         mex('-c', f, include_flag);   % puts .o and .mod in module_folder
%         [~, name] = fileparts(f);
%         o_files{end+1} = fullfile(module_folder, [name '.o']);
%     end
%     fprintf('\n');
% end
 
% 2. Collect all pieces for the final link
all_sources = [{main_f90}, extra_sources_here];
all_objects = [obj_files];   
all_o       = [module_o '/mt_ckd_h2o_module.o'];;

% Get flags from nf-config (run this in MATLAB)
[status, nf_libs] = system('nf-config --flibs');
nf_libs = strtrim(nf_libs);  % e.g. "-L/usr/lib/x86_64-linux-gnu -lnetcdff -lnetcdf ..."
[status, nf_inc]  = system('nf-config --fflags');
nf_inc = strtrim(nf_inc);    % for -I includes if needed

% Get flags automatically
[~, nf_libs]  = system('nf-config --flibs');   nf_libs  = strtrim(nf_libs);
[~, nf_inc]   = system('nf-config --fflags');  nf_inc   = strtrim(nf_inc);

nf_inc
nf_libs

% 3. Build mex command
mex_cmd = {mexF77, 'FFLAGS=-fPIC',all_o, ...
          all_sources{:}, all_objects{:}, ...    %% or your main gateway file name
           include_flag, ...
	   nf_libs, ...	                         %% adds -lnetcdff etc.
	   '-I', nf_inc, ...                     %% if needed for includes
           '-output', output_name};              

disp('put -fPIC into ~/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/build/makefile.common and recompile!!! gmake -f make_mt_ckd_h2o linuxGNUdbl')
disp('put -fPIC into ~/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/build/makefile.common and recompile!!! gmake -f make_mt_ckd_h2o linuxGNUdbl')
disp('put -fPIC into ~/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/build/makefile.common and recompile!!! gmake -f make_mt_ckd_h2o linuxGNUdbl')


% Optional: more control over compiler flags (uncomment if needed)
% mex_cmd{end+1} = '-v';                     % verbose
% mex_cmd{end+1} = 'FCFLAGS=-O2 -cpp';       % example gfortran flags
% mex_cmd{end+1} = '-largeArrayDims';        % if you need >2 GB arrays (very old MATLAB)

% 4. Run it
fprintf('Compiling MEX   →  %s.%s\n\n', output_name, mexext);
disp(strjoin(mex_cmd, '  '));
disp('ret to continue'); pause
disp(' ');
boo = strjoin(mex_cmd,' ');
boo = ['!' boo];
%keyboard_nowindow
eval(boo)

%mex(mex_cmd{:});

fprintf('\nDone.\n');
if exist([output_name '.' mexext], 'file')
    fprintf('Success →  <a href="matlab:which %s">%s.%s</a>\n', ...
        output_name, output_name, mexext);
else
    fprintf(2, '→ Compilation failed.\n');
end
