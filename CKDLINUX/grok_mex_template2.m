
% Quick version - known modules
mex_setup Fortran

mod_dir = fullfile(pwd, 'modules');
src_dir = fullfile(pwd, 'src');

% Compile modules first
mex('-c', '-I', mod_dir, '-outdir', mod_dir, fullfile(src_dir, 'mod_precision.f90'));
mex('-c', '-I', mod_dir, '-outdir', mod_dir, fullfile(src_dir, 'mod_constants.f90'));
mex('-c', '-I', mod_dir, '-outdir', mod_dir, fullfile(src_dir, 'mod_utils.f90'));

% Then the real mex file(s)
mex('-I', mod_dir, ...
    fullfile(src_dir, 'mex_main_driver.f90'), ...
    fullfile(mod_dir, '*.o'), ...
    '-outdir', pwd);
