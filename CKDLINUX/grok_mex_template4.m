% compile_fortran_mex_with_modules.m
% Compile a Fortran MEX-file when module sources (.f90) or .mod files are in another folder
% Sergio-friendly version — Feb 2025 style

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

main_f90      = 'mex_main.f90';               % your mexFunction file (must be in current folder)
module_folder = fullfile('..', 'modules');    % ←←← most important: where your modules live

output_name   = 'my_fortran_mex';             % name of the final .mexw64 / .mexa64 / .mexmaci64 file

% List module source files (only needed if you want this script to compile them)
% Leave empty {} if the .mod files are already compiled and up-to-date
module_sources = { ...
    'precisions.f90', ...
    'constants.f90', ...
    'my_math_utils.f90', ...
    'grid_tools.f90' ...
};

% Extra source files in the SAME folder as main_f90 (if any)
extra_sources_here = {};   % example: {'utils.f90', 'helpers.f90'}

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

% 2. Collect all pieces for the final link
all_sources = [{main_f90}, extra_sources_here];
all_objects = [obj_files];   % already full paths

% 3. Build mex command
mex_cmd = {'mex', all_sources{:}, all_objects{:}, ...
           include_flag, ...
           '-output', output_name};


% Optional: more control over compiler flags (uncomment if needed)
% mex_cmd{end+1} = '-v';                     % verbose
% mex_cmd{end+1} = 'FCFLAGS=-O2 -cpp';       % example gfortran flags
% mex_cmd{end+1} = '-largeArrayDims';        % if you need >2 GB arrays (very old MATLAB)

% 4. Run it
fprintf('Compiling MEX   →  %s.%s\n\n', output_name, mexext);
disp(strjoin(mex_cmd, '  '));
disp(' ');

mex(mex_cmd{:});

fprintf('\nDone.\n');
if exist([output_name '.' mexext], 'file')
    fprintf('Success →  <a href="matlab:which %s">%s.%s</a>\n', ...
        output_name, output_name, mexext);
else
    fprintf(2, '→ Compilation failed.\n');
end
