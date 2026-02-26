%{
mex for Matlab : write a script (not a makefile) to compile fortran
files using the mex compiler, but the modules are in another directory
%}

% compile_fortran_with_modules.m
% Compile Fortran mex files when modules (.mod + .o) are in another directory
% Sergio - Feb 2025 version

clear all;
clc;

%% ─── Configuration ────────────────────────────────────────────────────────
fortran_files = { ...
    'my_routine1.f90', ...
    'my_routine2.f90', ...
    'fast_math_stuff.f90' ...
};

module_dir    = fullfile(pwd, 'modules');          % where .mod and .o files go
source_dir    = fullfile(pwd, 'src');              % where .f90 files live
output_dir    = pwd;                               % where .mexw64 / .mexa64 will be created

% Common extra flags you might want/need
extra_mex_flags = [ ...
    '-O3 ' ...                 % or -Og, -Ofast, etc.
    '-g ' ...                  % keep if you want debugging symbols
    '-I"' module_dir '" ' ...  % tell compiler where to find .mod files
    ];

% If you need to link extra libraries
extra_link_flags = '';   % example: '-lblas -llapack'

%% ─── Preparation ──────────────────────────────────────────────────────────
if ~exist(module_dir, 'dir')
    mkdir(module_dir);
end

% Make sure we're looking at the right compiler
mex -setup Fortran

fprintf('Using Fortran compiler: %s\n', mex.getCompilerConfigurations('Fortran','Selected').Name);

%% ─── Compile modules first ────────────────────────────────────────────────
fprintf('\n=== Compiling modules =====================================\n');

% Find all module files (files containing "module" as first non-comment word)
% This is a simple heuristic — adjust if needed
module_files = dir(fullfile(source_dir, '*.f90'));
module_list = {};
for i = 1:length(module_files)
    fname = fullfile(source_dir, module_files(i).name);
    fid = fopen(fname, 'r');
    if fid == -1, continue; end
    
    first_meaningful_line = '';
    while isempty(first_meaningful_line) && ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line) || startsWith(line, '!'), continue; end
        first_meaningful_line = line;
    end
    fclose(fid);
    
    if startsWith(lower(first_meaningful_line), 'module') && ...
       ~startsWith(lower(first_meaningful_line), 'module procedure')
        module_list{end+1} = module_files(i).name;
        fprintf('  detected module: %s\n', module_files(i).name);
    end
end

% Compile modules → object files + .mod files in module_dir
for i = 1:length(module_list)
    src = fullfile(source_dir, module_list{i});
    [~, name] = fileparts(module_list{i});
    
    cmd = sprintf('mex %s -c "%s" -outdir "%s" %s', ...
        extra_mex_flags, src, module_dir, extra_link_flags);
    
    fprintf('→ %s\n', cmd);
    status = evalc('eval(cmd);');  % capture output
    
    if status ~= 0
        error('Failed to compile module: %s', module_list{i});
    end
end

%% ─── Compile the actual mex entry points ──────────────────────────────────
fprintf('\n=== Compiling mex gateway routines =========================\n');

for i = 1:length(fortran_files)
    src = fullfile(source_dir, fortran_files{i});
    [~, name] = fileparts(fortran_files{i});
    
    % Important: object files of modules + -I for .mod
    cmd = sprintf(['mex %s ' ...
                   '-I"%s" ' ...
                   '"%s" ' ...
                   '%s/*.o ' ...           % link all module objects
                   '-outdir "%s" ' ...
                   '%s'], ...
        extra_mex_flags, ...
        module_dir, ...
        src, ...
        module_dir, ...
        output_dir, ...
        extra_link_flags);
    
    fprintf('→ %s\n', strrep(cmd, pwd, '.'));  % shorten output
    
    try
        eval(cmd);
        fprintf('   → successfully created  %s.%s\n', name, mexext);
    catch ME
        fprintf(2, 'Compilation failed for %s\n', fortran_files{i});
        fprintf(2, '%s\n', ME.message);
        return
    end
end

fprintf('\nAll done.\n');
fprintf('Modules (.mod + .o) → %s/\n', module_dir);
fprintf('MEX files           → current folder\n');

