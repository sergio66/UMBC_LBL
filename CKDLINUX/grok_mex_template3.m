% compile_fortran_mex.m
% This script compiles Fortran source files into a MEX-file using the mex command.
% It assumes your main Fortran file (entry point for MEX) is in the current directory,
% and module files (.f90 sources or precompiled .mod) are in another directory.
% Adjust the paths and file names as needed.

% Define paths and files
main_file = 'your_main_mex_file.f90';  % The main Fortran file with the mexFunction or equivalent entry point
module_dir = '../path/to/modules';    % Directory where module source files or .mod files are located
output_mex = 'your_mex_output';       % Name of the output MEX-file (without extension)

% List additional source files if any (e.g., other .f90 files in current dir or elsewhere)
additional_sources = {'other_file1.f90', 'other_file2.f90'};  % Add full paths if not in current dir

% Compile modules first if they are source files (assuming they need to be compiled to .mod)
% If modules are already precompiled .mod files, skip this step and just use -I
module_sources = fullfile(module_dir, {'module1.f90', 'module2.f90'});  % List module source files
for mod_file = module_sources
    mex('-c', mod_file{1});  % Compile each module to object and .mod (outputs .o and .mod)
end

% Now compile the main MEX-file, including the module directory for .mod includes
% and linking any compiled objects
compiled_objects = fullfile(module_dir, {'module1.o', 'module2.o'});  % List compiled object files from modules

% Build the mex command arguments
mex_args = {main_file};
for src = additional_sources
    mex_args = [mex_args, src];
end
for obj = compiled_objects
    mex_args = [mex_args, obj{1}];
end
mex_args = [mex_args, '-I', module_dir, '-output', output_mex];

% Execute the mex compilation
mex(mex_args{:});

disp(['MEX-file ', output_mex, '.', mexext, ' compiled successfully.']);
