## https://www.mathworks.com/matlabcentral/answers/251028-can-t-find-mexopts-bat-file

## https://www.mathworks.com/help/matlab/matlab_external/create-fortran-source-mex-file.html 
##      subroutine timestwo(y_output, x_input)
##      real*8 x_input, y_output
##      y_output = 2.0 * x_input
##      return
##      end

## or
## copyfile(fullfile(matlabroot,'extern','examples','refbook','timestwo.F'),'.','f')

module load GCCcore/13.2.0

which mex
mex timestwo.F
mex -v timestwo.F
