function outvect = xloop(ziso,mass_iso,brd,strength,... 
                   centerfreq,wavenumber,tempr,numlines,sizewave,LVG);  

% see FORTRANLINUUX/loopg.F

% http://www.mathworks.com/matlabcentral/fileexchange/25934-fortran-95-interface-to-matlab-api-with-extras
% change      *.f                  to      *.F
% include     fintrf.h
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% this sorta works %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change      integer m_in,n_in    to      mwSize m_in,n_in
% change      integer nlhs,nrhs    to      integer*2 nlhs,nrhs
% change      mxCreateFull         to      mxCreateDoubleMatrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% this sorta works %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% this does work   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change      mxCreateFull         to      mxCreateDoubleMatrix      mwpointer
% change      mxGetPr              to      mwpointer
% leave       integer nlhs,nrhs????
% change      integer nlhs,nrhs    to      integer*2 nlhs,nrhs
% change      plhs,prhs            to      mwpointer
% change all pointers to mwPointer
% change all sizes    to mwSize
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% this does work   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fmex5 vhh1.f vhh1g.f FFLAGS='$FFLAGS -u -64 -mips4' LDFLAGS='$LDFLAGS
% -64 -mips4'

% the matlab call is : 
% and so the FOTRAN call will be
%    subroutine loop(outvect,ziso,mass_iso,brd,strength,
%                     centerfreq,wavenumber,tempr,numlines,sizewave,LVG); 

% outvect     = results vector                                   (sizewave x 1)
% ziso        = vector of mass isotope identifiers (1,2,3 ..)      (N x 1)
% mass_iso    = vector of isotope masses (eg for CO2 44,45, ....)  (20 x 1)
% brd         = vector of line broadening cm-1                     (N x 1)
% strength    = vector of line strengths                           (N x 1)
% centerfreq  = vector of line centers                             (N x 1)
% wavenumber    = vector over which to compute line shapes         (sizewave x 1)
% tempr       = layer temperature                                  (1 x 1)
% numlines    = number of line centers
% sizewave    = number of wavevector points to compute shapes over
% LVG         = -2 for w/o basement, -1 for lor, 0 for vanhuber, 1 for voigt, 2 for SpeedDependentVoigt

disp('calling loop')
which loop
outvect = loop(ziso,mass_iso,brd,strength,centerfreq,wavenumber,tempr,numlines,sizewave,LVG);
%whos ziso mass_iso brd strength centerfreq wavenumber tempr numlines sizewave LVG outvect
fprintf(1,'isnan(out) isinf(out) = %6i %6i \n',[length(find(isnan(outvect))) length(find(isinf(outvect)))])
disp('in xloop outvect f77 : pause'); pause
disp(' ')

disp('calling x2loop')
which x2loop
outvect2 = x2loop(ziso,mass_iso,brd,strength,centerfreq,wavenumber,tempr,numlines,sizewave,LVG);
%whos ziso mass_iso brd strength centerfreq wavenumber tempr numlines sizewave LVG outvect
fprintf(1,'isnan(out) isinf(out) = %6i %6i \n',[length(find(isnan(outvect))) length(find(isinf(outvect)))])
disp('in xloop outvectB matlab : pause'); pause
disp(' ')

semilogy(wavenumber,outvect,'b.-',wavenumber,outvect2,'r'); grid
disp('ret'); pause;

%{
wah = find(ziso == 1 & centerfreq >= min(wavenumber) & centerfreq <= max(wavenumber));
if length(wah) == 0
  disp('bah 1')
  wah = find(centerfreq >= min(wavenumber) & centerfreq <= max(wavenumber));
end
if length(wah) == 0
  disp('bah 2')
  wah = find(ziso == 1);
end  

bwah = strength(wah);
mwah = find(bwah == max(bwah),1);
gah = wah(mwah);
outvect1 = loop(ziso(gah),mass_iso(ziso(gah)),brd(gah),strength(gah),centerfreq(gah),wavenumber,tempr,1,sizewave,LVG);
junk = [ziso(gah) mass_iso(ziso(gah)) brd(gah) strength(gah) centerfreq(gah)];
fprintf(1,'iso mass brd strength v0 = %3i %8.6f %8.6e %8.6e %8.6f \n',junk);
semilogy(wavenumber,outvect,'b.-',wavenumber,outvect1,'r'); grid
disp(' ')

if length(strength) < 100
  outvect2 = zeros(size(outvect));
  for gah = 1 : -1 : length(strength)
    outvect1 = loop(ziso(gah),mass_iso(ziso(gah)),brd(gah),strength(gah),centerfreq(gah),wavenumber,tempr,1,sizewave,LVG);
    junk = [ziso(gah) mass_iso(ziso(gah)) brd(gah) strength(gah) centerfreq(gah)];
    fprintf(1,'iso mass brd strength v0 = %3i %8.6f %8.6e %8.6e %8.6f \n',junk);
    outvect2 = outvect2 + outvect1;
    semilogy(wavenumber,outvect,'b.-',wavenumber,outvect2,'k',wavenumber,outvect1,'r');
      grid; title([num2str(gah) '  ' num2str(centerfreq(gah))])
    junk = [length(find(isnan(outvect1))) length(find(isinf(outvect1)))];
    fprintf(1,'  line %3i num(isnan) num(isinf) %3i %3i \n',gah,junk)
    pause(0.1)
  end
end  
%}