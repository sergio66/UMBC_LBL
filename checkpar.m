function [z]=checkpar(fstep,ffin,fmed,nbox,fcor,fmax,fmin,xnear,xmed,xfar); 
%function [z]=checkpar(fstep,ffin,fmed,nbox,fcor,fmax,fmin,xnear,xmed,xfar); 
%checks to see that the parameters sent in make sense

% 
%  TYPE     VAR           DESCRIPTION              TYPICAL VALUE 
%---------------------------------------------------------------------- 
%integer   gasID          HITRAN gas ID                   3 
% 
%integer   fmin           minimum freq (cm-1)            605 
%integer   fmax           maximum freq (cm-1)            630 
% 
%real      ffin           fine point spacing (cm-1)      0.0005 
%real      fmed           medium point spacing (cm-1)    0.1 
%real      fcor           coarse point spacing (cm-1)    1.0 
% 
%real      fstep          wide mesh width size (cm-1)      1.0 
%real      xnear          near wing distance(cm-1)         1.0 
%real      xmed           med wing distance(cm-1)          2.0 
%real      xfar           far wing distance(cm-1)          25.0 
% 
%integer   nbox           boxcar sum size (odd integer)    1,5 
% 
%restrictions : 
%(1) xnear <= xmed <= xfar         
%(2) xnear >= fstep 
%(3) xmed/fmed  xnear/ffin  fstep/fmed   fstep/ffin        are integers 
%(4)fstep/(nbox*ffin)        fcor/ffin                     are integers 
%(5)(fmax-fmin)/fstep        (fmax-fmin)/fcor              are integers 
% 
%thus ffin,fmed,fcor are all 10/2^n   n <= 5 

z=1;

if (fmin > fmax)
  z=0;
  fprintf(1,'fmin, fmax = %8.5e %8.5e \n',fmin,fmax);
  error('fmin > fmax')
end

if (ffin > fmed)
  z=0;
  fprintf(1,'ffin, fmed = %8.5e %8.5e \n',ffin,fmed);
  error('ffin > fmed')
end

if (fmed > fcor)
  z=0;
  fprintf(1,'fmed, fcor = %8.5e %8.5e \n',fmed,fcor);
  error('fmed > fcor')
end

if (xfar < xnear)
  z=0;
  fprintf(1,'xfar, xnear = %8.5e %8.5e \n',xfar,xnear)
  error('need xfar >= xnear');
end

if (xfar < xmed)
  z=0;
  fprintf(1,'xfar, xmed = %8.5e %8.5e \n',xfar,xmed)
  error('need xfar >= xmed');
end

if (xmed < xnear)
  z=0;
  fprintf(1,'xmed, xnear = %8.5e %8.5e \n',xmed,xnear)
  error('need xmed >= xnear');
end

if (fstep > xnear)
  z=0;
  fprintf(1,'fstep, xnear = %8.5e %8.5e \n',fstep,xnear)
  error('need fstep <= xnear');
end

if (isint(mod(xmed,fmed)) ~= 1)
  z=0;
  fprintf(1,'xmed, fmed = %8.5e %8.5e \n',xmed,fmed)
  error('isint(xmed/fmed) ~= 1');
end

if (isint(mod(xnear,ffin)) ~= 1)
  z=0;
  fprintf(1,'xnear, ffin = %8.5f  %8.5f \n',xnear, ffin);
  error('isint(xnear/ffin) ~= 1');
end

if (isint(mod(fstep,ffin)) ~= 1)
  z=0;
  fprintf(1,'fstep, ffin = %8.5f  %8.5f \n',fstep, ffin);
  error('isint(fstep/ffin) ~= 1');
end

if (isint(mod(fstep,fmed)) ~= 1)
  z=0;
  fprintf(1,'fstep, fmed = %8.5f  %8.5f \n',fstep, fmed);
  error('isint(fstep/fmed) ~= 1');
end

if (isint(mod(fstep,(nbox*ffin))) ~= 1)
  z=0;
  fprintf(1,'fstep, nbox, ffin = %8.5f  %3i %8.5f \n',fstep, nbox, ffin);
  error('isint(fstep/(nbox*ffin)) ~= 1');
end

if (isint(mod(fcor,ffin)) ~= 1)
  z=0;
  fprintf(1,'fcor, ffin = %8.5f  %8.5f \n',fcor,ffin);
  error('isint(fcor/ffin) ~= 1');
end

if (isint(mod(fmax-fmin,fstep)) ~= 1)
  if (abs(mod(fmax-fmin,fstep)) > eps*2)
    z=0;
    format long e
    fprintf(1,'fmax, fmin, fstep, (fmax-fmin)/fstep = %8.2f %8.2f %8.6f %8.3f \n',...
            fmax,fmin,fstep,(fmax-fmin)/fstep)
    error('isint((fmax-fmin)/fstep) ~= 1');
  end
end

if (isint(mod(fmax-fmin,fcor)) ~= 1)
  if (abs(mod(fmax-fmin,fcor)) > eps*2)
    z=0;
    error('isint((fmax-fmin)/fcor) ~= 1');
  end
end

