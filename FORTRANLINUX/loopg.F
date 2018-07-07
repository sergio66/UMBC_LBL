c http://www.mathworks.com/matlabcentral/fileexchange/25934-fortran-95-interface-to-matlab-api-with-extras
c change      *.f                  to      *.F
c include     fintrf.h
c%%%%%%%%%%%%%%%%%%%%%%%%%%%% this sorta works %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c change      integer m_in,n_in    to      mwSize m_in,n_in
c change      integer nlhs,nrhs    to      integer*2 nlhs,nrhs
c change      mxCreateFull         to      mxCreateDoubleMatrix
c%%%%%%%%%%%%%%%%%%%%%%%%%%%% this sorta works %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%% this does work   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c change      mxCreateFull         to      mxCreateDoubleMatrix      mwpointer
c change      mxGetPr              to      mwpointer
c leave       integer nlhs,nrhs????
c change      integer nlhs,nrhs    to      integer*2 nlhs,nrhs
c change      plhs,prhs            to      mwpointer
c change all pointers to mwPointer
c change all sizes    to mwSize
c%%%%%%%%%%%%%%%%%%%%%%%%%%%% this does work   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c fmex5 vhh1.f vhh1g.f FFLAGS='$FFLAGS -u -64 -mips4' LDFLAGS='$LDFLAGS
c -64 -mips4'

c the matlab call is : 
c    outvect=loop(ziso,mass_iso,brd,strength,... 
c                   centerfreq,wavenumber,tempr,numlines,sizewave,LVG);  
c and so the FOTRAN call will be
c    subroutine loop(outvect,ziso,mass_iso,brd,strength,
c                     centerfreq,wavenumber,tempr,numlines,sizewave,LVG); 

c outvect     = results vector                                   (sizewave x 1)
c ziso        = vector of mass isotope identifiers (1,2,3 ..)      (N x 1)
c mass_iso    = vector of isotope masses (eg for CO2 44,45, ....)  (20 x 1)
c brd         = vector of line broadening cm-1                     (N x 1)
c strength    = vector of line strengths                           (N x 1)
c centerfreq  = vector of line centers                             (N x 1)
c wavenumber    = vector over which to compute line shapes         (sizewave x 1)
c tempr       = layer temperature                                  (1 x 1)
c numlines    = number of line centers
c sizewave    = number of wavevector points to compute shapes over
c LVG         = -2 for w/o basement, -1 for lor, 0 for vanhuber, 1 for voigt, 2 for SpeedDependentVoigt

#include "fintrf.h"
      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

      implicit none

      include '/home/sergio/SPECTRA/FORTRANFILES/max.inc'
c      include 'max.inc'

      mwpointer plhs(*),prhs(*)
      mwpointer mxGetPr,mxCreateDoubleMatrix

      integer*2 nlhs,nrhs   

      mwSize mxGetM,mxGetN,mx(10),nx(10)
      mwSize m_in,n_in
      mwSize k_in,o_in

      mwPointer outvectp
      real*8 raOutvect(MaxLen)

      mwPointer zisop,mass_isop,brdp,strengthp,centerfreqp
      mwPointer wavenumberp,temprp,numlinesp,sizewavep,lvgp
      real*8 raWavenumber(MaxLen),rtempr,rnumlines,rsizewave,rlvg
      real*8 raZiso(kMaxBandLines),raMass_iso(kMaxIsotopes)
      real*8 raBrd(kMaxBandLines),raStrenght(kMaxBandLines)
      real*8 raCenterFreq(kMaxBandLines)

      integer MM,jj,kk,ii

c      print *,'nrhs,nlhs = ',nrhs,nlhs                   
c       print *,'boooooooooooooooooooooooo 0',nrhs

      do ii = 1,nrhs
        k_in=mxGetM(prhs(ii))
        o_in=mxGetN(prhs(ii))
        mx(ii) = k_in
        nx(ii) = o_in
c        print *,'pointer lenghts',ii,mx(ii),nx(ii)
        end do

c       print *,'boooooooooooooooooooooooo 1'

c check for proper number of arguments
      if (nrhs .ne. 10) then
        print *,'nrhs,nlhs = ',nrhs,nlhs
        call mexErrMsgTxt('10 input args required')
        endif
      if (nlhs .ne. 1) then
        print *,'nrhs,nlhs = ',nrhs,nlhs
        call mexErrMsgTxt('1 output arg required')
        endif

c        print *,'nrhs,nlhs = ',nrhs,nlhs

c want to check sizes of input vectors
      jj=mxGetM(prhs(1)) 
      kk=mxGetN(prhs(1))
      do ii=1,5
        m_in=mxGetM(prhs(ii)) 
        n_in=mxGetN(prhs(ii))
        if ((m_in .gt. MaxLen) .or. (n_in .gt. MaxLen)) then
          print *,'loopg.F : checking array number ii = '
          print *,ii,m_in,n_in,MaxLen
          call mexErrMsgTxt('array size has to be smaller than MaxLen')
          endif
        if ((ii .ne. 2) .and. (ii .ne. 6)) then   
        !mass.ISO, wavenumbers arrays are different
          if (m_in .ne. jj) then
            print *,'checking array number ii = ',ii
            call mexErrMsgTxt('array 1 and this have diff dimensions')         
            endif
          if (n_in .ne. kk) then
            print *,'checking array number ii = ',ii
            call mexErrMsgTxt('array 1 and this have diff dimensions')         
            endif
          endif
        if ((m_in .ne. 1)  .and.  (n_in .ne. 1)) then
          call mexErrMsgTxt('args 1:6 are arrays')
          endif
        enddo

c want to set size of input mass isotope vector
      jj=mxGetN(prhs(2))
      kk=mxGetM(prhs(2))
      o_in=max(int(jj),int(kk))
c want to set size of input wavevector
      m_in=mxGetM(prhs(6))
      n_in=mxGetN(prhs(6))
      MM=max(int(m_in),int(n_in))
c want to set size of input line parameters
      jj=mxGetM(prhs(1))
      kk=mxGetN(prhs(1))
      k_in=max(int(jj),int(kk))

c     outvect=loop(ziso,mass_iso,brd,strength,... 
c         centerfreq,wavenumber,tempr,numlines,sizewave,LVG);    
      zisop        = mxGetPr(prhs(1))
      mass_isop    = mxGetPr(prhs(2))
      brdp         = mxGetPr(prhs(3))
      strengthp    = mxGetPr(prhs(4))
      centerfreqp  = mxGetPr(prhs(5))
      wavenumberp  = mxGetPr(prhs(6))
      temprp       = mxGetPr(prhs(7))
      numlinesp    = mxGetPr(prhs(8))
      sizewavep    = mxGetPr(prhs(9))
      lvgp         = mxGetPr(prhs(10))      

c copy right hand arguments to local arrays or variables       
c z = boxint3(y,v0,T,m,brd)
      call mxCopyPtrToReal8(zisop,raZiso,nx(1))
      call mxCopyPtrToReal8(mass_isop,raMass_iso,nx(2))
      call mxCopyPtrToReal8(brdp,raBrd,nx(3))
      call mxCopyPtrToReal8(strengthp,raStrenght,nx(4))
      call mxCopyPtrToReal8(centerfreqp,raCenterFreq,nx(5))
      call mxCopyPtrToReal8(wavenumberp,raWavenumber,nx(6))
      call mxCopyPtrToReal8(temprp,rtempr,nx(7))
      call mxCopyPtrToReal8(numlinesp,rnumlines,nx(8))
      call mxCopyPtrToReal8(sizewavep,rsizewave,nx(9))
      call mxCopyPtrToReal8(lvgp,rlvg,nx(10))

c create a matrix for return argument and assign pointers to the 
c output parameters 
      plhs(1) = mxCreateDoubleMatrix(m_in,n_in,0)
      outvectp = mxGetPr(plhs(1))

c     do the actual computations in a subroutine
c      do ii = 1,10
c        print *,'gateway',ii,raMass_iso(ii)
c        end do      

      if (abs(rlvg) .GT. 2) then
        print *,'need LVG = -2 -1 0 +1 +2, not ',rlvg
        stop
      end if 
      if (rlvg .EQ. 2) then
        print *,'need LVG = +2 (speed dep voigt) to be implemented'
        stop
      end if 

      call loop(raOutvect,raZiso,raMass_iso,raBrd,raStrenght,
     $         raCenterFreq,raWavenumber,rtempr,
     $         int(rnumlines),int(rsizewave),int(o_in),int(rlvg)) 

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raOutvect, outvectp, max(m_in,n_in))

      return
      end


