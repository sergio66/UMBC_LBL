c http://www.mathworks.com/matlabcentral/fileexchange/25934-fortran-95-interface-to-matlab-api-with-extras
c change      *.f                  to      *.F
c change      integer m_in,n_in    to      mwSize m_in,n_in
c change      integer nlhs,nrhs    to      integer*2 nlhs,nrhs
c change      mxCreateFull         to      mxCreateDoubleMatrix
c include     fintrf.h

#include "fintrf.h"

      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

c fmex5 vhh1.f vhh1g.f FFLAGS='$FFLAGS -u -64 -mips4' LDFLAGS='$LDFLAGS
c -64 -mips4'

c the matlab call is : 
c    outvect=loop(ziso,mass_iso,brd,strength,centerfreq,wavenumber,...
c                 tempr,numlines,sizewave,chichi,linetype,pf,ps);  
c and so the FOTRAN call will be
c    subroutine loop(outvect,ziso,mass_iso,brd,strength,
c        centerfreq,wavenumber,tempr,numlines,sizewave,chichi,linetype,pf,ps); 

c outvect     = results vector                                   (sizewave x 1)
c ziso        = vector of mass isotope identifiers (1,2,3 ..)      (N x 1)
c mass_iso    = vector of isotope masses (eg for CO2 44,45, ....)  (20 x 1)
c brd         = vector of line broadening cm-1                     (N x 1)
c strength    = vector of line strengths                           (N x 1)
c centerfreq  = vector of line centers                             (N x 1)
c wavenumber  = vector over which to compute line shapes         (sizewave x 1)
c tempr       = layer temperature                                  
c numlines    = number of line centers
c sizewave    = number of wavevector points to compute shapes over
c chichi      = +1 if cousin, -1 o/w
c linetype    = +1 if voigt/vanvleck, -1 if lorentz/vanvleck
c pf          = foreign pressure
c ps          = self pressure

      include 'max.inc'

      integer plhs(*),prhs(*)
      integer*2 nlhs,nrhs

      integer mxGetM,mxGetN
      integer mxGetPr,mxCreateDoubleMatrix

      integer outvectp,zisop,mass_isop,brdp,strengthp,centerfreqp
      integer wavenumberp,temprp,numlinesp,sizewavep
      integer chichip,linetypep,pfp,psp
      real*8 raOutVect(MaxLen),raZiso(MaxLen),rmass_isop(20)
      real*8 raBrd(MaxLen),raStrength(MaxLen),raCenterFreq(MaxLen)
      real*8 raWavenumber(MaxLen),rtemprp,rnumlinesp,rsizewavep
      real*8 rchichip,rlinetypep,rpfp,rpsp,mm,oo,kk

      mwSize m_in,n_in
      integer o_in,p_in,k_in,l_in,ii
             
c check for proper number of arguments
      if (nrhs .ne. 13) then
        call mexErrMsgTxt('13 input args required')
        endif
      if (nlhs .ne. 1) then
        call mexErrMsgTxt('1 output arg required')
        endif

c want to check sizes of input vectors
      do ii=1,6
        m_in=mxGetM(prhs(ii)) 
        n_in=mxGetN(prhs(ii))
        if ((m_in .gt. MaxLen) .or. (n_in .gt. MaxLen)) then
          print *,'loopco2g.F: checking array number ii = ',ii,m_in,n_in,MaxLen
          call mexErrMsgTxt('array size has to be smaller than MaxLen')
          endif        
        if ((m_in .ne. 1)  .and.  (n_in .ne. 1)) then
          call mexErrMsgTxt('args 1:6 are arrays')
          endif
        enddo

c want to set size of input mass isotope vector
      o_in=mxGetN(prhs(2))
      p_in=mxGetM(prhs(2))
      oo=max(o_in,p_in)
c want to set size of input wavevector
      m_in=mxGetM(prhs(6))
      n_in=mxGetN(prhs(6))
      mm=max(m_in,n_in)
c want to set size of input line parameters
      k_in=mxGetM(prhs(1))
      l_in=mxGetN(prhs(1))
      kk=max(k_in,l_in)
      do ii=3,5
        if (max(k_in,l_in) .ne. kk) then 
          call mexErrMsgTxt('input line parameters are not same sizes!')
          endif
        end do

c     outvect=loop(ziso,mass_iso,brd,strength,... 
c         centerfreq,wavenumber,tempr,numlines,sizewave);    
      zisop        = mxGetPr(prhs(1))
      mass_isop    = mxGetPr(prhs(2))
      brdp         = mxGetPr(prhs(3))
      strengthp    = mxGetPr(prhs(4))
      centerfreqp  = mxGetPr(prhs(5))
      wavenumberp  = mxGetPr(prhs(6))
      temprp       = mxGetPr(prhs(7))
      numlinesp    = mxGetPr(prhs(8))
      sizewavep    = mxGetPr(prhs(9))
      chichip      = mxGetPr(prhs(10))
      linetypep    = mxGetPr(prhs(11))
      pfp          = mxGetPr(prhs(12))
      psp          = mxGetPr(prhs(13))

c copy right hand arguments to local arrays or variables       
c z = boxint3(y,v0,T,m,brd)
      call mxCopyPtrToReal8(zisop,raZiso,int(kk))
      call mxCopyPtrToReal8(mass_isop,rmass_isop,int(oo))
      call mxCopyPtrToReal8(brdp,raBrd,int(kk))
      call mxCopyPtrToReal8(strengthp,raStrength,int(kk))
      call mxCopyPtrToReal8(centerfreqp,raCenterFreq,int(kk))
      call mxCopyPtrToReal8(wavenumberp,raWavenumber,int(mm))
      call mxCopyPtrToReal8(temprp,rtemprp,1)
      call mxCopyPtrToReal8(numlinesp,rnumlinesp,1)
      call mxCopyPtrToReal8(sizewavep,rsizewavep,1)
      call mxCopyPtrToReal8(chichip,rchichip,1)
      call mxCopyPtrToReal8(linetypep,rlinetypep,1)
      call mxCopyPtrToReal8(pfp,rpfp,1)
      call mxCopyPtrToReal8(psp,rpsp,1)

c create a matrix for return argument and assign pointers to the 
c output parameters 
      plhs(1) = mxCreateDoubleMatrix(m_in,n_in,0)
      outvectp = mxGetPr(plhs(1))

c   do the actual computations in a subroutine
      call loopco2(raOutVect,raZiso,rmass_isop,raBrd,raStrength,
     $                raCenterFreq,raWavenumber,rtemprp,
     $                int(rnumlinesp),int(rsizewavep),int(oo),
     $                int(rchichip),int(rlinetypep),rpfp,rpsp) 

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raOutVect, outvectp, int(mm))

      return
      end


