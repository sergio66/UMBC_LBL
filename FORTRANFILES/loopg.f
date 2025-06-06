      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

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
c LVG         = -1 for lor, 0 for vanhuber, 1 for voigt

      include 'max.inc'

      integer*8 plhs(*),prhs(*)
      integer nlhs,nrhs

      integer mxGetM,mxGetN
      integer*8 mxGetPr,mxCreateFull

      integer*8 outvectp,zisop,mass_isop,brdp,strengthp,centerfreqp
      integer*8 wavenumberp,temprp,numlinesp,sizewavep,lvgp
      real*8 raOutvect(MaxLen),raZiso(MaxLen),rmass_isop(20),rlvgp
      real*8 raBrd(MaxLen),raStrenght(MaxLen),raCenterFreq(MaxLen)
      real*8 raWavenumber(MaxLen),rtemprp,rnumlinesp,rsizewavep

      integer m_in,n_in,k_in,MM,o_in,jj,kk,ii
             
c check for proper number of arguments
      if (nrhs .ne. 10) then
        call mexErrMsgTxt('10 input args required')
        endif
      if (nlhs .ne. 1) then
        call mexErrMsgTxt('1 output arg required')
        endif

c want to check sizes of input vectors
      jj=mxGetM(prhs(1)) 
      kk=mxGetN(prhs(1))
      do ii=1,5
        m_in=mxGetM(prhs(ii)) 
        n_in=mxGetN(prhs(ii))
        if ((m_in .gt. MaxLen) .or. (n_in .gt. MaxLen)) then
          print *,'checking array number ii = ',ii
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
      o_in=max(jj,kk)
c want to set size of input wavevector
      m_in=mxGetM(prhs(6))
      n_in=mxGetN(prhs(6))
      MM=max(m_in,n_in)
c want to set size of input line parameters
      jj=mxGetM(prhs(1))
      kk=mxGetN(prhs(1))
      k_in=max(jj,kk)

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
      call mxCopyPtrToReal8(zisop,raZiso,int(k_in))
      call mxCopyPtrToReal8(mass_isop,rmass_isop,int(o_in))
      call mxCopyPtrToReal8(brdp,raBrd,int(k_in))
      call mxCopyPtrToReal8(strengthp,raStrenght,int(k_in))
      call mxCopyPtrToReal8(centerfreqp,raCenterFreq,int(k_in))
      call mxCopyPtrToReal8(wavenumberp,raWavenumber,int(MM))
      call mxCopyPtrToReal8(temprp,rtemprp,1)
      call mxCopyPtrToReal8(numlinesp,rnumlinesp,1)
      call mxCopyPtrToReal8(sizewavep,rsizewavep,1)
      call mxCopyPtrToReal8(lvgp,rlvgp,1)

c create a matrix for return argument and assign pointers to the 
c output parameters 
      plhs(1) = mxCreateFull(m_in,n_in,0)
      outvectp = mxGetPr(plhs(1))

c   do the actual computations in a subroutine
      call loop(raOutvect,raZiso,rmass_isop,raBrd,raStrenght,
     $         raCenterFreq,raWavenumber,rtemprp,
     $         int(rnumlinesp),int(rsizewavep),int(o_in),int(rlvgp)) 

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raOutvect, outvectp, max(m_in,n_in))

      return
      end


