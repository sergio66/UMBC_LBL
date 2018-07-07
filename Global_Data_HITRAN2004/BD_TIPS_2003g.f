      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

c /usr/local/matlab/bin/mex  BD_TIPS_2003.f BD_TIPS_2003g.f FFLAGS='$FFLAGS'  LDFLAGS='$LDFLAGS' FLIBS='$FLIBS'

c      Subroutine BD_TIPS_2003(
c     I MOL,		! HITRAN molecule number
c     I Temp,		! temperature in K
c     I ISO,		! isotopomer index
c     O gi,		! state independent degeneracy factor
c     O QT)		! total internal partition sum
c      implicit DOUBLE PRECISION (a-h,o-z)

c call this from MATLAB as qt = q2003(mol,temp,iso,gi)

      include 'max.inc'

      integer plhs(*),prhs(*)
      integer nlhs,nrhs

      integer mxGetM,mxGetN
      integer mxGetPr,mxCreateFull

      integer molp, tempp, isop, gip,  rz
      real*8  rmolp,rtempp,risop,rgip, rzp

      integer m_in,n_in
             
c check for proper number of arguments
c call this from MATLAB as qt = q2003(mol,temp,iso,gi)
      if (nrhs .ne. 4) then
        call mexErrMsgTxt('4 input args required')
        endif
      if (nlhs .ne. 1) then
        call mexErrMsgTxt('1 output arg required')
        endif

      molp  = mxGetPr(prhs(1))
      tempp = mxGetPr(prhs(2))
      isop  = mxGetPr(prhs(3))
      gip   = mxGetPr(prhs(4))

c copy right hand arguments to local arrays or variables       
      call mxCopyPtrToReal8(molp,  rmolp,  1)
      call mxCopyPtrToReal8(tempp, rtempp, 1)
      call mxCopyPtrToReal8(isop,  risop,  1)
      call mxCopyPtrToReal8(gip,   rgip,   1)

c create a matrix for return argument and assign pointers to the 
c output parameters z
      plhs(1) = mxCreateFull(1,1,0)
      zp      = mxGetPr(plhs(1))

c   do the actual computations in a subroutine
c call this from MATLAB as qt = q2003(mol,temp,iso,gi)
      call BD_TIPS_2003(rmolp,rtempp,risop,rgip,rzp)

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(rz, rzp, 1)

      return
      end


