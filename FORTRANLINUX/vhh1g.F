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
      include 'max.inc'

      integer plhs(*),prhs(*)
      integer*2 nlhs,nrhs

      integer mxGetM,mxGetN
      integer mxGetPr,mxCreateDoubleMatrix

      integer zp,yp,v0p,Tp,mp,brdp
      real*8 raY(MaxLen),raZ(MaxLen),rTp,rv0p,rmp,rbrdp,mx

      mwSize m_in,n_in
             
c check for proper number of arguments
c want to call the functios as z=boxint2(y,nbox)
      if (nrhs .ne. 5) then
        call mexErrMsgTxt('5 input args required')
        endif
      if (nlhs .ne. 1) then
        call mexErrMsgTxt('1 output arg required')
        endif

c want to check sizes of input wavevector array "y"
      m_in=mxGetM(prhs(1)) 
      n_in=mxGetN(prhs(1))
      if ((m_in .gt. MaxLen) .or. (n_in .gt. MaxLen)) then
        call mexErrMsgTxt('array y has to be smaller than MaxLen')
        endif        
      if ((m_in .ne. 1)  .and.  (n_in .ne. 1)) then
        call mexErrMsgTxt('array y needs to be (1,ylen) or (ylen,1)')
        endif
      mx=max(m_in,n_in)
       
      yp    = mxGetPr(prhs(1))
      v0p   = mxGetPr(prhs(2))
      Tp    = mxGetPr(prhs(3))
      mp    = mxGetPr(prhs(4))
      brdp  = mxGetPr(prhs(5))

c copy right hand arguments to local arrays or variables       
c z = boxint3(y,v0,T,m,brd)
      call mxCopyPtrToReal8(yp, raY, int(mx))
      call mxCopyPtrToReal8(v0p, rv0p, 1)
      call mxCopyPtrToReal8(Tp, rTp, 1)
      call mxCopyPtrToReal8(mp, rmp, 1)
      call mxCopyPtrToReal8(brdp, rbrdp, 1)

c create a matrix for return argument and assign pointers to the 
c output parameters z = boxint3(y,nbox,zlen)
      plhs(1) = mxCreateDoubleMatrix(m_in,n_in,0)
      zp    = mxGetPr(plhs(1))

c   do the actual computations in a subroutine
      call vhh1(raZ,raY,rv0p,rTp,rmp,rbrdp,int(mx))

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raZ, zp, int(mx))

      return
      end


