c http://www.mathworks.com/matlabcentral/fileexchange/25934-fortran-95-interface-to-matlab-api-with-extras
c change      *.f                  to      *.F
c change      integer m_in,n_in    to      mwSize m_in,n_in
c change      integer nlhs,nrhs    to      integer*2 nlhs,nrhs
c change      mxCreateFull         to      mxCreateDoubleMatrix
c include     fintrf.h

#include "fintrf.h"

      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

      include 'max.inc'

      integer plhs(*),prhs(*)
      integer*2 nlhs,nrhs

      integer mxGetM,mxGetN
      integer mxGetPr,mxCreateDoubleMatrix

      integer zp,yp,nboxp
      real*8 raY(MaxLenJMH),raZ(MaxLenJMH),rzlenp,rnboxp,mx

      integer ii
      mwSize m_in,n_in,len_in
             
c check for proper number of arguments
c want to call the functios as z=boxint2(y,nbox)
      if (nrhs .ne. 2) then
        print *,'nlhs,nrhs = ',nlhs,nrhs
        call mexErrMsgTxt('2 input args required')
        endif
      if (nlhs .ne. 1) then
        print *,'nlhs,nrhs = ',nlhs,nrhs
        call mexErrMsgTxt('1 output arg required')
        endif

c want to check sizes of input array "y"
      m_in=mxGetM(prhs(1)) 
      n_in=mxGetN(prhs(1))
      if ((m_in .gt. MaxLenJMH) .or. (n_in .gt. MaxLenJMH)) then
        print *,'boxin2g.F = ',m_in,n_in,MaxLen
        call mexErrMsgTxt('array size has to be smaller than MaxLenJMH')
        endif        
      if ((m_in .ne. 1)  .and.  (n_in .ne. 1)) then
        call mexErrMsgTxt('input param y needs to be 1d array')
        endif
       
      yp    = mxGetPr(prhs(1))
      nboxp = mxGetPr(prhs(2))

c copy right hand arguments to local arrays or variables       
c z = boxint3(y,nbox)
c note that in reality nbox and zlen are integers
      call mxCopyPtrToReal8(yp, raY, int(max(n_in,m_in)))
      call mxCopyPtrToReal8(nboxp, rnboxp, 1)
      mx=max(m_in,n_in)
      rzlenp=mx/rnboxp
      ii=int(rzlenp)
      
      if (abs(ii-rzlenp) .ge. 1e-5) then
        print *,'boxing2g : ii,rzlenp',ii,rzlenp
        call mexErrMsgTxt('need len(inputarray)/boxp = integer')
        endif

c create a matrix for return argument and assign pointers to the 
c output parameters z = boxint3(y,nbox,zlen)
      if (m_in .eq. 1) then
        len_in = int(n_in/rnboxp)
        plhs(1) = mxCreateDoubleMatrix(m_in,len_in,0)
      elseif (n_in .eq. 1) then
        len_in = int(m_in/rnboxp)
        plhs(1) = mxCreateDoubleMatrix(len_in,n_in,0)
      else
        call mexErrMsgTxt('need inputarray to be 1d')
        endif
      zp    = mxGetPr(plhs(1))

c   do the actual computations in a subroutine
      call boxint2_jmh(raZ,raY,int(rnboxp),int(rzlenp))

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raZ, zp, int(mx/rnboxp))

      return
      end


