      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

      include 'max.inc'

      integer*8 plhs(*),prhs(*)
      integer nlhs,nrhs

      integer mxGetM,mxGetN
      integer*8 mxGetPr,mxCreateFull

      integer*8 zp,yp,nboxp
      real*8 raY(MaxLen),raZ(MaxLen),rzlenp,rnboxp,mx

      integer m_in,n_in,ii
             
c check for proper number of arguments
c want to call the functios as z=boxint2(y,nbox)
      if (nrhs .ne. 2) then
        call mexErrMsgTxt('2 input args required')
        endif
      if (nlhs .ne. 1) then
        call mexErrMsgTxt('1 output arg required')
        endif

c want to check sizes of input array "y"
      m_in=mxGetM(prhs(1)) 
      n_in=mxGetN(prhs(1))
      if ((m_in .gt. MaxLen) .or. (n_in .gt. MaxLen)) then
        call mexErrMsgTxt('array size has to be smaller than MaxLen')
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
        call mexErrMsgTxt('need len(inputarray)/boxp = integer')
        endif

c create a matrix for return argument and assign pointers to the 
c output parameters z = boxint3(y,nbox,zlen)
      if (m_in .eq. 1) then
        plhs(1) = mxCreateFull(m_in,int(n_in/rnboxp),0)
      elseif (n_in .eq. 1) then
        plhs(1) = mxCreateFull(int(m_in/rnboxp),n_in,0)
      else
        call mexErrMsgTxt('need inputarray to be 1d')
        endif
      zp    = mxGetPr(plhs(1))

c   do the actual computations in a subroutine
      call boxint2(raZ,raY,int(rnboxp),int(rzlenp))

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raZ, zp, int(mx/rnboxp))

      return
      end


