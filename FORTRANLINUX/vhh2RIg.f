      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

c fmex5 vhh1.f vhh1g.f FFLAGS='$FFLAGS -u -64 -mips4' LDFLAGS='$LDFLAGS
c -64 -mips4'

      include 'max.inc'

      integer plhs(*),prhs(*)
      integer nlhs,nrhs

      integer mxGetM,mxGetN
      integer mxGetPr,mxCreateFull

      integer zRealp,zImagp,yp,v0p,Tp,mp,brdp
      real*8 raY(MaxLen),raZReal(MaxLen),raZImag(MaxLen)
      real*8 rTp,rv0p,rmp,rbrdp,mx

      integer m_in,n_in,ii,m_out,n_out
             
c check for proper number of arguments
c want to call the functios as z=boxint2(y,nbox)
      if (nrhs .ne. 5) then
        call mexErrMsgTxt('5 input args required')
        endif
      if (nlhs .ne. 2) then
        call mexErrMsgTxt('2 output arg required')
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
       
      m_out=mxGetM(prhs(2)) 
      n_out=mxGetN(prhs(2))

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

c create a matrix for return arguments and assign pointers to the 
c output parameters z = boxint3(y,nbox,zlen)
      plhs(1) = mxCreateFull(m_in,n_in,0)
      zRealp    = mxGetPr(plhs(1))
      plhs(2) = mxCreateFull(m_in,n_in,0)
      zImagp    = mxGetPr(plhs(2))

c   do the actual computations in a subroutine
      call vhh2RI(raZReal,raZImag,raY,rv0p,rTp,rmp,rbrdp,int(mx))

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raZReal, zRealp, int(mx))
      call mxCopyReal8ToPtr(raZImag, zImagp, int(mx))
      return
      end


