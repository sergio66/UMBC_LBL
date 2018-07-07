      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

c fmex5 vhh1.f vhh1g.f FFLAGS='$FFLAGS -u -64 -mips4' LDFLAGS='$LDFLAGS
c -64 -mips4'


c      subroutine cousin(z,v,v0,w_tot,T,pfor,pco2,n_in)
c z    = results array
c v    = frequency array
c v0   = center freq
c T    = temperature
c pfor = foreign press
c pco2 = co2 press
c n_in = no of input points

      include 'max.inc'

      integer plhs(*),prhs(*)
      integer nlhs,nrhs

      integer mxGetM,mxGetN
      integer mxGetPr,mxCreateFull

      integer zp,vp,v0p,Tp,pforp,pco2p,w_totp
      real*8 raV(MaxLen),raZ(MaxLen),rTp,rv0p,rpforp,rpco2p
      real*8 mx,rw_totp

      integer m_in,n_in
             
c check for proper number of arguments
c want to call the functios as z=cousin(v,v0,w_tot,T,pfor,pco2)
      if (nrhs .ne. 6) then
        call mexErrMsgTxt('6 input args required')
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
       
      vp    = mxGetPr(prhs(1))
      v0p   = mxGetPr(prhs(2))
      w_totp= mxGetPr(prhs(3))
      Tp    = mxGetPr(prhs(4))
      pforp = mxGetPr(prhs(5))
      pco2p = mxGetPr(prhs(6))

c copy right hand arguments to local arrays or variables       
c z = boxint3(y,v0,T,m,brd)
      call mxCopyPtrToReal8(vp, raV, int(mx))
      call mxCopyPtrToReal8(v0p, rv0p, 1)
      call mxCopyPtrToReal8(w_totp, rw_totp, 1)
      call mxCopyPtrToReal8(Tp, rTp, 1)
      call mxCopyPtrToReal8(pforp, rpforp, 1)
      call mxCopyPtrToReal8(pco2p, rpco2p, 1)

c create a matrix for return argument and assign pointers to the 
c output parameters z = boxint3(y,nbox,zlen)
      plhs(1) = mxCreateFull(m_in,n_in,0)
      zp    = mxGetPr(plhs(1))

c   do the actual computations in a subroutine
      call cousin(raZ,raV,rv0p,rw_totp,rTp,rpforp,rpco2p, int(mx))

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raZ, zp, int(mx))

      return
      end


