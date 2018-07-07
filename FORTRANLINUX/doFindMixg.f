      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

c fmex5 vhh1.f vhh1g.f FFLAGS='$FFLAGS -u -64 -mips4' LDFLAGS='$LDFLAGS
c -64 -mips4'

c      subroutine doFindMix(Y_1st,trans_ampl,W_matrix,freqq,n_in)  
c Y_1st= results array
c trans= frequency array
c W    = center freq
c freqq= temperature
c n_in = no of input points

      include 'max.inc'

      integer plhs(*),prhs(*)
      integer nlhs,nrhs

      integer mxGetM,mxGetN
      integer mxGetPr,mxCreateFull

      integer trp,Wp,frp,zp
      real*8 raZ(MaxPQR),raTr(MaxPQR),raFr(MaxPQR),raaW(MaxPQR,MaxPQR)

      integer m_in,n_in,x_in,y_in,a_in,b_in,nn,bb
             
c check for proper number of arguments
c want to call the functios as z=doFindMix(trans_ampl,W_matrix,freqq)
      if (nrhs .ne. 3) then
        print *,'doFindMixg.f : nrhs = ',nrhs
        call mexErrMsgTxt('3 input args required')
        endif
      if (nlhs .ne. 1) then
        call mexErrMsgTxt('1 output arg required')
        endif

c want to check sizes of input wavevector array "trans"
      m_in=mxGetM(prhs(1)) 
      n_in=mxGetN(prhs(1))
      if ((m_in .ne. 1)  .and.  (n_in .ne. 1)) then
        call mexErrMsgTxt('transampl needs to be (1,ylen) or (ylen,1)')
        endif
      if ((m_in .gt. MaxPQR) .or. (n_in .gt. MaxPQR)) then 
        call mexErrMsgTxt('PQR lines have to be less than MaxPQR') 
        endif                
      nn=max(m_in,n_in)

c want to check sizes of input matrix "W"
      x_in=mxGetM(prhs(2)) 
      y_in=mxGetN(prhs(2))
      if (x_in .ne. y_in) then
        call mexErrMsgTxt('matrix W needs to be square')
        endif
      if ((x_in .gt. MaxPQR) .or. (y_in .gt. MaxPQR)) then 
        call mexErrMsgTxt('PQR lines have to be less than MaxPQR') 
        endif                

c want to check sizes of input wavevector array "freqq"
      a_in=mxGetM(prhs(3)) 
      b_in=mxGetN(prhs(3))
      if ((a_in .ne. 1)  .and.  (b_in .ne. 1)) then
        call mexErrMsgTxt('freqq needs to be (1,ylen) or (ylen,1)')
        endif
      if ((a_in .gt. MaxPQR) .or. (b_in .gt. MaxPQR)) then 
        call mexErrMsgTxt('PQR lines have to be less than MaxPQR') 
        endif                
      bb=max(a_in,b_in)
       
      trp   = mxGetPr(prhs(1))
      Wp    = mxGetPr(prhs(2))
      frp   = mxGetPr(prhs(3))

c copy right hand arguments to local arrays or variables       
c z = boxint3(y,v0,T,m,brd)
      call mxCopyPtrToReal8(trp, raTr, int(nn))
      call mxCopyPtrToReal8(Wp,  raaW, MaxPQR*MaxPQR)
      call mxCopyPtrToReal8(frp, raFr, int(nn))

c create a matrix for return argument and assign pointers to the 
c output parameters z = boxint3(y,nbox,zlen)
      plhs(1) = mxCreateFull(a_in,b_in,0)
      zp    = mxGetPr(plhs(1))

c   do the actual computations in a subroutine
      call doFindMix(raZ,raTr,raaW,raFr,int(nn))

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raZ, zp, int(nn))

      return
      end


