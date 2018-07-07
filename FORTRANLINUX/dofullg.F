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

c the .m call to this is all the arguments below except z=dofull2(...,birn)
c      subroutine dofull(z,v,v0,w_tot,
c     $           temperature,duration,pressure_for,pressure_self,
c     $           arg1,arg2,diagL,ksm,bsm,birn,
c     $           lenf,no_lines) 
c
c z              = output array
c ksm            = go to correct units
c bsm            = band strength mult ~ 1
c v              = input wavevector
c w_tot          = lorentz broadening coeffs
c arg1,arg2      = from the eigenmatrix; related to population densities
c diagL          = from the eigenmatrix; related to population densities
c v0             = line centers
c temperature,duration = for birnbaum
c pfro,pself           = for cousin
c lenf,no_lines  = length(f),length(v0)
c birn           = -1 (cousin),0 (no chi),1 (birnbaum)
c
c      integer lenf,no_lines,birn 
c      real*8 z(lenf),ksm,bsm,f(lenf) 
c      real*8 w_tot(no_lines),v0(no_lines) 
c      complex*8 arg1(no_lines),arg2(no_lines),diagL(no_lines) 

      include 'max.inc'
 
      integer*2 plhs(*),prhs(*)
      integer nlhs,nrhs

      integer mxGetM,mxGetN
      integer mxGetPr,mxGetPi,mxCreateDoubleMatrix

      integer zp,vp,v0p,Tp,tau2p,w_totp,ksmp,bsmp,pforp,pselfp,birnp,
     $  arg1pr,arg2pr,diagLpr,arg1pi,arg2pi,diagLpi
      real*8 raV(MaxLen),raZ(MaxLen),rTp,raV0(MaxPQR),rtau2p,
     $          rpforp,rpselfp,rksmp,rbsmp,rbirnp,ra_Wtot(MaxPQR)
      complex*16 carg1p(MaxPQR),carg2p(MaxPQR),cdiagLp(MaxPQR)

      mwSize k_in,l_in,m_in,n_in,nnin,llin
             
c check for proper number of arguments
c want to call the functios as z=birn(v,v0,wtot,T,tau2)
      if (nrhs .ne. 13) then
        call mexErrMsgTxt('13 input args required')
        endif
      if (nlhs .ne. 1) then
        call mexErrMsgTxt('1 output arg required')
        endif

c want to check sizes of input wavevector array "y"
      m_in=mxGetM(prhs(1)) 
      n_in=mxGetN(prhs(1))
      if ((m_in .gt. MaxLen) .or. (n_in .gt. MaxLen)) then
        call mexErrMsgTxt('array v has to be smaller than MaxLen')
        endif        
      if ((m_in .ne. 1)  .and.  (n_in .ne. 1)) then
        call mexErrMsgTxt('array v needs to be (1,ylen) or (ylen,1)')
        endif
      nnin=max(m_in,n_in)
c want to check sizes of input linecenters "v0"
      k_in=mxGetM(prhs(2)) 
      l_in=mxGetN(prhs(2))
      if ((k_in .gt. MaxPQR) .or. (k_in .gt. MaxPQR)) then
        call mexErrMsgTxt('PQR lines have to be smaller than MaxPQR')
        endif        
      if ((k_in .ne. 1)  .and.  (l_in .ne. 1)) then
        call mexErrMsgTxt('array v0 needs to be (1,ylen) or (ylen,1)')
        endif
      llin=max(k_in,l_in)

      vp    = mxGetPr(prhs(1))
      v0p   = mxGetPr(prhs(2))
      w_totp= mxGetPr(prhs(3))
      Tp    = mxGetPr(prhs(4))
      tau2p = mxGetPr(prhs(5))
      pforp = mxGetPr(prhs(6))
      pselfp= mxGetPr(prhs(7))

      arg1pr = mxGetPr(prhs(8))
      arg2pr = mxGetPr(prhs(9))
      diagLpr= mxGetPr(prhs(10))
      arg1pi = mxGetPi(prhs(8))
      arg2pi = mxGetPi(prhs(9))
      diagLpi= mxGetPi(prhs(10))

      ksmp  = mxGetPr(prhs(11))
      bsmp  = mxGetPr(prhs(12))
      birnp = mxGetPr(prhs(13))

c copy right hand arguments to local arrays or variables        
      call mxCopyPtrToReal8(vp, raV, int(nnin))
      call mxCopyPtrToReal8(v0p, raV0, int(llin))
      call mxCopyPtrToReal8(w_totp, ra_Wtot, int(llin))
      call mxCopyPtrToReal8(Tp, rTp, 1)
      call mxCopyPtrToReal8(tau2p, rtau2p, 1)
      call mxCopyPtrToReal8(pforp,  rpforp,  1)
      call mxCopyPtrToReal8(pselfp, rpselfp, 1)
      call mxCopyPtrToComplex16(arg1pr, arg1pi, carg1p, int(llin))
      call mxCopyPtrToComplex16(arg2pr, arg2pi, carg2p, int(llin))
      call mxCopyPtrToComplex16(diagLpr, diagLpi, cdiagLp, int(llin))
      call mxCopyPtrToReal8(ksmp, rksmp, 1)
      call mxCopyPtrToReal8(bsmp, rbsmp, 1)
      call mxCopyPtrToReal8(birnp, rbirnp, 1)

c create a matrix for return argument and assign pointers to the 
c output parameters z = boxint3(y,nbox,zlen)
      plhs(1) = mxCreateDoubleMatrix(m_in,n_in,0)
      zp    = mxGetPr(plhs(1))

c fullmix4=dofull(f,freqq_shift,w_tot,temperature,duration,pressure_for,... 
c                 pressure_self,arg1,arg2,diag(L),K_scale_mixing,bsm,birn); 

      if (abs(rbirnp + 1.0) .le. 1e-3) then
        call mexErrMsgTxt('full mix, cousin incompatible in dofullg.f')
        endif

c   do the actual computations in a subroutine
      call dofull(raZ,raV,raV0,ra_Wtot,rTp,rtau2p,rpforp,rpselfp,
     $              carg1p,carg2p,cdiagLp,rksmp,rbsmp,rbirnp,
     $              int(nnin),int(llin))

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raZ, zp, int(nnin))

      return
      end


