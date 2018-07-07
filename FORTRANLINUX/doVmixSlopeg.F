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
c      subroutine doVmix(z,v,v0,w_tot,temperature,duration,pressure_for, 
c     $      pressure_self,mass,strenqt,ymix,kvm,klm,bsm,birnbirn,NIFNIF,
c     $      ratio,slope,endpt,lenf,no_lines) 
c 
c z              = output array
c kvm,klm        = go to correct units for voigt/lorentz
c bsm            = band strength mult ~ 1
c v              = input wavevector
c w_tot          = lorentz broadening coeffs
c arg1,arg2      = from the eigenmatrix; related to population densities
c diagL          = from the eigenmatrix; related to population densities
c v0             = line centers
c temperature,duration = for birnbaum
c pfro,pself           = for cousin
c mass                 = mass
c lenf,no_lines  = length(f),length(v0)
c birn           = -1 (cousin),0 (no chi),1 (birnbaum)
c ratio          = mixing/lorentz
c slope/endpt    = how the ratio changes away from endpts

      include 'max.inc' 

      integer plhs(*),prhs(*)
      integer*2 nlhs,nrhs

      integer mxGetM,mxGetN
      integer*8 mxGetPr,mxCreateDoubleMatrix

      integer*8 zp,vp,v0p,Tp,tau2p,w_totp,klmp,kvmp,bsmp,
     $          durationp,pforp,pselfp,birnp,massp,strenqtp,ymixp,
     $          nifp,ratiop,slopep,endptp
      real*8 raV(MaxLen),raZ(MaxLen),rTp,raV0(MaxPQR),rtau2p,rdurationp,
     $          rpforp,rpselfp,rkvmp,rklmp,rbsmp,rbirnp,ra_Wtot(MaxPQR),
     $          raStrenqt(MaxPQR),raYmix(MaxPQR),rmassp,rratiop,
     $          rslopep,rendptp,rnifp

      mwSize k_in,l_in,m_in,n_in,nnin,llin,ii
             
c check for proper number of arguments
      if (nrhs .ne. 18) then
        call mexErrMsgTxt('18 input args required')
        endif
      if (nlhs .ne. 1) then
        call mexErrMsgTxt('1 output arg required')
        endif

c want to check sizes of input wavevector array "y"
      m_in=mxGetM(prhs(1)) 
      n_in=mxGetN(prhs(1))
      if ((m_in .gt. MaxLen)  .or. (n_in .gt. MaxLen)) then
        call mexErrMsgTxt('array v needs to be smaller than MaxLen')
        endif
      if ((m_in .ne. 1)  .and.  (n_in .ne. 1)) then
        call mexErrMsgTxt('array v needs to be (1,ylen) or (ylen,1)')
        endif
      nnin=max(m_in,n_in)
c want to check sizes of input linecenters "v0"
      k_in=mxGetM(prhs(2)) 
      l_in=mxGetN(prhs(2))
      if ((k_in .gt. MaxPQR)  .or. (l_in .gt. MaxPQR)) then
        call mexErrMsgTxt('array v0 needs to be smaller than MaxPQR')
        endif
      if ((k_in .ne. 1)  .and.  (l_in .ne. 1)) then
        call mexErrMsgTxt('array v0 needs to be (1,ylen) or (ylen,1)')
        endif
      llin=max(k_in,l_in)

      vp        =  mxGetPr(prhs(1))
      v0p       =  mxGetPr(prhs(2))
      w_totp    =  mxGetPr(prhs(3))
      Tp        =  mxGetPr(prhs(4))
      tau2p     =  mxGetPr(prhs(5))
      pforp     =  mxGetPr(prhs(6))
      pselfp    =  mxGetPr(prhs(7))
      massp     =  mxGetPr(prhs(8))
      strenqtp  =  mxGetPr(prhs(9))
      ymixp     =  mxGetPr(prhs(10))
      kvmp      =  mxGetPr(prhs(11))
      klmp      =  mxGetPr(prhs(12))
      bsmp      =  mxGetPr(prhs(13))
      birnp     =  mxGetPr(prhs(14))
      nifp      =  mxGetPr(prhs(15))
      ratiop    =  mxGetPr(prhs(16))
      slopep    =  mxGetPr(prhs(17))
      endptp    =  mxGetPr(prhs(18))

c copy right hand arguments to local arrays or variables        
      call mxCopyPtrToReal8(vp, raV, int(nnin))
      call mxCopyPtrToReal8(v0p, raV0, int(llin))
      call mxCopyPtrToReal8(w_totp, ra_Wtot, int(llin))
      call mxCopyPtrToReal8(Tp, rTp, 1)
      call mxCopyPtrToReal8(tau2p, rtau2p, 1)
      call mxCopyPtrToReal8(pforp,  rpforp,  1)
      call mxCopyPtrToReal8(pselfp, rpselfp, 1)
      call mxCopyPtrToReal8(massp, rmassp, 1)
      call mxCopyPtrToReal8(strenqtp, raStrenqt, int(llin))
      call mxCopyPtrToReal8(ymixp, raYmix, int(llin))
      call mxCopyPtrToReal8(kvmp, rkvmp, 1)
      call mxCopyPtrToReal8(klmp, rklmp, 1)
      call mxCopyPtrToReal8(bsmp, rbsmp, 1)
      call mxCopyPtrToReal8(birnp, rbirnp, 1)
      call mxCopyPtrToReal8(nifp, rnifp, 1)
      call mxCopyPtrToReal8(ratiop, rratiop, 1)
      call mxCopyPtrToReal8(slopep, rslopep, 1)
      call mxCopyPtrToReal8(endptp, rendptp, 1)

c create a matrix for return argument and assign pointers to the 
c output parameters z = boxint3(y,nbox,zlen)
      plhs(1) = mxCreateDoubleMatrix(m_in,n_in,0)
      zp    = mxGetPr(plhs(1))

c   do the actual computations in a subroutine
      call doVmixSlope(raZ,raV,raV0,ra_Wtot,rTp,rtau2p,
     $    rpforp,rpselfp,rmassp,raStrenqt,raYmix,rkvmp,rklmp,
     $    rbsmp,rbirnp,rnifp,rratiop,rslopep,rendptp,
     $    int(nnin),int(llin))

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raZ, zp, int(nnin))

      return
      end


