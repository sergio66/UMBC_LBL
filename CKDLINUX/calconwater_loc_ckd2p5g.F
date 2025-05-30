c http://www.mathworks.com/matlabcentral/fileexchange/25934-fortran-95-interface-to-matlab-api-with-extras
c change      *.f                  to      *.F
c change      integer m_in,n_in    to      mwSize m_in,n_in
c change      integer nlhs,nrhs    to      integer*2 nlhs,nrhs
c change      mxCreateFull         to      mxCreateDoubleMatrix
c include     fintrf.h

#include "fintrf.h"

      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

      include '/home/sergio/SPECTRA/FORTRANFILES/max.inc'

      mwpointer plhs(*),prhs(*)
      mwpointer mxGetPr,mxCreateDoubleMatrix

      integer nlhs,nrhs

      mwsize mxGetM,mxGetN,mx(13),nx(13)
      mwSize m_in,n_in
      mwsize o_in,p_in
             
      mwPointer ip,np,fp,fsp,nlap,tp,pp,ppp,ap,wckdp,wlp,zp
      mwPointer selfp,forp

      real*8 raT(kMaxLayer),raP(kMaxLayer),raPP(kMaxLayer)
      real*8 raA(kMaxLayer),raZ(MaxLen),whichlayer
      real*8 idgas,nfreq,raFreq(MaxLen),fstep,nlay,ckd
      real*8 self,for

c check for proper number of arguments
c       SUBROUTINE CALCON23( CON, IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P,  
c     $    PARTP, AMNT, ckd, selfmult,formult, WHICHLAYER) 
c       the call is con = CALCON**( IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P,  
c     $    PARTP, AMNT, WHICHCKD, selfmult,formult, WHICHLAYER) 
c  freq= raFreq= array of wavenumbers = Freq(*)
c T,P,PartP,amnt = arrays of layer variables <= 100 == T(*),P(*) etc
c con = raa(kMaxPts) == CON(*)
c WHICHCKD = integer = 00,21,23
c WHICHLAYER = which gas amt etc to use
c selfmult = parameter to multiply self part with
c formult  = parameter to multiply for part with
      if (nrhs .ne. 13) then
        call mexErrMsgTxt('13 input args required')
        endif
      if (nlhs .ne. 1) then
        call mexErrMsgTxt('1 output arg required')
        endif

      do ii = 1,nrhs
        o_in=mxGetM(prhs(ii))
        p_in=mxGetN(prhs(ii))
        mx(ii) = o_in
        nx(ii) = p_in
        end do

c want to check sizes of input array "freq"
      m_in=mxGetM(prhs(3)) 
      n_in=mxGetN(prhs(3))
      if ((m_in .gt. MaxLen) .or. (n_in .gt. MaxLen)) then
        call mexErrMsgTxt('wavevector has to be smaller than MaxLen')
        endif        
      if ((m_in .ne. 1)  .and.  (n_in .ne. 1)) then
        call mexErrMsgTxt('input wavevector needs to be 1d array')
        endif

c want to check sizes of input array "T,P,PP,A"
      o_in=mxGetM(prhs(6)) 
      p_in=mxGetN(prhs(6))
      if ((o_in .gt. kMaxLayer) .or. (p_in .gt. kMaxLayer)) then
        call mexErrMsgTxt('p,pp,t,a has to be smaller than kMaxLayer')
        endif        
      if ((o_in .ne. 1)  .and.  (p_in .ne. 1)) then
        call mexErrMsgTxt('p,pp,t,a wavevector needs to be 1d array')
        endif

      ip    = mxGetPr(prhs(1))
      np    = mxGetPr(prhs(2))
      fp    = mxGetPr(prhs(3))
      fsp   = mxGetPr(prhs(4))
      nlap  = mxGetPr(prhs(5))
      tp    = mxGetPr(prhs(6))
      pp    = mxGetPr(prhs(7))
      ppp   = mxGetPr(prhs(8))
      ap    = mxGetPr(prhs(9))
      wckdp = mxGetPr(prhs(10))
      selfp = mxGetPr(prhs(11))
      forp  = mxGetPr(prhs(12))
      wlp   = mxGetPr(prhs(13))

c copy right hand arguments to local arrays or variables       
c z = calcon**()
c note that in reality nbox and zlen are integers
      call mxCopyPtrToReal8(ip, idgas, 1)
      call mxCopyPtrToReal8(np, nfreq, 1)
      call mxCopyPtrToReal8(fp, raFreq, int(max(mx(3),nx(3))))
      call mxCopyPtrToReal8(fsp, fstep, 1)
      call mxCopyPtrToReal8(nlap, nlay, 1)
      call mxCopyPtrToReal8(tp, raT,   int(max(o_in,p_in)))
      call mxCopyPtrToReal8(pp, raP,   int(max(o_in,p_in)))
      call mxCopyPtrToReal8(ppp, raPP, int(max(o_in,p_in)))
      call mxCopyPtrToReal8(ap, raA,   int(max(o_in,p_in)))
      call mxCopyPtrToReal8(wckdp, ckd ,1)
      call mxCopyPtrToReal8(selfp, self ,1)
      call mxCopyPtrToReal8(forp, for ,1)
      call mxCopyPtrToReal8(wlp, whichlayer,1)
      
c create a matrix for return argument and assign pointers to the 
c output parameters z = boxint3(y,nbox,zlen)
      plhs(1) = mxCreateDoubleMatrix(m_in,n_in,0)
      zp    = mxGetPr(plhs(1))

c   do the actual computations in a subroutine

c these are MT-CKD 2.5 or mods of that
      if (int(ckd) .eq. 25) then
        print *,'the MT_CKD version = ',int(ckd)
c        print *,int(idgas),int(nfreq),fstep,int(nlay),self,for,
c     $          int(whichlayer)
c        print *,'calling calcon_mtckd_25_loc .....'
        call calcon_mtckd_25_loc(raZ,int(idgas),int(nfreq),raFreq,
     $   fstep,int(nlay),raT,raP,raPP,raA,self,for,int(whichlayer))
      else
        print *,'the CKD version = ',int(ckd)
        print *,'calconwater_loc_ckd2p5 code needs <25>'
        call mexErrMsgTxt('please retry')
        end if

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raZ, zp, max(m_in,n_in))

      return
      end


