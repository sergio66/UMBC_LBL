      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

      include '../FORTRANFILES/max.inc'

      integer*8 plhs(*),prhs(*)
      integer nlhs,nrhs

      integer mxGetM,mxGetN
      integer*8 mxGetPr,mxCreateFull

      integer*8 ip,np,fp,fsp,nlap,tp,pp,ppp,ap,wckdp,wlp,zp
      real*8 raT(kMaxLayer),raP(kMaxLayer),raPP(kMaxLayer)
      real*8 raA(kMaxLayer),raZ(MaxLen),whichlayer
      real*8 idgas,nfreq,raFreq(MaxLen),fstep,nlay,ckd

      integer m_in,n_in,o_in,p_in
             
c check for proper number of arguments
c       SUBROUTINE CALCON23( CON, IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P,  
c     $    PARTP, AMNT, ckd, whichlayer) 
c       the call is con = CALCON**( IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P,  
c     $    PARTP, AMNT, WHICHCKD, WHICHLAYER) 
c  freq= raFreq= array of wavenumbers = Freq(*)
c T,P,PartP,amnt = arrays of layer variables <= 100 == T(*),P(*) etc
c con = raa(kMaxPts) == CON(*)
c WHICH = integer = 00,21,23
c WHICHLAYER = which gas amt etc to use

      if (nrhs .ne. 11) then
        call mexErrMsgTxt('11 input args required')
        endif
      if (nlhs .ne. 1) then
        call mexErrMsgTxt('1 output arg required')
        endif

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
      wlp   = mxGetPr(prhs(11))

c copy right hand arguments to local arrays or variables       
c z = calcon**()
c note that in reality nbox and zlen are integers
      call mxCopyPtrToReal8(ip, idgas, 1)
      call mxCopyPtrToReal8(np, nfreq, 1)
      call mxCopyPtrToReal8(fp, raFreq, int(max(n_in,m_in)))
      call mxCopyPtrToReal8(fsp, fstep, 1)
      call mxCopyPtrToReal8(nlap, nlay, 1)
      call mxCopyPtrToReal8(tp, raT,   int(max(o_in,p_in)))
      call mxCopyPtrToReal8(pp, raP,   int(max(o_in,p_in)))
      call mxCopyPtrToReal8(ppp, raPP, int(max(o_in,p_in)))
      call mxCopyPtrToReal8(ap, raA,   int(max(o_in,p_in)))
      call mxCopyPtrToReal8(wckdp, ckd ,1)
      call mxCopyPtrToReal8(wlp, whichlayer,1)
      
c create a matrix for return argument and assign pointers to the 
c output parameters z = boxint3(y,nbox,zlen)
      plhs(1) = mxCreateFull(m_in,n_in,0)
      zp    = mxGetPr(plhs(1))

c   do the actual computations in a subroutine

cccccccccc all the water stuff is handled by calconwater, calconwaterg
cccccc      if (int(idgas) .eq. 1) then
cccccc        if (int(ckd) .eq. 00) then
ccccccc         print *,'your CKD version = ',int(ckd)
cccccc          call calcon00(raZ,int(idgas),int(nfreq),raFreq,fstep,
cccccc     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))
cccccc        elseif (int(ckd) .eq. 21) then
ccccccc         print *,'your CKD version = ',int(ckd)
cccccc          call calcon21(raZ,int(idgas),int(nfreq),raFreq,fstep,
cccccc     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))
cccccc        elseif (int(ckd) .eq. 23) then
ccccccc         print *,'your CKD version = ',int(ckd)
cccccc          call calcon23(raZ,int(idgas),int(nfreq),raFreq,fstep,
cccccc     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))
cccccc        else
cccccc          print *,'your CKD version = ',int(ckd)
cccccc          print *,'code needs 0,21,23'
cccccc          call mexErrMsgTxt('please retry')
cccccc          end if

c********************************************************
c id=7 ==> O2, id=22 ==> N2
c you might need to go to calcon.f look at the include files there; these
c turn on the appropriate common blocks
c make sure that you don't call the same common block more than once!!!
c********************************************************
c old O2 continuum
c      if (int(idgas) .eq. 7) then
c        call calconOXY(raZ,int(idgas),int(nfreq),raFreq,fstep,
c     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))

c New O2 continuum from LBLRTMv5.10, but only for N2//O2 = 79/21
c      if (int(idgas) .eq. 7) then
c        call calconOXYNEW_LBLRTM(raZ,int(idgas),int(nfreq),raFreq,fstep,
c     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))

c New O2 continuum from LBLRTMv5.10 with arb mix ratio
      if (int(idgas) .eq. 7) then
        call calconOXYNEW(raZ,int(idgas),int(nfreq),raFreq,fstep,
     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))

c********************************************************
c old N2 continuum
c      elseif (int(idgas) .eq. 22) then
c        call calconNIT(raZ,int(idgas),int(nfreq),raFreq,fstep,
c     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))

c New N2 continuum from LBLRTMv5.10, but only for N2//O2 = 79/21
c      elseif (int(idgas) .eq. 22) then
c        call calconNITNEW_LBLRTM(raZ,int(idgas),int(nfreq),raFreq,fstep,
c     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))

c New N2 continuum from LBLRTMv5.10 with arb mix ratio
      elseif (int(idgas) .eq. 22) then
        call calconNITNEW(raZ,int(idgas),int(nfreq),raFreq,fstep,
     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))
c********************************************************

      else
        print *,'your GasID = ',int(idgas)
        print *,'code needs 7,22 (oxygen,nitrogen)'
        call mexErrMsgTxt('please retry')
        endif

c copy output which is stored in local array to matrix output
      call mxCopyReal8ToPtr(raZ, zp, max(m_in,n_in))

      return
      end


