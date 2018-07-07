c http://www.mathworks.com/matlabcentral/fileexchange/25934-fortran-95-interface-to-matlab-api-with-extras
c change      *.f                  to      *.F
c include     fintrf.h
c%%%%%%%%%%%%%%%%%%%%%%%%%%%% this sorta works %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c change      integer m_in,n_in    to      mwSize m_in,n_in
c change      integer nlhs,nrhs    to      integer*2 nlhs,nrhs
c change      mxCreateFull         to      mxCreateDoubleMatrix
c%%%%%%%%%%%%%%%%%%%%%%%%%%%% this sorta works %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%% this does work   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c change      mxCreateFull         to      mxCreateDoubleMatrix      mwpointer
c change      mxGetPr              to      mwpointer
c leave       integer nlhs,nrhs????
c change      integer nlhs,nrhs    to      integer*2 nlhs,nrhs
c change      plhs,prhs            to      mwpointer
c change all pointers to mwPointer
c change all sizes    to mwSize
c%%%%%%%%%%%%%%%%%%%%%%%%%%%% this does work   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "fintrf.h"
      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

      implicit none

      include '/home/sergio/SPECTRA/FORTRANFILES/max.inc'
      
      mwpointer plhs(*),prhs(*)
      mwpointer mxGetPr,mxCreateDoubleMatrix

      integer nlhs,nrhs

      mwSize mxGetM,mxGetN,mx(11),nx(11)
      mwSize m_in,n_in,one
      mwSize o_in,p_in

      mwPointer ip,np,fp,fsp,nlap,tp,pp,ppp,ap,wckdp,wlp
      real*8 raT(kMaxLayer),raP(kMaxLayer),raPP(kMaxLayer)
      real*8 raA(kMaxLayer),whichlayer
      real*8 idgas,nfreq,raFreq(MaxLen),fstep,nlay,ckd

      mwPointer zp
      real*8 raZ(MaxLen)

      integer ii,iVers
             
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

c      print *,'nlhs nrhs = ',nlhs,nrhs
      if (nrhs .ne. 11) then
        call mexErrMsgTxt('11 input args required')
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
      wlp   = mxGetPr(prhs(11))

c copy right hand arguments to local arrays or variables       
c z = calcon**()
c note that in reality nbox and zlen are integers
      call mxCopyPtrToReal8(ip, idgas, nx(1))
      call mxCopyPtrToReal8(np, nfreq, nx(2))
      call mxCopyPtrToReal8(fp, raFreq, int(max(mx(3),nx(3))))
      call mxCopyPtrToReal8(fsp, fstep, mx(4))
      call mxCopyPtrToReal8(nlap, nlay, mx(5))
      call mxCopyPtrToReal8(tp, raT,   int(max(o_in,p_in)))
      call mxCopyPtrToReal8(pp, raP,   int(max(o_in,p_in)))
      call mxCopyPtrToReal8(ppp, raPP, int(max(o_in,p_in)))
      call mxCopyPtrToReal8(ap, raA,   int(max(o_in,p_in)))
      call mxCopyPtrToReal8(wckdp, ckd ,1)
      call mxCopyPtrToReal8(wlp, whichlayer,1)
      
c create a matrix for return argument and assign pointers to the 
c output parameters z = boxint3(y,nbox,zlen)
      plhs(1) = mxCreateDoubleMatrix(m_in,n_in,0)
      zp    = mxGetPr(plhs(1))

c   do the actual computations in a subroutine      
c********************************************************
c id=7 ==> O2, id=22 ==> N2
c you might need to go to calcon.f look at the include files there; these
c turn on the appropriate common blocks
c make sure that you don't call the same common block more than once!!!

ccccc
      iVers = 1     !!!! this is LBL    v1.0  (orig)
      ! used for orig database 2000 2002
ccccc
      iVers = 2     !!!! this is LBLRTM v5.10 but only  N2//O2 = 79/21
      !unused
ccccc
      iVers = 3     !!!! this is LBLRTM v5.10 for arb N2, O2
      ! used for HITRAN2UMBCLBL/MAKEIR_ALL_H08 database (july 2010)
ccccc
      iVers = 4     !!!! this is LBLRTM v12.1 for arb N2, O2
ccccc

      iVers = 1     !!!! this is LBL    v1.0  (orig)
      iVers = 3     !!!! this is LBLRTM v5.10 for arb N2, O2
      iVers = 4     !!!! this is LBLRTM v12.1 for arb N2, O2

c********************************************************
      if (int(idgas) .eq. 7) then
c        if (iVers .EQ. 1) THEN
c          ! works!!!! used for 2000, 2002 databases
c          ! old O2 continuum
c          call calconOXY(raZ,int(idgas),int(nfreq),raFreq,fstep,
c     $              int(nlay),raT,raP,raPP,raA,int(whichlayer))
c        elseif (iVers .EQ. 2) THEN         
c          ! New O2 continuum from LBLRTMv5.10, but only for N2//O2 = 79/21
c          call calconOXYNEW_LBLRTM(raZ,int(idgas),
c     $              int(nfreq),raFreq,fstep,
c     $              int(nlay),raT,raP,raPP,raA,int(whichlayer))
c        elseif (iVers .EQ. 3) THEN         
c         ! used for HITRAN2UMBCLBL/MAKEIR_ALL_H08 database (july 2010)
c          !New O2 continuum from LBLRTMv5.10 with arb mix ratio
c          call calconOXYNEW(raZ,int(idgas),int(nfreq),raFreq,fstep,
c     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))
c        elseif (iVers .EQ. 4) THEN         
c          !New O2 continuum from LBLRTMv12.1 with arb mix ratio
          call calconOXYNIT_12_1(raZ,int(idgas),int(nfreq),raFreq,fstep,
     $            int(nlay),raT,raP,raPP,raA,int(whichlayer))
c        end if
c********************************************************
      elseif (int(idgas) .eq. 22) then
c        if (iVers .EQ. 1) THEN
c          ! works!!!! used for 2000, 2002 databases
c          ! old N2 continuum
c          call calconNIT(raZ,int(idgas),int(nfreq),raFreq,fstep,
c     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))
c        elseif (iVers .EQ. 2) THEN
c          ! New N2 continuum from LBLRTMv5.10, but only for N2//O2 = 79/21
c          call calconNITNEW_LBLRTM(raZ,
c     $                  int(idgas),int(nfreq),raFreq,fstep,
c     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))
c        elseif (iVers .EQ. 3) THEN
c         ! used for HITRAN2UMBCLBL/MAKEIR_ALL_H08 database (july 2010)
c          ! New N2 continuum from LBLRTMv5.10 with arb mix ratio
c          call calconNITNEW(raZ,int(idgas),int(nfreq),raFreq,fstep,
c     $                  int(nlay),raT,raP,raPP,raA,int(whichlayer))
c        elseif (iVers .EQ. 4) THEN         
c          !New N2 continuum from LBLRTMv12.1 with arb mix ratio
          print *,'New N2 continuum from LBLRTMv12.1 with arb mix ratio'
          call calconOXYNIT_12_1(raZ,int(idgas),int(nfreq),raFreq,fstep,
     $            int(nlay),raT,raP,raPP,raA,int(whichlayer))
c        end if
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


