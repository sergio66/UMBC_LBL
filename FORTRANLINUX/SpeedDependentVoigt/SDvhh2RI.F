      subroutine SDvhh2RI(wr,wi,v,v0,Temp,m,brd,n_in)

c References
c
c Boone, C. D., K. A. Walker, and P. F. Bernath, Speed-dependent Voigt
c profile for water vapor in infrared remote sensing
c applications. JQSRT, 105, 525-532, 2007.
c
c Boone, C. D., K. A. Walker, and P. F. Bernath (2011): An efficient
c analytical approach for calculating line mixing in atmospheric remote
c sensing applications. JQSRT, 112, 980-989.
c
c Humlicek, J., Optimized computation of the Voigt and complex
c probability functions. JQSRT, 27(4), 437-444, 1982.
c
c code from
c Efficient computation of some speed-dependent isolated line profiles
c H. Tran n, N.H. Ngo, J.-M. Hartmann JQSRT 129 (2013) 199–203
c code at http://dx.doi.org/10.1016/j.jqsrt.2013.06.015.

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c subroutine voigt1(wr,wi,0,Temp,m,brd,n_in)
c gets the real and imag parts of the fcn
c v    = frequency array
c v0   = center freq
c T    = temperature
c m    = molecular mass (amu)
c brd  = broadening

      include '/home/sergio/SPECTRA/FORTRANLINUX/max.inc' 

      real*8 wr(MaxLen),wi(MaxLen),v(MaxLen),v0,Temp,m,brd
      integer n_in

c local variables
       real*8 factor,c2
       integer ii
       complex U,T,VTP(MaxLen),VTM(MaxLen)

c local variables

      real*8 k,c_light,amu,mass,r2,alpha_doppler,g0
      real*8 repwid,Y,X,S1V,S2V,fact,c2_0

      IF (n_in .gt. MaxLen) THEN
        print *,'in vhh2RI.f, n_in .gt. MaxLen',n_in,MaxLen
        STOP
        END IF

c  SORT THE (X,Y) PAIRS INTO THE 4 REGIONS OF THE HUMLIcEK 
c  EXPRESSIONS OF THE VOIGT LINE PROFILE.
c

c do the doppler widths first 
      k=1.380658e-23
      c_light=2.99792458e8         !ms-1
      amu=1.6605402e-27            !nucleon mass/kg
      mass=m                       !change to kg

c alpha_doppler=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light)
c r2=2*log(2)*k/amu
      r2=11526.218
      alpha_doppler=v0/c_light*sqrt(r2*Temp/mass)
      repwid=0.8325546/alpha_doppler

c do the g0 factor 
c  SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.5641895 
c g0=sqrt(log(2)/pi)/alpha_doppler
      g0= 0.83255461 * 0.5641895 / alpha_doppler

c------------------------------  do the v-vi part ---------------------------
      CALL loop_pCqSDHC(v0,alpha_doppler/2,Gam0,Gam2,Shift0,Shift2,anuVC,eta,
     &    n_in,v,vtm)
c------------------------------  do the v+vi part ---------------------------
      CALL loop_pCqSDHC_plus(v0,alpha_doppler/2,Gam0,Gam2,Shift0,Shift2,anuVC,eta,
     &    n_in,v,vtp)
c --------------- now do the VHH part --------------------------------
c from the voigt subroutine we need some adjustment factors
cC
cC  SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.5641895
cC       REPWID = 0.83255461/DOPWID
c       H0 = REPWID*0.5641895*STRPAR   ---> STRPAR is multiplied outside
c       Y = WIDPAR*REPWID
cC
cC  CALCULATE THE COMPLEX PROBABILITY FUNCTION
cC
c       CALL VOIVEC(NUM,NL,NH,NSTP,FF,WNUM,REPWID,Y,VT)
cC 
cC  COMPUTE ABSORPTION DUE TO VOIGT LINE SHAPE
cC
c       DO 10 IP=NL,NH,NSTP
c         XABS(IP) = H0*REAL(VT(IP))
c 10        CONTINUE       
cC

      c2=1.4387863        !K/ cm-1  from Genln2 manual 
 
      c2_0=0.5*c2*v0/Temp 
      c2_0=v0*tanh(c2_0) 
      fact=repwid*0.5641895/c2_0 
 
      c2=0.5*c2/Temp 
 
      do 270 ii=1,n_in 
        factor=fact*v(ii)*tanh(c2*v(ii)) 
c        if (ii .eq. 1) print *,ii,factor,vtp(ii),vtm(ii)
        vtp(ii)=factor*(vtp(ii)+vtm(ii)) 
        wr(ii)=real(vtp(ii)) 
        wi(ii)=aimag(vtp(ii)) 
 270    continue 

c      print *,alpha_doppler,repwid,fact,c2

      return
      end

c************************************************************************
c************************************************************************
c************************************************************************
