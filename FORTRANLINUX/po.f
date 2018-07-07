      subroutine voigt2R(wr,wi,v,v0,Temp,m,brd,n_in)

      real*8 wr(MaxLen),v(MaxLen),v0,Temp,m,brd
      integer n_in

      real*8 wi(MaxLen)

c local variables
       real*8 factor,c2
       integer ii
       complex U,T,VTP(MaxLen)

c local variables

      real*8 k,c_light,amu,mass,r2,alpha_doppler,g0
      real*8 repwid,Y,X,S1V,S2V

      IF (n_in .gt. MaxLen) THEN
        print *,'in voigt2R.f, n_in .gt. MaxLen',n_in,MaxLen
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
c define arrays for the new Voigt fcn 
      Y=brd/alpha_doppler*0.83255461

       do 20 II=1,n_in
         X =  (V(II) - V0)*REPWID
         S1V = abs(X) + Y
         S2V = (0.195*abs(X)) - 0.176
         T = cmplx(Y,-X)
c
c  FOR REGION 1 OF HUMLIcEK
c
         if (S1V .ge. 15.0) then
           VTP(II) = T*0.5641896/(0.5+(T*T))
c
c  REGION 2 OF HUMLIcEK
c
         elseif (S1V .ge. 5.5) then
           U = T*T
           VTP(II) = T*(1.410474 + U*.5641896)/(.75 + U*(3.+U))
c
c  REGION 3 OF HUMLIcEK
c
         elseif (Y .ge. S2V) then
           VTP(II) =
     1     (16.4955+T*(20.20933+T*(11.96482+T*(3.778987+
     2     T*.5642236))))/
     3     (16.4955+T*(38.82363+T*(39.27121+
     4     T*(21.69274+T*(6.699398+T)))))
c
c  REGION 4 OF HUMLIcEK
c
         else
           U = T*T
           VTP(II)=cexp(U)-T*(36183.31-U*(3321.9905-
     1     U*(1540.787-U*(219.0313-U*
     2     (35.76683-U*(1.320522-U*.56419))))))/
     3     (32066.6-U*(24322.84-U*
     4     (9022.228-U*(2186.181-U*(364.2191-
     5     U*(61.57037-U*(1.841439-U)))))))
         endif
c
       U = cmplx(1, -1)
c
   20  continue
c
c
      c2=1.4387863        !K/ cm-1  from Genln2 manual

      do 270 ii=1,n_in
        factor=repwid*0.5641895
        vtp(ii)=factor*vtp(ii)
        wr(ii)=real(vtp(ii))
        wi(ii)=aimag(vtp(ii))
 270    continue
        
       return
       end

c************************************************************************
