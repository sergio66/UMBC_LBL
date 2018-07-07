      subroutine vhh1(wr,v,v0,Temp,m,brd,n_in)
c subroutine voigt1(wr,v,v0,Temp,m,brd,n_in)
c this is the Van Huber lineshape : using the more accurate dave tobin approx
c v    = frequency array
c v0   = center freq
c T    = temperature
c m    = molecular mass (amu)
c brd  = broadening

      include 'max.inc'

      real*8 wr(MaxLen),v(MaxLen),v0,Temp,m,brd

c local variables
      real*8 res1(MaxLen),res2(MaxLen),factor
      real*8 k,c_light,amu,mass,r2,alpha_doppler,g0
      real*8 X1(MaxLen),T(6),C(6),S(6),d1,d2,d3,d4,d,r,c2
      real*8 wr1(MaxLen),wr2(MaxLen)
      real*8 xr1(MaxLen),xr2(MaxLen),yr1,yr2,y1,y2,y3

      integer ii,jj,n_in,region1(MaxLen),region2(MaxLen),len1,len2
      integer iTruth

      IF (n_in .gt. MaxLen) THEN
        print *,'in vhh1.f, n_in .gt. MaxLen',n_in,MaxLen
        STOP
        END IF


      T(1) = 0.314240376
      T(2) = 0.947788391
      T(3) = 1.59768264
      T(4) = 2.27950708
      T(5) = 3.02063703
      T(6) = 3.8897249

      C(1) = 1.01172805
      C(2) = -0.75197147
      C(3) = 1.2557727E-2
      C(4) = 1.00220082E-2
      C(5) = -2.42068135E-4
      C(6) = 5.00848061E-7

      S(1) = 1.393237
      S(2) = 0.231152406
      S(3) = -0.155351466
      S(4) = 6.21836624E-3
      S(5) = 9.19082986E-5
      S(6) = -6.27525958E-7

c do the doppler widths first 
      k=1.380658e-23
      c_light=2.99792458e8         !ms-1
      amu=1.6605402e-27            !nucleon mass/kg
      mass=m                       !change to kg

c alpha_doppler=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light)
c r2=2*log(2)*k/amu
      r2=11526.218
      alpha_doppler=v0/c_light*sqrt(r2*Temp/mass)

c do the g0 factor 
c  SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.5641895 
c g0=sqrt(log(2)/pi)/alpha_doppler
      g0= 0.83255461 * 0.5641895 / alpha_doppler


c************************** this is the X1= (v-v0)/alpha_doppler*0.8325546
c define arrays for the new Voigt fcn 
      Y1=brd/alpha_doppler*0.83255461
      do 30 ii=1,n_in
c        wr(ii)=0.0
c        wi(ii)=0.0
        X1(ii)= (v(ii)-v0)/alpha_doppler*0.83255461
 30          CONTINUE

      len1=0
      len2=0
      iTruth = +1           !assume all elements of |xr2| < 12 
      do 20 ii=1,n_in
        if ((Y1 .gt. 0.85) .OR. (abs(X1(ii))  .lt. (18.1*Y1+1.65))) then
          len1=len1+1
          region1(len1)=ii
          xr1(len1)=X1(ii)
        else
          len2=len2+1
          region2(len2)=ii
          xr2(len2)=X1(ii)
          if (abs(xr2(len2)) .GT. 12.0) then 
            iTruth = -1 
            end if 
          end if
 20     CONTINUE
      yr1=Y1
      yr2=Y1

c Do region 1 first
      if (len1 .gt. 0 ) then
        do 40 ii=1,len1
          wr1(ii)=0.0
c          wi1(ii)=0.0
 40       CONTINUE
        y1=yr1+1.5
        y2=y1*y1

        do 50 jj=1,6 
          do 60 ii=1,len1
            r=xr1(ii)-T(jj)
            d=1/(r*r+y2)
            d1=y1*d
            d2=r*d
            r=xr1(ii)+T(jj)
            d=1/(r*r+y2)
            d3=y1*d
            d4=r*d
            wr1(ii)=wr1(ii)+C(jj)*(d1+d3)-S(jj)*(d2-d4)
c            wi1(ii)=wi1(ii)+C(jj)*(d2+d4)+S(jj)*(d1-d3)
 60         CONTINUE
 50       CONTINUE
        end if

c Now do region 2

      if (len2 .gt. 0) then
        if (iTruth .GT. 0) then   !all |xr2(ii)| < 12.0 
          do 70 ii=1,len2 
            wr2(ii)=exp(-xr2(ii)*xr2(ii)) 
c            wi2(ii)=0.0 
 70         CONTINUE 
        else                      !some  |xr2(ii)| > 12.0 
          do 75 ii=1,len2 
            wr2(ii)=0.0 
c            wi2(ii)=0.0 
 75         CONTINUE 
          end if 

        y1=yr2+1.5
        y2=y1*y1
        y3=yr2+3

        do 80 jj=1,6 
          do 90 ii=1,len2
            r=xr2(ii)-T(jj)
            d=1/(r*r+y2)
            d1=y1*d
            d2=r*d
            wr2(ii)=wr2(ii)+
     $              yr2*(C(jj)*(r*d2-1.5*d1)+S(jj)*(y3*d2))/(r*r+2.25)
            r=xr2(ii)+T(jj)
            d=1/(r*r+y2)
            d3=y1*d
            d4=r*d
            wr2(ii)=wr2(ii)+
     $              yr2*(C(jj)*(r*d4-1.5*d3)-S(jj)*(y3*d4))/(r*r+2.25)
c            wi2(ii)=wi2(ii)+C(jj)*(d2+d4)+S(jj)*(d1-d3)
 90         CONTINUE
 80       CONTINUE
        end if
 
      if (len1 .gt. 0) then
        do 100 ii=1,len1
c          wr(region1(ii))=wr1(ii)*g0
c          wi(region1(ii))=wi1(ii)*g0
          res1(region1(ii))=wr1(ii)*g0
 100      CONTINUE
        end if

      if (len2 .gt. 0) then
        do 110 ii=1,len2
c          wr(region2(ii))=wr2(ii)*g0
c          wi(region2(ii))=wi2(ii)*g0
          res1(region2(ii))=wr2(ii)*g0
 110      CONTINUE
        end if

c************************** this is the X2= (-v-v0)/alpha_doppler*0.83255461;
c define arrays for the new Voigt fcn 
      Y1=brd/alpha_doppler*0.83255461
      do 130 ii=1,n_in
c        wr(ii)=0.0
c        wi(ii)=0.0
        X1(ii)= (-v(ii)-v0)/alpha_doppler*0.83255461
 130          CONTINUE

      len1=0
      len2=0
      iTruth = +1           !assume all elements of |xr2| < 12 
      do 120 ii=1,n_in
        if ((Y1 .gt. 0.85) .OR. (abs(X1(ii))  .lt. (18.1*Y1+1.65))) then
          len1=len1+1
          region1(len1)=ii
          xr1(len1)=X1(ii)
        else
          len2=len2+1
          region2(len2)=ii
          xr2(len2)=X1(ii)
          if (abs(xr2(len2)) .GT. 12.0) then 
            iTruth = -1 
            end if 
          end if
 120     CONTINUE
       yr1=Y1
       yr2=Y1


c Do region 1 first
      if (len1 .gt. 0 ) then
        do 140 ii=1,len1
          wr1(ii)=0.0
c          wi1(ii)=0.0
 140       CONTINUE
        y1=yr1+1.5
        y2=y1*y1

        do 150 jj=1,6 
          do 160 ii=1,len1
            r=xr1(ii)-T(jj)
            d=1/(r*r+y2)
            d1=y1*d
            d2=r*d
            r=xr1(ii)+T(jj)
            d=1/(r*r+y2)
            d3=y1*d
            d4=r*d
            wr1(ii)=wr1(ii)+C(jj)*(d1+d3)-S(jj)*(d2-d4)
c            wi1(ii)=wi1(ii)+C(jj)*(d2+d4)+S(jj)*(d1-d3)
 160         CONTINUE
 150       CONTINUE
        end if

c Now do region 2

      if (len2 .gt. 0) then
        if (iTruth .GT. 0) then   !all |xr2(ii)| < 12.0 
          do 170 ii=1,len2 
            wr2(ii)=exp(-xr2(ii)*xr2(ii)) 
c            wi2(ii)=0.0 
 170         CONTINUE 
        else                      !some  |xr2(ii)| > 12.0 
          do 175 ii=1,len2 
            wr2(ii)=0.0 
c            wi2(ii)=0.0 
 175        CONTINUE 
          end if 

        y1=yr2+1.5
        y2=y1*y1
        y3=yr2+3

        do 180 jj=1,6 
          do 190 ii=1,len2
            r=xr2(ii)-T(jj)
            d=1/(r*r+y2)
            d1=y1*d
            d2=r*d
            wr2(ii)=wr2(ii)+
     $              yr2*(C(jj)*(r*d2-1.5*d1)+S(jj)*(y3*d2))/(r*r+2.25)
            r=xr2(ii)+T(jj)
            d=1/(r*r+y2)
            d3=y1*d
            d4=r*d
            wr2(ii)=wr2(ii)+
     $              yr2*(C(jj)*(r*d4-1.5*d3)-S(jj)*(y3*d4))/(r*r+2.25)
c            wi2(ii)=wi2(ii)+C(jj)*(d2+d4)+S(jj)*(d1-d3)
 190         CONTINUE
 180       CONTINUE
        end if
 
      if (len1 .gt. 0) then
        do 200 ii=1,len1
c          wr(region1(ii))=wr1(ii)*g0
c          wi(region1(ii))=wi1(ii)*g0
          res2(region1(ii))=wr1(ii)*g0
 200      CONTINUE
        end if

      if (len2 .gt. 0) then
        do 210 ii=1,len2
c          wr(region2(ii))=wr2(ii)*g0
c          wi(region2(ii))=wi2(ii)*g0
          res2(region2(ii))=wr2(ii)*g0
 210      CONTINUE
        end if

c************************** this is the end
      c2=1.4387863        !K/ cm-1  from Genln2 manual

      do 270 ii=1,n_in
        factor=v(ii)*tanh(c2*v(ii)/2/Temp)/(v0*tanh(c2*v0/2/Temp))
        wr(ii)=factor*(res1(ii)+res2(ii))
 270          CONTINUE
        
 
      return
      end

