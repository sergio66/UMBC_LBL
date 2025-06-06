      subroutine vhh(z,v,v0,T,m,brd,n_in)
c subroutine vhh(z,v,v0,T,m,brd,n_in)
c this is the van huber lineshape : using the polynomial approx
c v    = frequency array
c v0   = center freq
c T    = temperature
c m    = molecular mass (amu)
c brd  = broadening
c z    = results array

      real*8 z(n_in),v(n_in),v0,T,m,brd

      real*8 k,c_light,amu,mass,r2,alpha_doppler,g0
      real*8 X1(n_in),X2(n_in),Z1(n_in),Z2(n_in)
      real*8 Y1,a(4),b(4),c(4),d(4),temp,c2,factor
      real*8 c2_0,fact

      integer ii,jj,n_in

c do the doppler widths first 
      k=1.380658e-23
      c_light=2.99792458e8         !ms-1
      amu=1.6605402e-27            !nucleon mass/kg
      mass=m                       !change to kg

c alpha_doppler=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light); 
c r2=2*log(2)*k/amu; 
      r2=11526.218
      alpha_doppler=v0/c_light*sqrt(r2*T/mass)

c do the g0 factor 
c  SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.5641895 
c g0=sqrt(log(2)/pi)/alpha_doppler;
      g0= 0.83255461 * 0.5641895 / alpha_doppler

c define arrays for the new Voigt fcn 
      Y1=brd/alpha_doppler*0.83255461
      do 30 ii=1,n_in
        Z1(ii)=0.0
        Z2(ii)=0.0
        X1(ii)= (v(ii)-v0)/alpha_doppler*0.83255461
        X2(ii)= (-v(ii)-v0)/alpha_doppler*0.83255461
 30     CONTINUE


      a(1)= -1.21498213255730
      a(2)= -1.35094358543273
      a(3)= -1.21498213255730 
      a(4)= -1.35094358543273


      b(1)= 1.23588765343593
      b(2)=  0.37861161238627
      b(3)= -1.23588765343593
      b(4)= -0.37861161238627

      c(1)= -0.30849709263498
      c(2)=  0.59059188440886
      c(3)= -0.30849709263498
      c(4)=  0.59059188440886

      d(1)= 0.02098588080036
      d(2)= -1.18584332504043
      d(3)= -0.02098588080036
      d(4)=  1.18584332504043


      do 40 ii=1,n_in
        do 50 jj=1,4
          temp=(c(jj)*(y1-a(jj))+d(jj)*(X1(ii)-b(jj)))
          temp=temp/((y1-a(jj))**2 + (X1(ii)-b(jj))**2)
          Z1(ii)=Z1(ii)+temp

          temp=(c(jj)*(y1-a(jj))+d(jj)*(X2(ii)-b(jj)))
          temp=temp/((y1-a(jj))**2 + (X2(ii)-b(jj))**2)
          Z2(ii)=Z2(ii)+temp
 50       CONTINUE
 40     CONTINUE

      c2=1.4387863        !K/ cm-1  from Genln2 manual

      c2_0=0.5*c2*v0/Temp 
      c2_0=v0*tanh(c2_0) 
      fact=g0/c2_0 
 
      c2=0.5*c2/Temp 
 
      do 270 ii=1,n_in 
        factor=fact*v(ii)*tanh(c2*v(ii)) 
        z(ii)=factor*(Z1(ii)+z2(ii)) 
 270    continue 

      return
      end


