      subroutine lorentz_vhh(wr,v,v0,Temp,mass,brd,n_in) 
c subroutine voigt1(wr,v,v0,Temp,m,brd,n_in) 
c this is the lorentz lineshape 
c v    = frequency array 
c v0   = center freq 
c T    = temperature 
c m    = molecular mass (amu) 
c brd  = broadening 

      include 'max.inc'
 
      real*8 wr(n_in),v(n_in),v0,Temp,mass,brd 
      integer ii,n_in 
      real*8 brd2 
        
      brd2=brd*brd 
      do ii=1,n_in 
        wr(ii)=0.31830988*brd/(brd2 + (v(ii)-v0)**2) 
        wr(ii)=wr(ii)*v(ii)/v0 
        end do 
 
      return 
      end 
 
 
