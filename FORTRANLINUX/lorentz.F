      subroutine lorentz(wr,v,v0,Temp,m,brd,n_in) 
c subroutine lorentz(wr,v,v0,Temp,m,brd,n_in) 
c this is the lorentz lineshape
c v    = frequency array 
c v0   = center freq 
c T    = temperature 
c m    = molecular mass (amu) 
c brd  = broadening 
 
      real*8 wr(n_in),v(n_in),v0,Temp,m,brd 
      integer ii,n_in 
      real*8 brd2,pi_brd 
 
      brd2=brd*brd 
      pi_brd=0.31830988*brd 
 
c  1/pi = 0.318... 
      do ii=1,n_in 
        wr(ii)=pi_brd/(brd2 + (v(ii)-v0)*(v(ii)-v0)) 
        enddo 
      return 
      end 
 

c************************************************************************
      subroutine lorentz_basement(wr,v,v0,Temp,m,brd,n_in) 
c subroutine lorentz(wr,v,v0,Temp,m,brd,n_in) 
c this is the lorentz lineshape, except that to use CKD defn of water 
c continuum, subtract off the basement term

c v    = frequency array 
c v0   = center freq 
c T    = temperature 
c m    = molecular mass (amu) 
c brd  = broadening 
 
      real*8 wr(n_in),v(n_in),v0,Temp,m,brd 
      integer ii,n_in 
      real*8 brd2,pi_brd,basmnt
 
      brd2=brd*brd 
      pi_brd=0.31830988*brd 

c      1/PI = 0.31830988
c      625.0 = 25.0**2
       basmnt = pi_brd/(625.0 + brd2)
 
c  1/pi = 0.318... 
      do ii=1,n_in 
        wr(ii)=pi_brd/(brd2 + (v(ii)-v0)*(v(ii)-v0)) - basmnt
        enddo 

      return 
      end 
 
c************************************************************************
      subroutine lorentz_base2(wr,v,v0,Temp,m,brd,n_in,basmnt) 
c subroutine lorentz(wr,v,v0,Temp,m,brd,n_in) 
c this is the lorentz lineshape, except that to use CKD defn of water 
c continuum, subtract off the basement term, which has been precomputed using
c the VHH lineshape. also note, outside of 25 cm-1, wr is explicitly set to 0

c v    = frequency array 
c v0   = center freq 
c T    = temperature 
c m    = molecular mass (amu) 
c brd  = broadening
c basement = precomputed basement term 
 
      real*8 wr(n_in),v(n_in),v0,Temp,m,brd
      integer ii,n_in 
      real*8 brd2,pi_brd,basmnt
 
      do ii=1,n_in 
        if (abs(v(ii)-v0) .le. 25.0) then
          wr(ii)=pi_brd/(brd2 + (v(ii)-v0)*(v(ii)-v0)) - basmnt
        else
          wr(ii)=0.0
          end if
        enddo 

      return 
      end 
 
c************************************************************************

