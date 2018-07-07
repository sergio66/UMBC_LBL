      subroutine dofullNEWBIRN(z,v,v0,w_tot,temperature,duration,
     $             pressure_for,pressure_self,arg1,arg2,diagL,ksm,bsm,
     $             birn,lenf,no_lines)
 
c this is the code that is called
c same as dofull.f, except that birn is [ind "line" *  birn] instead of
c                                       sum(lines)*(birn*lor)/lor
      include 'max.inc'

      integer lenf,no_lines

      real*8 z(MaxLen),ksm,bsm,v(MaxLen),birn
      real*8 w_tot(MaxPQR),v0(MaxPQR)
      complex*16 arg1(MaxPQR),arg2(MaxPQR),diagL(MaxPQR)

      real*8 chi(MaxLen)
      complex*16 temp(MaxLen)

      real*8 temperature,duration,pressure_for,pressure_self,ksmbsm

      integer i,j

      ksmbsm=ksm*bsm

c      if (lenf .gt. MaxLen) then
c        print *,'in dofullnewbirn.f need MaxLen > lenf'
c        stop
c        end if

c      if (no_lines .gt. MaxPQR) then
c        print *,'in dofullfullnewbirn.f need MaxPQR > no_lines'
c        stop
c        end if

      do j=1,lenf
        z(j)=0.0
        end do

      if (abs(birn-1.0) .le. 0.01) then !have to do birnbaum
        do j=1,lenf
          chi(j)=1.0
          end do
      else
        do j=1,lenf
          chi(j)=1.0
          end do
        endif

      do i=1,no_lines
        do j=1,lenf
          temp(j)=arg1(i)*arg2(i)/(v(j)-diagL(i))
          end do
  
        if (abs(birn-1.0) .le. 0.01) then   !birnbaum 
          call birnbaum(chi,v,v0(i),w_tot(i),temperature,duration,
     $                   lenf)
          do j=1,lenf
            z(j)=z(j)+dimag(temp(j))*chi(j)
            end do
        elseif (abs(birn+1.0) .le. 0.01) then   !cousin 
          !cousin  .... not possible in FULL mix
          call cousin(chi,v,v0(i),w_tot(i),temperature,
     $          pressure_for,pressure_self,lenf) 
          do j=1,lenf
            z(j)=z(j)+dimag(temp(j))*chi(j)
            end do
        else                                     !no chi factor, or birnbaum
          do j=1,lenf
            z(j)=z(j)+dimag(temp(j))
            end do
          endif

        enddo 

      do j=1,lenf 
        z(j)=z(j)*v(j)*ksmbsm 
        end do

      return
      end

c ***********************************************************************
      include 'lorentz.f'
      include 'cousin1.f'
      include 'birnbaum.f'
