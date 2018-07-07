      subroutine dofull(z,v,v0,w_tot,temperature,duration,pressure_for,
     $    pressure_self,arg1,arg2,diagL,ksm,bsm,birn,lenf,no_lines)

c doe thast is called is dofullNEWBIRN.f
c this code is NOT called. also see that below we do 
c k(i)=sum(full mix lines) * birn * sumlorbirn/sumlor

      include 'max.inc'

      integer lenf,no_lines

      real*8 z(lenf),ksm,bsm,v(lenf),birn
      real*8 w_tot(no_lines),v0(no_lines)
      complex*16 arg1(no_lines),arg2(no_lines),diagL(no_lines)

      real*8 chi(MaxLen),sumlor(MaxLen),sumlorbirn(MaxLen),lor(MaxLen)
      complex*16 temp(MaxLen)

      real*8 temperature,duration,pressure_for,pressure_self,ksmbsm

      integer i,j

      ksmbsm=ksm*bsm

      if (lenf .gt. MaxLen) then
        print *,'in dofull2.f need MaxLen > lenf'
        stop
        end if

      do j=1,lenf
        z(j)=0.0
        end do

      if (abs(birn-1.0) .le. 0.01) then !have to do birnbaum
        do j=1,lenf
          chi(j)=1.0
          sumlor(j)=0.0
          sumlorbirn(j)=0.0
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
  
        if (abs(birn-1.0) .le. 0.01) then        !birnbaum 
          call birnbaum(chi,v,v0(i),w_tot(i),temperature,duration,
     $                   lenf)
          call lorentz(lor,v,v0,temperature,1.0,w_tot,lenf)
          do j=1,lenf
            sumlor(j)=sumlor(j)+lor(j)
            sumlorbirn(j)=sumlorbirn(j)+lor(j)*chi(j)
            end do
        elseif (abs(birn+1.0) .le. 0.01) then   
          !cousin  .... not possible in FULL mix
          call cousin(chi,v,v0(i),w_tot(i),temperature,
     $          pressure_for,pressure_self,lenf) 
          endif
 
        if (abs(birn+1.0) .le. 0.01) then   !cousin 
          do j=1,lenf
            z(j)=z(j)+imag(temp(j))*chi(j)
            end do
        else                                     !no chi factor, or birnbaum
          do j=1,lenf
            z(j)=z(j)+imag(temp(j))
            end do
          endif

        enddo 

      if (abs(birn-1.0) .le. 0.01) then        !birnbaum 
        do j=1,lenf 
          z(j)=z(j)*v(j)*ksmbsm*sumlorbirn(j)/sumlor(j) 
          end do 
      else           !no chi function
        do j=1,lenf 
          z(j)=z(j)*v(j)*ksmbsm 
          end do
        endif

      return
      end

c ***********************************************************************
      include 'lorentz.f'
      include 'cousin1.f'
      include 'birnbaum.f'
