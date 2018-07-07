      subroutine doVmixSimple(z,v,v0,w_tot,temperature,duration,
     $      pressure_for,pressure_self,mass,strenqt,ymix,kvm,klm,bsm,
     $      birn,NIFNIF,IOIO,theRATIO,lenf,no_lines)

c       !no need to do the v/v0(i)*tanh(v/T)/tanh(v0(i)/T) because this
c       !is already done by the vhh mex file

      include 'max.inc'

      integer lenf,no_lines

      real*8 z(MaxLen),kvm,klm,bsm,v(MaxLen),birn,NIFNIF,IOIO,mass
      real*8 w_tot(MaxPQR),v0(MaxPQR),theRATIO
      real*8 strenqt(MaxPQR),ymix(MaxPQR)

      real*8 chi(MaxLen)
      real*8 wr(MaxLen),wi(MaxLen),temp(MaxLen)

      real*8 temperature,duration,pressure_for,pressure_self,ksmbsm

      real*8 v0min,v0max

      integer i,j,V_or_L,nio,IO

      nio=INT(ioio)

c      print *,'in subroutine MEX doVmixSimple : theRATIO = ',theRATIO
c      print *,'in subroutine MEX doVmixSimple : nio = ',nio

      V_or_L=int(NIFNIF)
      if (V_or_L .eq. 1) then
        ksmbsm=kvm*bsm
      else
        ksmbsm=klm*bsm
        endif

      if (lenf .gt. MaxLen) then
        print *,'in doVmix.f need MaxLen > lenf',MaxLen,lenf
        stop
        endif

      if (no_lines .gt. MaxPQR) then
        print *,'in doVmixSimple.f need MaxPQR > no_lines'
        stop
        endif

      do j=1,lenf
        z(j)=0.0
        end do

      if (abs(int(birn)) .le. 0.01) then
        do j=1,lenf
          chi(j)=1.0
          end do
        endif

      v0min = +1.0e6
      v0max = -1.0e6
      do i=1,no_lines
        if (v0(i) .lt. v0min) then
          v0min = v0(i)
          endif
        if (v0(i) .gt. v0max) then
          v0max = v0(i)
          endif
        end do
c      print *,'v0min,v0max = ',v0min,v0max

      IO=+1                    !assume we are doing mixing               
      chi(1)=0.0
      do i=1,no_lines
        chi(1)=chi(1)+abs(ymix(i))
        end do
      if (chi(1) .le. 1e-6) then     !all mix coeffs were zero ==> no mixing
        IO = -1
        end if
      if ((IO .gt. 0) .and. (abs(int(birn)+1.0) .le. 0.01)) then   
        !cousin and first order=WRONG
        print *,'cousin and first order incompatible!!'
        stop
        endif
      chi(1)=0.0

      do i=1,no_lines
        if (IO .gt. 0) then !have to do first order mixing
          if (V_or_L .eq. 1) then
c this is the fast complex voigter
            call vhh2RI(wr,wi,v,v0(i),temperature,mass,w_tot(i),lenf)
            do j=1,lenf
              !ymix=zeros if IO set to 0
              temp(j)=strenqt(i)*(wr(j)+ymix(i)*wi(j))
              end do
          elseif (V_or_L .eq. -1) then
            do j=1,lenf
              temp(j) = 
     $          strenqt(i)/3.141592653589793 * (v(j)/v0(i))*(w_tot(i)+
     $          (v(j)-v0(i))*ymix(i))/((v(j)-v0(i))**2+(w_tot(i))**2)
              end do
          else
            print *,'need to do either voigt or lorentz'
            stop
            end if
        else           !no mixing
          if (V_or_L .eq. 1) then
c this is the fast complex voigter
            call vhh2RI(wr,wi,v,v0(i),temperature,mass,w_tot(i),lenf)
            do j=1,lenf
              !ymix=zeros if IO set to 0
              temp(j)=strenqt(i)*(wr(j))
              end do
          elseif (V_or_L .eq. -1) then
            do j=1,lenf
              temp(j) = 
     $           strenqt(i)/3.141592653589793 *(v(j)/v0(i))*w_tot(i)
     $                 /((v(j)-v0(i))**2+(w_tot(i))**2)
              end do
          else
            print *,'need to do either voigt or lorentz'
            stop
            end if
          end if

        if (nio .gt. 0) then !simple linemixing
          do j=1,lenf
            if (abs(v(j)-v0(i)) .gt. 15.0) then   !orig code
              temp(j)=theRATIO*temp(j)
              endif
            end do
          endif

c this is new : to account for weak PR linemixing in between lines!
c        if (nio .gt. 0) then !simple linemixing in between lines          
c          do j=1,lenf
c            if ((v(j) .ge. v0min) .and. (v(j) .le. v0max) .and. 
c     $           (abs(v(j)-v0(i)) .gt. w_tot(i)*5)) then
c              temp(j)=temp(j)*0.5
c              endif
c            end do
c          endif

c this is computing either birnbaum or cousin chi fcn
       if (abs(int(birn)-1.0) .le. 0.01) then        !birnbaum 
          call birnbaum(chi,v,v0(i),w_tot(i),temperature,duration,
     $                   lenf)
        elseif (abs(int(birn)+1.0) .le. 0.01) then   !cousin 
          call cousin(chi,v,v0(i),w_tot(i),temperature,
     $          pressure_for,pressure_self,lenf) 
          endif

        if (abs(birn) .le. 1e-6) then           !no chi factor 
          do j=1,lenf
            z(j)=z(j)+temp(j)
            end do
        else                                    !include chi factor
          do j=1,lenf
            z(j)=z(j)+temp(j)*chi(j)
            end do
          end if

        enddo

c both voigt and lorentz already have the v/vo(j) factor 
      do j=1,lenf 
        z(j)=z(j)*ksmbsm 
        end do

      return
      end

c************************************************************************
      include 'vhh2RI.f' 
      include 'cousin1.f' 
      include 'birnbaum.f' 
