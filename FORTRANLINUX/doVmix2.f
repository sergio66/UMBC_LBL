c      subroutine doVmix2(z,v,v0,w_tot,temperature,duration,pressure_for, 
c     $      pressure_self,mass,strenqt,ymix,kvm,klm,bsm,birnbirn,NIFNIF, 
c     $      c0,cspan,lenf,no_lines) 

      subroutine doVmix2(z,v,v0,w_tot,temperature,duration,pressure_for,
     $      pressure_self,mass,strenqt,ymix,kvm,klm,bsm,birn,NIFNIF,
     $      c0,cspan,lenf,no_lines)

c same as doVmix.f except that it has additional 2 parameters for
c the gaussian chi function ychi = c0 * exp(-((vo(i)-v)/cspan)**2)
c this chi function has c0 > 0, and is between c0 and 1 (ie the chifunction
c INTENSIFIES the lineshape. this is so that our mix code gives results that
c are compatible with AIRS and JOHN JOHN data)

c       !no need to do the v/v0(i)*tanh(v/T)/tanh(v0(i)/T) because this
c       !is already done by the vhh mex file

      include 'max.inc'

      integer lenf,no_lines

      real*8 z(MaxLen),kvm,klm,bsm,v(MaxLen),birn,NIFNIF,mass,c0,cspan
      real*8 w_tot(MaxPQR),v0(MaxPQR)
      real*8 strenqt(MaxPQR),ymix(MaxPQR)

      real*8 chi(MaxLen)
      real*8 wr(MaxLen),wi(MaxLen),temp(MaxLen)

      real*8 temperature,duration,pressure_for,pressure_self,ksmbsm,newchi

      integer i,j,V_or_L,IO

      V_or_L=int(NIFNIF)
      if (V_or_L .eq. 1) then
        ksmbsm=kvm*bsm
      else
        ksmbsm=klm*bsm
        endif

      print *,'doc,c0,cspan = ',duration,c0,cspan
      print *,'T,pf,ps = ',temperature,pressure_for,pressure_self
      print *,'lenf,no_lines = ',lenf,no_lines

      if (lenf .gt. MaxLen) then
        print *,'in doVmix2.f need MaxLen > lenf'
        stop
        endif

      if (no_lines .gt. MaxPQR) then
        print *,'in doVmix2.f need MaxPQR > no_lines'
        stop
        endif

      IO=+1                    !assume we are doing mixing               
      chi(1)=0.0
      do i=1,no_lines
        chi(1)=chi(1)+abs(ymix(i))
        end do
      if (chi(1) .le. 1.0e-6) then     !all mix coeffs were zero ==> no mixing
        IO = -1
        endif
      if ((IO .gt. 0) .and. (abs(birn+1.0) .le. 0.01)) then    !birn=-1
        !cousin and first order=WRONG
        print *,'cousin and first order incompatible!!'
        stop
        endif
      chi(1)=1.0

      do j=1,lenf
        z(j)=0.0
        end do

      if (abs(int(birn)) .le. 0.01) then
        do j=1,lenf
          chi(j)=1.0
          end do
        endif

      print *,'in doVmix2.f, IO,V_or_L, = ',IO,V_or_L

      do i=1,no_lines
        if (IO .gt. 0) then !have to do first order mixing
          if (V_or_L .eq. 1) then
c this is the fast complex voigter : use voigt fcn for mixing
c make sure the chifunction is >= 1
            call vhh2RI(wr,wi,v,v0(i),temperature,mass,w_tot(i),lenf)
            do j=1,lenf
              !ymix=zeros if IO set to -1
              temp(j)=strenqt(i)*(wr(j)+ymix(i)*wi(j))
              end do
          elseif (V_or_L .eq. -1) then
c do lorentz : use lorentz for mixing
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
              !ymix=zeros if IO set to -1
              temp(j)=strenqt(i)*(wr(j))
              end do
          elseif (V_or_L .eq. -1) then
            print *,'start line number ',i, ' of total ',no_lines
            do j=1,lenf
              temp(j) = 
     $            strenqt(i)/3.141592653589793 *(v(j)/v0(i))*w_tot(i)
     $                 /((v(j)-v0(i))**2+(w_tot(i))**2)
              end do
            print *,'end line number ',i, ' of total ',no_lines
          else
            print *,'need to do either voigt or lorentz'
            stop
            endif
          endif

       if (abs(birn-1.0) .le. 0.01) then        !birn=+1 ==> birnbaum
          call birnbaum(chi,v,v0(i),w_tot(i),temperature,duration,
     $                   lenf)
        elseif (abs(birn+1.0) .le. 0.01) then   !birn=-1 ==> cousin 
          call cousin(chi,v,v0(i),w_tot(i),temperature,
     $          pressure_for,pressure_self,lenf) 
          endif

        if (abs(birn) .le. 1e-6) then           !no chi factor 
          do j=1,lenf
            z(j)=z(j)+temp(j)
            end do
        else                                    !include chi factor
          !also put in the c0 * exp(-((v-v0)/sigma ** 2))
          do j=1,lenf
            newchi = (v(j)-v0(i))/cspan
            newchi = exp(-(newchi * newchi) * 1.0d0)
            newchi = max(1.0d0, c0*1.0d0 * newchi)
            !z(j)=z(j)+temp(j)*chi(j)            !only birnbaum
            !z(j)=z(j)+temp(j)*chi(j)*newchi     !birnbaum times our modify
            z(j)=z(j)+temp(j)*newchi            !only our modify
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
