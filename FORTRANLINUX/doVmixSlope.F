      subroutine doVmixSlope(z,v,v0,w_tot,temperature,duration,
     $      pressure_for,pressure_self,mass,strenqt,ymix,kvm,klm,
     $      bsm,birn,NIFNIF,ratio,slope,endpt,lenf,no_lines)

      include 'max.inc'

      integer lenf,no_lines

      real*8 z(lenf),kvm,klm,bsm,v(lenf),birn,NIFNIF,mass
      real*8 w_tot(no_lines),v0(no_lines)
      real*8 strenqt(no_lines),ymix(no_lines)
      real*8 rr,ratio,slope,endpt

      real*8 chi(MaxLen),c2,c2_0,fact,factor,repwid
      real*8 wr(MaxLen),wi(MaxLen),temp(MaxLen)

      real*8 lor(MaxLen),sumlor(MaxLen),sumlorbirn(MaxLen)

      real*8 temperature,duration,pressure_for,pressure_self,ksmbsm

      integer i,j,nif

      c2=1.4387863        !K/ cm-1  from Genln2 manual  

      nif=int(NIFNIF)
      if (nif .eq. 1) then
        ksmbsm=kvm*bsm
      else
        ksmbsm=klm*bsm
        endif

      if (lenf .gt. MaxLen) then
        print *,'in doVmix2.f need MaxLen > lenf'
        stop
        endif

      do j=1,lenf
        z(j)=0.0
        end do

      if (abs(int(birn)) .le. 0.01) then         !birn= 0 ==> no chi fcn
        do j=1,lenf
          chi(j)=1.0
          end do
      elseif (abs(birn-1.0) .le. 0.0) then  !birn=+1  ==> use birnbaum
        do j=1,lenf
          sumlor(j)=0.0
          sumlorbirn(j)=0.0
          end do
c      elseif (abs(birn+1.0) .le. 0.0) then  !birn=-1  ==> use cousin
c        do j=1,lenf
c          chi(j)=1.0
c          end do
        endif


      do i=1,no_lines
        !nif ===== 0 always here ==> do lorentz*ratio

        c2_0=0.5*c2*v0(i)/temperature
        c2_0=v0(i)*tanh(c2_0)  
        fact=1/c2_0  
        c2=0.5*c2/temperature

        do j=1,lenf
          rr=ratio-(endpt-v(j))*slope
          temp(j)=rr*strenqt(i)/3.141592653589793 *w_tot(i)/
     $               ((v(j)-v0(i))**2+(w_tot(i))**2)
          factor=v(j)/v0(i)               !only have the v/vo correction 
c          factor=fact*v(j)*tanh(c2*v(j)) !have complete VHH correction 
          temp(j)=factor*temp(j)  
          end do
 
       if (abs(birn-1.0) .le. 0.01) then        !birnbaum 
          call birnbaum(chi,v,v0(i),w_tot(i),temperature,duration,
     $                   lenf)
          call lorentz(lor,v,v0(i),temperature,1.0,w_tot(i),lenf)
          do j=1,lenf
            sumlor(j)=sumlor(j)+lor(j)
            sumlorbirn(j)=sumlorbirn(j)+lor(j)*chi(j)
            enddo
        elseif (abs(birn+1.0) .le. 0.01) then   !cousin 
          call cousin(chi,v,v0(i),w_tot(i),temperature,
     $          pressure_for,pressure_self,lenf) 
          endif
 
       if (abs(birn) .le. 0.01) then        !no chi         
         do j=1,lenf
           z(j)=z(j)+temp(j)
           end do
       elseif (abs(birn-1.0) .le. 0.01) then        !birnbaum
         do j=1,lenf
           z(j)=z(j)+temp(j)
           end do
       elseif (abs(birn+1.0) .le. 0.01) then        !cousin
         do j=1,lenf
           z(j)=z(j)+temp(j)*chi(j)
           end do
         endif

        enddo

c both voigt and lorentz already have the v/vo(j) factor 
      if (abs(birn-1.0) .ge. 0.01) then !birn=-1 or 0 ==> cousin or no chi
        do j=1,lenf 
          z(j)=z(j)*ksmbsm 
          end do
      else
        do j=1,lenf 
          z(j)=z(j)*ksmbsm*sumlorbirn(j)/sumlor(j)
          end do
        endif

      return
      end

c ***********************************************************************
      include 'lorentz.f'
      include 'cousin1.f'
      include 'birnbaum.f'
