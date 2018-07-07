      subroutine doVmixRatio(z,v,v0,w_tot,temperature,duration,
     $      pressure_for,pressure_self,mass,strenqt,ymix,kvm,klm,
     $      bsm,birn,NIFNIF,ratio,lenf,no_lines)

      include 'max.inc'

      integer lenf,no_lines

      real*8 z(MaxLen),kvm,klm,bsm,v(MaxLen),birn,NIFNIF,mass
      real*8 w_tot(MaxPQR),v0(MaxPQR),one
      real*8 strenqt(MaxPQR),ymix(MaxPQR),ratio

      real*8 chi(MaxLen),c2,c2_0,fact,factor
      real*8 wr(MaxLen),wi(MaxLen),temp(MaxLen)

      real*8 lor(MaxLen),sumlor(MaxLen),sumlorbirn(MaxLen)

      real*8 temperature,duration,pressure_for,pressure_self,ksmbsm

      integer i,j,nif

      one=1.0d0
      c2=1.4387863        !K/ cm-1  from Genln2 manual  

      nif=int(NIFNIF)
      if (nif .eq. 1) then
        ksmbsm=kvm*bsm
      else
        ksmbsm=klm*bsm
        endif

      if (lenf .gt. MaxLen) then
        print *,'in doVmixRatio.f need MaxLen > lenf'
        stop
        endif

      if (no_lines .gt. MaxPQR) then
        print *,'in doVmixRatio.f need MaxPQR > no_lines'
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

c      print *,'k=klor*ratio : ratio,birn = ',ratio,birn
      do i=1,no_lines
        c2_0=0.5*c2*v0(i)/temperature
        c2_0=v0(i)*tanh(c2_0)  
        fact=1/c2_0  
        c2=0.5*c2/temperature

        if (nif .eq. 1) then 
c this is the fast complex voigter : use voigt fcn for mixing 
c note this ALREADY does the van huber shape
          call vhh2RI(wr,wi,v,v0(i),temperature,mass,w_tot(i),lenf) 
          do j=1,lenf 
            temp(j)=ratio*strenqt(i)*wr(j)
            end do 
        elseif (nif .eq. -1) then 
c do lorentz : use lorentz for mixing 
          do j=1,lenf 
            temp(j)=ratio*strenqt(i)/3.141592653589793 * 
     $               w_tot(i)/((v(j)-v0(i))**2+(w_tot(i))**2)
ccc            factor=v(j)/v0(i)               !only have the v/vo correction 
               factor=1.0       !lorentz should only have this
ccc   this was always commented out
ccc   this is the complete VHH correction 
ccc         factor=fact*v(j)*tanh(0.5*c2*v(j)/temperature) 
            temp(j)=factor*temp(j)
            end do
          end if
 
       if (abs(birn-1.0) .le. 0.01) then        !birnbaum 
          call birnbaum(chi,v,v0(i),w_tot(i),temperature,duration,
     $                   lenf)
          call lorentz(lor,v,v0(i),temperature,one,w_tot(i),lenf)
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
      if (abs(birn-1.0) .ge. 0.01) then   !birn=-1 or 0 ==> cousin or no chi

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
      include 'vhh2RI.f'
      include 'lorentz.f'
      include 'cousin1.f'
      include 'birnbaum.f'
