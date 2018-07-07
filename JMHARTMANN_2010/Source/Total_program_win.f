c copied from loop_code_umbclnl_fast.f and also from
c /home/sergio/SPECTRA/JMHARTMANN/LM_PQR_CO2_2.0/Source2010/Total_program_win.f

C-------------------------------------------------
      program CalculAbsCO2
C-------------------------------------------------
C-------------------------------------------------
C	Siplest program for computing the Absorption 
C	spectrum of CO2 using our subroutines.
C	The user can easily change the input for the
C	final spectrum.
C
C	INPUT VARIABLES 
C --------------
C	 sgmin [cm-1]
C	 sgmax [cm-1]
C        xCO2   : CO2 volume mixing ratio (Input)  [no unit]
C        xWV    : WV volume mixing  ration (Input) [no unit VMR]
C	 T [K]
C	 p [atm]  total pressure

c these are new from Sergio, for compatibility with run8*.m codes
C        pp   [atm] partial pressure CO2
C        ppw  [atm] partial pressure WV
C        qamt [kilomoles/cm2] gas amount in Genln2 units
C        fout fname for output file
C        lenArray [integer] length of array
C  
C     MixFull = Switch to full diagonalization line-mixing (+1)
c                first order (0) or
c                voigt only (-1)
C     rdmmult = Distance from line center at which you
C               can shift Voigt to lorentz (save CPU time)
C               values lower than 30 should be avoided in
C               order to minimize error
C
C	RESULTS 
C --------------
C     AbsV : Absorption Coefficient neglecting LineMixing
C            (assuming Voigt Line-Shapes) (Cm-1)
C     AbsY : Absorption Coefficient predicted using the First
C            Order Line-Mixing Approximation (Cm-1)
C     AbsW : Absorption Coefficient predicted using Full
C            diagonalization Line-Mixing (Cm-1)
C-------------------------------------------------
C
      implicit none
      include 'parameters.inc'
      integer*4 MixFullX,lenArray,lenArrayX
      integer*4 i,nsig,it
      real*8 Temp,ptot,pCO2,pWV,xCO2,xWV,rdmult,qamt,pathlength
      real*8 sgmin,sgmax,dsg,sig,mgc,stotmax
C Results (Absorption Coefficients)
      real*8 AbsV(nSigmx)
      real*8 AbsY(nSigmx)
      real*8 AbsW(nSigmx)
      character*1 i1,i2
      character*2 cr
      real*8 AbsVtest,AbsYtest,sgTest
      character*80 fout
      real rsgmin,rsgmax,rdsg,rptot,rpCO2,rpWV,rTemp,rqamt
      integer iType,iTest
C----------
C

C Input quantities specification
      sgmin =  600.d0    !!start chunk
      sgmax = 1000.d0    !!stop chunk
      dsg   =  0.01d0    !!output AFTER boxcar

      Temp = 300.d0
      pTot = 1.d0
      xCO2 = 3.5d-4
      xWV  = 1.0d-3

c      print *,'units : cm-1  cm-1  cm-1  atm   atm   atm  K    kmoles/cm2 len iType'
c      print *,'Enter : sgmin sgmax dsg   pTot  pCO2  pWV  Temp qamt len iType:'
      read *,rsgmin,rsgmax,rdsg,rptot,rpCO2,rpWV,rTemp,rqamt,lenArray,iType
c      print *,rsgmin,rsgmax,rdsg,rptot,rpCO2,rpWV,rTemp,rqamt,lenArray,iType

      sgmin = anint(rsgmin)*1.0d0
      sgmax = anint(rsgmax)*1.0d0
      dsg   = rdsg*1.0d0
      if (iType .eq. 20) then
        dsg = 5.00d-4        !!! this is for the 605-2830 cm-1
      else
        print *,'huh dunno which chunk this is!!!!'
        stop
        end if
      ptot  = rptot*1.0d0
      pCO2  = rpCO2*1.0d0
      pWV   = rpWV *1.0d0     !!!! oops had forgotten this, fixed Jan 2011
      Temp  = rTemp*1.0d0
      qamt  = rqamt*1.0d0

      sgmin = dble(anint(rsgmin))
      sgmax = dble(anint(rsgmax))
      dsg   = dble(rdsg)
      if (iType .eq. 20) then
        dsg = 5.00d-4        !!! this is for the 605-2830 cm-1
      else
        print *,'huh dunno which chunk this is!!!!'
        stop
        end if
      ptot  = dble(rptot)
      pCO2  = dble(rpCO2)
      pWV   = dble(rpWV)    !!!! oops had forgotten this, fixed Jan 2011
      Temp  = dble(rTemp)
      qamt  = dble(rqamt)

      sgmin = sgmin - 2*dsg    !!! eg 705.000 starts at 705-2*0.0005
      sgmax = sgmax + (-3)*dsg !!! eg 730.000 ends   at 730+(-5+2)*0.0005
      lenArrayX=dInt( ((SgMax-SgMIn)/DSg)+0.5d0 )+1
      lenArrayX=dInt( ((SgMax-SgMIn)/DSg)+0.0d0 )+1
c      print *,lenArray,lenArrayX !!compare incoming vs computed array length

      xCO2 = pCO2/pTot   !!!vmr
      xWV  = pWV/pTot    !!!vmr

      mgc = 8.314674269981136
      pathlength = qamt*1.0e9*mgc*Temp/(101325*pCO2)

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> three comment lines
c      print *,' single : '
c      write(*,123) rsgmin,rsgmax,rdsg,-1.00D0
c      write(*,123) rptot,rpCO2,rpWV,rTemp
c      print *,' double, after sgmin = rsgmin - 2*dsg etc : '
c      write(*,223) sgmin,sgmax,dsg,-1.00D0
c      write(*,223) ptot,pCO2,pWV,Temp
c      write(*,223) xCO2,xWV,pathlength,-1.0000D0
c      write(*,*) lenarray
c      print *,' '
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> three comment lines

c Losch = (kAtm2mb*100/kBoltzmann/273 * 1e-6) = 2.6895e+19
c      xCO2 = xCO2/2.6895e+19
c      xCO2 = xCO2/1.0e12

c      MixFull=.true.
      MixFullX = 1

      rdmult=30.d0

      iTest = -1
      if (iTest .GT. 0) THEN
        !start testing
        lenArray = 4501
        SgMin = 4750.0D0
        SgMax = 5200.0D0
        Dsg = 0.1D0
        Ptot = 54.0D0
        xCO2 = 0.01D0
        xWV = 0.0D0
        Temp = 218.15D0
        pathlength = 1.0

        lenArray = 89000*5+5
        SgMin = 605.00 - 0.0005*2
        SgMax = 2830 + (-3)*0.0005
        Dsg = 0.0005D0
        Ptot = 54.0D0
        xCO2 = 0.01D0
        xWV = 0.0D0
        Temp = 218.15D0
        pathlength = 1.0
        lenArray=dInt( ((SgMax-SgMIn)/DSg)+0.0d0 )+1
        
        print *,SgMin,SgMax,DSg,Ptot,Ptot*xCO2,Ptot*xWV,Temp,1.000
        end if

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> three comment lines
      write(*,223) sgmin,sgmax,dsg,dble(lenarray*1.0)
      write(*,223) ptot,pCO2,pWV,Temp
      write(*,223) xCO2,xWV,pathlength,dble(lenarray*1.0)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> three comment lines

C --------- 
C Call at the beginning of a program

c Total band intensity cut-off (user supplied)
	sTotMax=0.d0
C
c Call routine to select band according to WN range and cut-off
       call DetBand(sgmin,sgmax,sTotMax)
c read Relaxation matrix files
       call ReadW

C Call for each atmospheric (p,T,VMR) layer
C
	 call CompAbs(Temp,Ptot,xCO2,xWV,SgMin,SgMax,DSg,
     &          AbsV,AbsY,AbsW,MixFullX,lenArray)   !! no rdmult
C     rdmult : Distance (in multirple of Doppler width) at which
C              you can shift Voigt calculation to Lorentz

C ---------

c       open (unit=55,file='./Test.dat',status='unknown')

C Write the results in 'Test.dat'
c recall JM Hartmann's output units are abs coeffs in per cm, so we need to
c multiply this by the path length
      nSig=dInt( ((SgMax-SgMIn)/DSg)+0.5d0 )+1
      nSig=lenArray
      do i=1,Nsig
	sig=sgmin+dfloat(i-1)*Dsg
        if(MixFullX .EQ. 1)then
       	write(*,123) sig,AbsV(i)*pathlength,AbsY(i)*pathlength,
     !                    AbsW(i)*pathlength
        else
       	write(*,124) sig,AbsV(i)*pathlength,AbsY(i)*pathlength
        endif
      enddo

c      close(55)

      stop

c 123  Format(f10.4,3(1pe12.3))
 123  Format(e26.18,3(1pe26.18))
 124  Format(f10.4,2(1pe12.3))
 125  Format(4(1pe12.3))
 223  Format(d26.18,3(1pd26.18))

 1000 Format(1x,'************ PROBLEM !!!! ******************',
     &     /,1x,'Your results differ from the reference',
     &     /,1x,'by a percent value grater than 1% limit',
     &     /,1x,'check you compiler and other error source',//)

      end
C-------------------------------------------------

C
C******************************************************************************
