C-------------------------------------------------
      program CalculAbsCO2
C-------------------------------------------------
C-------------------------------------------------
C same as loop_code,f except it only computes fullinemix or voigt or first
c order, as necessary. Use "integer MixFullX" instead of "logical MixFull"

C	Siplest program for computing the Absorption 
C	spectrum of CO2 using our subroutines.
C	The user can easily change the input for the
C	final spectrum.
C
C	INPUT VARIABLES 
C --------------
C	 sgmin [cm-1]  start wavenumber
C	 sgmax [cm-1]  stop  wavenumber
C	 dsg   [cm-1]  delta wavenumber AFTER boxcar integration
C	 xCO2 [no unit VMR]
C        xWV  [no unit VMR]
C	 T    [K]   temperature
C	 p    [atm] total pressure

c these are new from Sergio, for compatibility with run8*.m codes
C        pp     [atm] CO2 partial pressure
C        ppWV   [atm] WV partial pressure
C        qamt [kilomoles/cm2] gas amount in Genln2 units
C        fout fname for output file
C        lenArray [integer] length of array
C  
C     MixFullX = Switch to full diagonalization line-mixing (+1)
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
C            (assuming Voigt Line-Shapes) (Cm-1) * pathlength(cm) = no units
C     AbsY : Absorption Coefficient predicted using the First
C            Order Line-Mixing Approximation (Cm-1) * pathlength(cm) = no units
C     AbsW : Absorption Coefficient predicted using Full
C            diagonalization Line-Mixing (Cm-1) * pathlength(cm) = no units
C-------------------------------------------------
C
      implicit none
      include 'parameters.inc'
      integer*4 MixFullX,lenArray
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

cc      print *,'units : cm-1  cm-1  cm-1  atm   atm   atm K    kilomoles/cm2 '
cc      print *,'Enter : sgmin sgmax dsg   pTot  pCO2  pWV Temp qamt    : '
      read *,sgmin,sgmax,dsg,pTot,pCO2,pWV,Temp,qamt
cc      print *,'Enter name of output file : '
      read *,fout
cc      print *,'Enter   MixFullX     LenArray : '
      read *,MixFullX, lenArray
c      MixFullX = 2
c      lenArray = 50000+2+3
      print  *,sgmin,sgmax,dsg,pTot,pCO2,pWV,pWV,Temp,qamt
      print *,fout
      print *,MixFullX, lenArray

C      dsg  = dsg/nbox       !!!change to fine resolution
C      sgmin = sgmin - 2*dsg !!!so that eg 705.000 starts at 705-2*0.0005 
C      sgmax = sgmax - 3*dsg !!!so that eg 730.000 ends   at 730+(-5+2)*0.0005 
      xCO2 = pCO2/pTot   !!!vmr
      xWV  = pWV/pTot   !!!vmr

      mgc = 8.314674269981136
      pathlength = qamt*1.0e9*mgc*Temp/(101325*pCO2)

c Losch = (kAtm2mb*100/kBoltzmann/273 * 1e-6) = 2.6895e+19
c      xCO2 = xCO2/2.6895e+19
c      xCO2 = xCO2/1.0e12

c      MixFull=.true.
      rdmult=30.d0

C --------- 
C Call at the beginning of a program
C
c Total band intensity cut-off (user supplied)
	sTotMax=0.d0

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

C Write the results in 'Test.dat' .. changed to fout
       print *,fout
       open (unit=55,file=fout,status='unknown')

c recall JM Hartmann's output units are abs coeffs in per cm, so we need to 
c multiply this by the path length

      nSig=dInt( ((SgMax-SgMIn)/DSg)+0.5d0 )+1
      nSig=dInt( ((SgMax-SgMIn)/DSg))+1
      nSig=lenArray
      do i=1,Nsig
	sig=sgmin+(i-1)*Dsg
        if(MixFullX .EQ. -1) then
          write(55,125) sig,AbsV(i)*pathlength
        elseif(MixFullX .EQ. 0) then
          write(55,125) sig,AbsY(i)*pathlength
        elseif(MixFullX .EQ. +1) then
          write(55,125) sig,AbsW(i)*pathlength
        elseif(MixFullX .EQ. +2) then
          write(55,125) sig,AbsV(i)*pathlength,AbsY(i)*pathlength,AbsW(i)*pathlength
        endif
      enddo

      close(55)

      stop

 123  Format(f10.4,3(1pe12.3))
 124  Format(f10.4,2(1pe12.3))
 125  Format(f10.4,1(1pe12.3))

 1000 Format(1x,'************ PROBLEM !!!! ******************',
     &     /,1x,'Your results differ from the reference',
     &     /,1x,'by a percent value grater than 1% limit',
     &     /,1x,'check you compiler and other error source',//)

      end
C-------------------------------------------------

