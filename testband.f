! g77 -o testband.x testband.f
! see http://home.pcisys.net/~bestwork.1/CalcAbs/CalcAbsHitran.html

!     Program Reads HiTran data; calculates a table of absorption cross sections
!     HITRAN unit of line strength S of cm-1/(molecules cm-2)
!     Lorentz lineshape L has units of 1/gamma = 1/cm-1
!     So S L has units 1/(molecules cm-2) = cm2/molecule, and then we sum over all lines
!     ** Output units : cross section / molecule     [cm2/molecule] ** 
!
!      so to get eg [m2/mole] multiply by (100*100)^(-1) *Na  Na=avogadro
!      ie multiply by 6.023e23/(10000) = 6.023e19
!      ha = load('CO2abs.csv'); 
!      plot(ha(:,1),ha(:,2)*1e-4*6e23,'r');
!        xlabel('Wavenumber cm-1'); ylabel('cross section m2/mole')

c************************************************************************
       PROGRAM CalcXn 
       PARAMETER (MAXLINES=100000,NVALUES=5) 
       REAL*8 VALX(NVALUES),VAL(MAXLINES,NVALUES),DELTAKAPPA
       REAL*8 KAPPA,KAPPAHI,KAPPALO
       REAL*8 S,V,ILNZ,PI,et 
       PARAMETER  (PI=3.1415926535897932385D0) 

!      Assign starting and stopping points in the data range. 
!      Absorbance values are calculated at interval 2*DELTAKAPPA and 
!      start with KAPPA+DELTAKAPPA 
!      The output will be in units of cm^2/molecule. Essentially 
!      the cross section for absorption over the spectral bandwidth 2*DELTAKAPPA

       DELTAKAPPA=.0603/2        
       GID = 1; 
       GID = 2;
       print *,'ENTER GID (1,2,3,4,5,6,9,12) : '
       read(5,*) GID
       
!      #1: Input file
       IF (GID .EQ. 1) THEN
         OPEN (1,FILE='/asl/s1/motteler/read_hitr06/test2012/g1.dat')
	 KAPPALO=1250.; KAPPAHI=1750.;
 	 KAPPALO=1250.; KAPPAHI=1750.;
       ELSEIF (GID .EQ. 2) THEN
         OPEN (1,FILE='/asl/s1/motteler/read_hitr06/test2012/g2.dat')
         KAPPALO=0500.; KAPPAHI=0750.;
         KAPPALO=0850.; KAPPAHI=1150.;
         KAPPALO=2150.; KAPPAHI=2450.;	 	 	 
       ELSEIF (GID .EQ. 3) THEN
         OPEN (1,FILE='/asl/s1/motteler/read_hitr06/test2012/g3.dat')
         KAPPALO=875.; KAPPAHI=1175.;	 	 
       ELSEIF (GID .EQ. 4) THEN
         OPEN (1,FILE='/asl/s1/motteler/read_hitr06/test2012/g4.dat')
         KAPPALO=2100.; KAPPAHI=2300.;	 	 	 
       ELSEIF (GID .EQ. 5) THEN
         OPEN (1,FILE='/asl/s1/motteler/read_hitr06/test2012/g5.dat')
         KAPPALO=1900.; KAPPAHI=2300.;	 	 	 
       ELSEIF (GID .EQ. 6) THEN
         OPEN (1,FILE='/asl/s1/motteler/read_hitr06/test2012/g6.dat')
         KAPPALO=1250.; KAPPAHI=1350.;	 	 	 
       ELSEIF (GID .EQ. 9) THEN
         OPEN (1,FILE='/asl/s1/motteler/read_hitr06/test2012/g9.dat')
         KAPPALO=1250.; KAPPAHI=1450.;	 	 	 
       ELSEIF (GID .EQ. 12) THEN
         OPEN (1,FILE='/asl/s1/motteler/read_hitr06/test2012/g12.dat')
         KAPPALO=1200.; KAPPAHI=1400.;	 	 	 	 
       ELSE
         print *,'invalid GID ',GID
         STOP
       END IF
       OPEN (3,FILE='CO2abs.dat')        
       NLINES=1 
       DO WHILE(.TRUE.) 
!           ia,ib unused 
         READ(1,2,END=6) ia,ib,(VALX(J),J=1,NVALUES)
         IF ((ia .EQ. GID).and.(ib .eq. 1).and.(valx(1).GE.KAPPALO) 
     &                    .and. (valx(1) .LE. KAPPAHI)) THEN
           WRITE(3,*)ia,' ',ib,(VALX(J),J=1,NVALUES)
	   DO J=1,NVALUES
	     VAL(NLINES,J) = VALX(J)
	   END DO
           NLINES=NLINES+1
         END IF
       END DO 
 6     CONTINUE
 2     FORMAT(I2,I1,F12.6,2E10.3,2F5.4) 
 
!      #2: output file 
       OPEN (2,FILE='CO2abs.csv') 
!      Integrate each line from kappa-deltakappa to kappa+deltakappa 
!      Divide the result by 2*deltakappa to get the average intensity. 
!      Assign result to the midway value [kappa-deltakappa,kappa+deltakappa]
!      ie to kappa
       DO 10 KAPPA=KAPPALO+DELTAKAPPA,KAPPAHI+DELTAKAPPA,2*DELTAKAPPA 
         S=0. 
         DO 20 LINE=1,NLINES 
           V=   ILNZ(KAPPA+DELTAKAPPA,VAL(LINE,1),VAL(LINE,4)) 
           V= V-ILNZ(KAPPA-DELTAKAPPA,VAL(LINE,1),VAL(LINE,4))
!	   IF (ABS(KAPPA - (KAPPAHI+KAPPALO)/2.0) .LE. 0.1) print *,V
           S=S+VAL(LINE,2)*V 
 20      CONTINUE 
       S=S/(2*DELTAKAPPA) 
       WRITE(2,'(E16.6,",",E16.6)')KAPPA,S 
 10   CONTINUE 
      WRITE(6,*) NLINES ,' Lines',(KAPPAHI-KAPPALO)/2/DELTAKAPPA

      STOP 
      END

c************************************************************************
!     ILNZ returns the indefinite integral of the Lorentz line shape function 
      REAL*8 FUNCTION ILNZ(K,K0,G) 
      REAL*8 K,K0,G 
      ILNZ=DATAN2((K-K0),G)/3.141592658979D0 
      END

      
