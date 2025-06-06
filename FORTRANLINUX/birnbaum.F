      subroutine birnbaum(z,v,v0,w_tot,T,tau2,n_in)
c subroutine cousin(z,v,v0,w_tot,T,tau2,n_in)
c this is the birnbaum lineshape using lookup tables
c z    = results array
c v    = frequency array
c v0   = center freq
c T    = temperature
c tau2 = duration of collision
c n_in = no of input points

c see birn_lookup.m
c since w_tot is immaterial
c this does a scan of tau2 for 7 temps
c for ii=1:10
c   pow=1;                                                                 
c   doc=1e-3*(1+(ii-1)*3);                                                 
c   y100(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),100,doc);
c   y150(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),150,doc);
c   y200(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),200,doc);      
c   y250(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),250,doc);
c   y300(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),300,doc);
c   y350(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),350,doc);
c   y400(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),400,doc);
c   end

      include 'max.inc'

      real*8 z(MaxLen),v(MaxLen),v0,T,tau2,w_tot,dd
      integer n_in

c **********************************************************************
      integer npts,ndoc,ntemp
!      parameter(npts=2001,ntemp=7,ndoc=10)
      parameter(npts=2001,ntemp=29,ndoc=10)

      integer         i,j,iTau,iTemp,nl,nh,nstp,IOUN,iRead
      real*8          tautau,fdif,fchi,df,chil,chiu,chif,chip
      real*8          vbc,vtc,dvc,temp(ntemp),doc(ndoc)

      real*8 Tmax

c      real*8 chi100(npts,ndoc),chi150(npts,ndoc),chi200(npts,ndoc)
c      real*8 chi250(npts,ndoc),chi300(npts,ndoc),chi350(npts,ndoc)
c      real*8 chi400(npts,ndoc),chi450(npts.ndoc),chi500(npts,ndoc)
c      real*8 chi550(npts,ndoc),chi600(npts,ndoc),chi650(npts,ndoc)
c      real*8 chi700(npts,ndoc),chi750(npts,ndoc),chi800(npts,ndoc)
c      real*8 chi850(npts,ndoc),chi900(npts,ndoc),chi950(npts,ndoc)
c      real*8 chi1000(npts,ndoc),chi1050(npts.ndoc),chi1100(npts,ndoc)
c      real*8 chi1150(npts,ndoc),chi1200(npts,ndoc),chi1250(npts,ndoc)
c      real*8 chi1300(npts,ndoc),chi1350(npts.ndoc),chi1400(npts,ndoc)
c      real*8 chi1450(npts,ndoc),chi1500(npts,ndoc)
      real*8 chiT(ntemp,npts,ndoc)
      real*8 chi111,chi211,chi121,chi221,chi112,chi212,chi122,chi222

c we have 10 d of c parameters tau2=1e-3*(1,3,5,...,19) (ndoc=10)
c we have 07 temperatures : 100,150,200,250,300,350,400  (ntemp=7)
c we have 29 temperatures : 100:50:1500                  (ntemp=29)
c the frequency spacing of the y's is 0.25
c the dnu's range from -250 to +250 in steps of 0.25    (npts=2001)
      data temp/ 100.0,  150.0,  200.0,  250.0,  300.0,  350.0,  400.0,
     +           450.0,  500.0,  550.0,  600.0,  650.0,  700.0,  750.0, 
     +           800.0,
     +           850.0,  900.0,  950.0, 1000.0, 1050.0, 1100.0, 1150.0, 
     +          1200.0, 1250.0, 1300.0, 1350.0, 1400.0, 1450.0, 1500.0/

      data doc/1.0, 4.0, 7.0, 10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0/
      data vbc,vtc,dvc/-250.0,+250.0,0.25/

c******************************* read in data section ***********************

c orig Tmax = 400, now changed to 1500 
      Tmax = 400.0 
      Tmax = 1500.0
      if (n_in .gt. MaxLen) THEN
        print *, 'in birnbaum.f n_in .gt. MaxLen'
        CALL MEXERRMSGTXT('POOEY')
      END IF

      if ((tau2 .lt. 1e-3) .or. (tau2 .gt. 28e-3)) then
        print *,'ooops!!! tau2 = ',tau2
        call mexErrMsgTxt('tau has to be between 1-28 exp-3')
      endif
      if ((T .lt. 100.0) .or. (T .gt. Tmax)) then  
c used to be max 400 K, now changed to max 1500 K        
        print *,'ooops!!! t = ',T
        call mexErrMsgTxt('temp has to be between 100-1500')
      endif

c lookup table produced using birn_lookupNEW.m, which is a script file
c previous version read in the text file      include 'birn_lookup.dat'

c now we load in bbiirrnn.dat 
      IOUN=10

c this was prior to SOMETIME
c this uses the symmetric functions that Scott Hannon thinks to use
c see birn_lookupNEW2.m
c      open(unit=IOUN,FILE='../FORTRANLINUX/bbiirrnn2.dat',STATUS='OLD',
c     $                                   FORM='UNFORMATTED')

c this was prior to Feb 2020
c this uses the asymmetric functions that Dave Tobin used in his thesis
c see birn_lookupNEW.m
c       open(unit=IOUN,FILE='../FORTRANLINUX/bbiirrnn.dat',STATUS='OLD',
c     $                                   FORM='UNFORMATTED')

c this is after Feb 2020
c this uses the asymmetric functions that Dave Tobin used in his thesis
c see birn_lookupNEW.m
      iRead = 0
      open(unit=IOUN,FILE='../FORTRANLINUX/bbiirrnn_HITEMP.dat',
     $                      STATUS='OLD',FORM='UNFORMATTED')

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi100(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi150(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi200(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi250(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi300(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi350(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi400(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi450(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi500(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi550(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi600(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi650(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi700(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi750(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi800(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi850(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi900(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi950(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi1000(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi1050(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi1100(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi1150(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi1200(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi1250(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi1300(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi1350(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi1400(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi1450(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      iRead = iRead + 1
      read(IOUN) dd
      if (abs(dd-temp(iRead)) .ge. 1e-2) then
        print *,'hmm : current temp ',dd, ' should be ',temp(iRead)
        call mexErrMsgTxt('pooey')
      end if
      do j=1,ndoc
c        read(IOUN) (chi1500(i,j),i=1,npts)
        read(IOUN) (chiT(iRead,i,j),i=1,npts)
      end do

      close(IOUN)

      if (iRead .ne. ntemp) then
        print *,'hmm : should have read ',ntemp,' Temps, not ',iread
        call mexErrMsgTxt('iRead NE NTEMP')
      end if 
c******************************* end of data section ***********************

      nl=1
      nh=n_in
      nstp=1

c find which temps (100,150,200,...,max(temp)) bracket the needed temp
      j=-1
      tautau=T
      iTemp=1
      i=1
 10   continue
      if ((j .lt. 0) .and. (tautau .le. temp(i))) then
        j=1
        iTemp=i-1
      elseif (i .lt. ntemp) then
        i=i+1
        goto 10
      else
        print *,'sorry : temp is outside 100:150:...:1500'
        call mexErrMsgTxt('pooey')
        end if
      if (iTemp .lt. 1) iTemp = 1
      
c find which doc's (1,3,5,7,9,11,...,19) bracket the needed doc
      j=-1
      tautau=tau2/1e-3
      iTau=1
      i=1
 20   continue
      if ((j .lt. 0) .and. (tautau .le. doc(i))) then
        j=1
        iTau=i-1
      elseif (i .lt. ndoc) then
        i=i+1
        goto 20
      else
c        print *,'sorry : doc is outside 1:3:5:7:9:11:..:19 e-3'
        print *,'sorry : doc is outside 1:4:7:..:28 e-3'
        call mexErrMsgTxt('pooey')
        end if
      if (iTau .lt. 1) iTau = 1

c      print *,'n_in = ',n_in
c      print *,'temp = ',T,'bracketed by ',iTemp,iTemp+1  
c      print *,'tau = ',tautau,'bracketed by ',iTau,iTau+1  

c chi(npts,5,7) === chi(freq,doc,tempr)

c      do i=1600,1610
c        print *,i,chi100(i,1)
c        end do

      do i=nl,nh,nstp 
        fdif = (v(i) - v0) 
        j = min((npts-2),int(fdif/dvc)) + 1 
        j = j + 1000   !because we have chi100(1) ==> dnu = -250
        fchi = vbc + (j - 1)*dvc 
        df=(fdif-fchi)/dvc

        chi111=chiT(iTemp,j,iTau)
        chi211=chiT(iTemp,j+1,iTau)
        chi121=chiT(iTemp,j,iTau+1)
        chi221=chiT(iTemp,j+1,iTau+1)
        chi112=chiT(iTemp+1,j,iTau)
        chi212=chiT(iTemp+1,j+1,iTau)
        chi122=chiT(iTemp+1,j,iTau+1)
        chi222=chiT(iTemp+1,j+1,iTau+1)

c       do the interpolation in Tau, at temp iTemp
        chil=chi111 + (chi211-chi111)*df
        chiu=chi121 + (chi221-chi121)*df
        chif = chil + 
     +         (chiu-chil)*(tautau-doc(iTau))/(doc(iTau+1)-doc(iTau)) 

c       do the interpolation in Tau, at temp iTemp+1
        chil=chi112 + (chi212-chi112)*df
        chiu=chi122 + (chi222-chi122)*df
        chip = chil + 
     +         (chiu-chil)*(tautau-doc(iTau))/(doc(iTau+1)-doc(iTau)) 
 
c        do the interpolation in temp
         z(i)=chif + 
     +        (chip-chif)*(T-temp(iTemp))/(temp(iTemp+1)-temp(iTemp)) 

         enddo 

       return
       end

