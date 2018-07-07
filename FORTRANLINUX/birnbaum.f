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
      parameter(npts=2001,ntemp=7,ndoc=10)

      integer         i,j,iTau,iTemp,nl,nh,nstp,IOUN
      real*8          tautau,fdif,fchi,df,chil,chiu,chif,chip
      real*8          vbc,vtc,dvc,temp(ntemp),doc(ndoc)

      real*8 chi100(npts,ndoc),chi150(npts,ndoc),chi200(npts,ndoc)
      real*8 chi250(npts,ndoc),chi300(npts,ndoc),chi350(npts,ndoc)
      real*8 chi400(npts,ndoc)
      real*8 chi111,chi211,chi121,chi221,chi112,chi212,chi122,chi222

c we have 10 d of c parameters tau2=1e-3*(1,3,5,...,19) (ndoc=10)
c we have 7 temperatures : 100,150,200,250,300,350,400  (ntemp=7)
c the frequency spacing of the y's is 0.25
c the dnu's range from -250 to +250 in steps of 0.25    (npts=2001)
      data temp/100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0/
      data doc/1.0, 4.0, 7.0, 10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0/
      data vbc,vtc,dvc/-250.0,+250.0,0.25/

c******************************* read in data section ***********************

      if (n_in .gt. MaxLen) THEN
        print *, 'in birnbaum.f n_in .gt. MaxLen'
        CALL MEXERRMSGTXT('POOEY')
        END IF

      if ((tau2 .lt. 1e-3) .or. (tau2 .gt. 28e-3)) then
        print *,'ooops!!! tau2 = ',tau2
        call mexErrMsgTxt('tau has to be between 1-28 exp-3')
        endif
      if ((T .lt. 100.0) .or. (T .gt. 400.0)) then
        print *,'ooops!!! t = ',T
        call mexErrMsgTxt('temp has to be between 100-400')
        endif

c lookup table produced using birn_lookupNEW.m, which is a script file
c previous version read in the text file      include 'birn_lookup.dat'

c now we load in bbiirrnn.dat 
      IOUN=10
c this uses the asymmetric functions that Dave Tobin used in his thesis
c see birn_lookupNEW.m
       open(unit=IOUN,FILE='../FORTRANLINUX/bbiirrnn.dat',STATUS='OLD',
     $                                   FORM='UNFORMATTED')

c this uses the symmetric functions that Scott Hannon thinks to use
c see birn_lookupNEW2.m
c      open(unit=IOUN,FILE='../FORTRANLINUX/bbiirrnn2.dat',STATUS='OLD',
c     $                                   FORM='UNFORMATTED')

      read(IOUN) dd
      if (abs(dd-temp(1)) .ge. 1e-2) then
        print *,'hmm : current temp should be 100, not ',dd
        call mexErrMsgTxt('pooey')
        end if
      do j=1,ndoc
        read(IOUN) (chi100(i,j),i=1,npts)
        end do

      read(IOUN) dd
      if (abs(dd-temp(2)) .ge. 1e-2) then
        print *,'hmm : current temp should be 150, not ',dd
        call mexErrMsgTxt('pooey')
        end if
      do j=1,ndoc
        read(IOUN) (chi150(i,j),i=1,npts)
        end do

      read(IOUN) dd
      if (abs(dd-temp(3)) .ge. 1e-2) then
        print *,'hmm : current temp should be 200, not ',dd
        call mexErrMsgTxt('pooey')
        end if
      do j=1,ndoc
        read(IOUN) (chi200(i,j),i=1,npts)
        end do

      read(IOUN) dd
      if (abs(dd-temp(4)) .ge. 1e-2) then
        print *,'hmm : current temp should be 250, not ',dd
        call mexErrMsgTxt('pooey')
        end if
      do j=1,ndoc
        read(IOUN) (chi250(i,j),i=1,npts)
        end do

      read(IOUN) dd
      if (abs(dd-temp(5)) .ge. 1e-2) then
        print *,'hmm : current temp should be 300, not ',dd
        call mexErrMsgTxt('pooey')
        end if
      do j=1,ndoc
        read(IOUN) (chi300(i,j),i=1,npts)
        end do

      read(IOUN) dd
      if (abs(dd-temp(6)) .ge. 1e-2) then
        print *,'hmm : current temp should be 350, not ',dd
        call mexErrMsgTxt('pooey')
        end if
      do j=1,ndoc
        read(IOUN) (chi350(i,j),i=1,npts)
        end do

      read(IOUN) dd
      if (abs(dd-temp(7)) .ge. 1e-2) then
        print *,'hmm : current temp should be 400, not ',dd
        call mexErrMsgTxt('pooey')
        end if
      do j=1,ndoc
        read(IOUN) (chi400(i,j),i=1,npts)
        end do

      close(IOUN)
c******************************* end of data section ***********************

      nl=1
      nh=n_in
      nstp=1

c find which temps (100,150,200,...,400) bracket the needed temp
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
        print *,'sorry : temp is outside 100:150:...:400'
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

        if (iTemp .eq. 1) then 
          chi111=chi100(j,iTau)
          chi211=chi100(j+1,iTau)
          chi121=chi100(j,iTau+1)
          chi221=chi100(j+1,iTau+1)
          chi112=chi150(j,iTau)
          chi212=chi150(j+1,iTau)
          chi122=chi150(j,iTau+1)
          chi222=chi150(j+1,iTau+1)
        elseif (iTemp .eq. 2) then 
          chi111=chi150(j,iTau)
          chi211=chi150(j+1,iTau)
          chi121=chi150(j,iTau+1)
          chi221=chi150(j+1,iTau+1)
          chi112=chi200(j,iTau)
          chi212=chi200(j+1,iTau)
          chi122=chi200(j,iTau+1)
          chi222=chi200(j+1,iTau+1)
        elseif (iTemp .eq. 3) then 
          chi111=chi200(j,iTau)
          chi211=chi200(j+1,iTau)
          chi121=chi200(j,iTau+1)
          chi221=chi200(j+1,iTau+1)
          chi112=chi250(j,iTau)
          chi212=chi250(j+1,iTau)
          chi122=chi250(j,iTau+1)
          chi222=chi250(j+1,iTau+1)
        elseif (iTemp .eq. 4) then 
          chi111=chi250(j,iTau)
          chi211=chi250(j+1,iTau)
          chi121=chi250(j,iTau+1)
          chi221=chi250(j+1,iTau+1)
          chi112=chi300(j,iTau)
          chi212=chi300(j+1,iTau)
          chi122=chi300(j,iTau+1)
          chi222=chi300(j+1,iTau+1)
        elseif (iTemp .eq. 5) then 
          chi111=chi300(j,iTau)
          chi211=chi300(j+1,iTau)
          chi121=chi300(j,iTau+1)
          chi221=chi300(j+1,iTau+1)
          chi112=chi350(j,iTau)
          chi212=chi350(j+1,iTau)
          chi122=chi350(j,iTau+1)
          chi222=chi350(j+1,iTau+1)
        elseif (iTemp .eq. 6) then 
          chi111=chi350(j,iTau)
          chi211=chi350(j+1,iTau)
          chi121=chi350(j,iTau+1)
          chi221=chi350(j+1,iTau+1)
          chi112=chi400(j,iTau)
          chi212=chi400(j+1,iTau)
          chi122=chi400(j,iTau+1)
          chi222=chi400(j+1,iTau+1)
          endif

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

