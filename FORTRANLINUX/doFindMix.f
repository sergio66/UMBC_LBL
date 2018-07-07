      subroutine doFindMix(Y_1st,trans_ampl,W_matrix,freqq,n_in) 
 
      include 'max.inc' 

      integer n_in
      real*8 Y_1st(MaxPQR),freqq(MaxPQR),trans_ampl(MaxPQR)
      real*8 W_matrix(MaxPQR,MaxPQR)

      integer m,n,nnn

      if (MaxPQR .lt. n_in) then
        print *,'in doFindMix need MaxPQR > n_in'
        stop
        end if

      nnn=int(n_in)

c      print *,'n_in,nnn = ',n_in,nnn

      do n=1,nnn
        Y_1st(n)=0.0
        end do

      do n=1,nnn
        do m=1,nnn
          if (m .ne. n) then
            Y_1st(n)=Y_1st(n)+trans_ampl(m)*(W_matrix(m,n))/
     $                                    (freqq(n)-freqq(m))
            endif 
          enddo
        enddo 

      do n=1,nnn
        Y_1st(n)=2*Y_1st(n)/trans_ampl(n)
        end do

      return
      end
 
