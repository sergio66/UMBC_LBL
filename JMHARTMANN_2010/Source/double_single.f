      real*4 ra,rb
      real*8 da,db

      print *,'enter ra,rb : '
      read *,ra,rb
      da = dble(ra)
      db = dble(rb)
      print *,'real input, converted to double : '
      write(*,123) ra,rb
      write(*,223) da,db
      print *,' '

      print *,'enter da,db : '
      read *,da,db
      ra = sngl(da)
      rb = sngl(db)
      print *,'double input, converted to single : '
      write(*,123) ra,rb
      write(*,223) da,db
      print *,' '

      ra = 3.1415e0
      rb = 06.0e-4
      ra = anint(605.00001)
      rb = anint(2829.99999)
      da = dble(ra)
      db = dble(rb)
      print *,'set single, converted to double : '
      write(*,123) ra,rb
      write(*,223) da,db
      print *,' '

      da = 3.1415d0
      db = 06.0D-4
      da = dnint(605.00001d0)
      db = dnint(2829.99999d0)
      ra = sngl(da)
      rb = sngl(db)
      print *,'set double, converted to single : '
      write(*,123) ra,rb
      write(*,223) da,db
      print *,' '

 123  Format(2(1pe26.18))
 224  Format(2(1pd26.18))
 223  Format(d26.18,d28.18)
      stop
      end
