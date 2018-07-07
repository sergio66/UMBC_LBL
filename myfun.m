        function F = myfun(x,c)
        F = [ 2*x(1) - x(2) - exp(c*x(1))
              -x(1) + 2*x(2) - exp(c*x(2))];
