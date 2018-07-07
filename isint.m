function y = isint(x)
%ISINT - Is Integer
%function returns 1 if the input is an integer and 0 otherwise.

if mod(x,1) == 0
   y = 1;
else
   y = 0;
   end

