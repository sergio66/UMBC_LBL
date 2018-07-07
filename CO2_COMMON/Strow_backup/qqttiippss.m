function [a1,b1,c1,d1]=qqttiippss(isotope)
%this produces the qtips coeffs from qtips.pas in the HITRAN distr

if (isotope == 1) 
        a1=-0.2199475485E+01; b1=0.9675055715E+00; 
        c1=-0.8082711378E-03; d1=0.2803987451E-05; 
elseif (isotope == 2) 
        a1=-.38840E+01; b1=.19263E+01; 
        c1=-.16058E-02; d1=.58202E-05; 
elseif (isotope == 3) 
        a1=-.47289E+01; b1=.20527E+01; 
        c1=-.17421E-02; d1=.60748E-05; 
  end 
  
