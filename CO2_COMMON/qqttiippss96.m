function [a1,b1,c1,d1]=qqttiippss(isotope)
%this produces the qtips coeffs from qtips.pas in the HITRAN distr

if (isotope == 1) 
        a1=-0.2199475485E+01; b1=0.9675055715E+00; 
        c1=-0.8082711378E-03; d1=0.2803987451E-05; 
elseif (isotope == 2) 
        a1=-.20631E+01; b1=.18873E+01; 
        c1=-.13669E-02; d1=.54032E-05; 
elseif (isotope == 3) 
        a1=-.29175E+01; b1=.20114E+01; 
        c1=-.14786E-02; d1=.55941E-05; 
  end 
  
