function [which_isotope] = iso_list(gasid);
%this is from the latest 
%/salsify/packages/HITRAN/HITRAN96/SOFTWARE/GENERIC/BD_QT.FOR

if gasid ==  1; g = [ 1  1  6  6 ]; 		 end  	% H2O
if gasid ==  2; g = [ 1  2  1  6  2  12  1  6 ] ;end   	% CO2
if gasid ==  3; g = [ 1  1  1 6  6]; 		 end 	% O3
if gasid ==  4; g = [ 9  6  6  9 54 ]; 		 end	% N2O
if gasid ==  5; g = [ 1  2  1  6 2  12];         end	% CO
if gasid ==  6; g = [ 1  2  3 ]; 		 end	% CH4
if gasid ==  7; g = [ 1  1  6 ];		 end	% O2  
if gasid ==  8; g = [ 3 2 3  ]; 		 end	% NO
if gasid ==  9; g = [ 1  1 ]; 			 end  	% SO2
if gasid == 10; g = [ 3 ];			 end 	% NO2
if gasid == 11; g = [ 3  2 ];			 end	% NH3
if gasid == 12; g = [ 6 ];			 end	% HNO3
if gasid == 13; g = [ 2 2 3 ];  		 end	% OH
if gasid == 14; g = [ 4 ];			 end	% HF
if gasid == 15; g = [ 8  8 ];  			 end	% HCL
if gasid == 16; g = [ 8  8 ]; 			 end	% HBR
if gasid == 17; g = [ 12 ];			 end	% HI
if gasid == 18; g = [ 4  4 ]; 			 end	% CLO
if gasid == 19; g = [ 1  1  2  1 ];		 end	% OCS
if gasid == 20; g = [ 1  2  1 ];		 end	% H2CO
if gasid == 21; g = [ 8  8 ];			 end	% HOCL
if gasid == 22; g = [ 1 ];			 end	% N2 
if gasid == 23; g = [ 6  12  4 ];		 end	% HCN
if gasid == 24; g = [ 4  4 ];			 end	% CH3CL
if gasid == 25; g = [ 1	];			 end  	% H2O2
if gasid == 26; g = [ 1  8 ];			 end	% C2H2
if gasid == 27; g = [ 64 ];			 end	% C2H6
if gasid == 28; g = [ 2	];			 end	% PH3
if gasid == 29; g = [ 1	];			 end 	% COF2
if gasid == 30; g = [ 1	];			 end	% SF6
if gasid == 31; g = [ 1 1 4];	        	 end	% H2S
if gasid == 32; g = [ 4	];			 end	% HCOOH

if gasid == 1  		% H20 isotopes 161  181  171  162
  fprintf(1,'The isotopes are  161,181,171,162 \n');
  which_isotope=input('enter which isotope (0 for all)');

elseif gasid == 2      %CO2 sotopes 626 636 628 627 638 637 828 728
  fprintf(1,'The isotopes are   626,636,628,627,638,637,828,728 \n');
  which_isotope=input('enter which isotope (0 for all)');

elseif gasid == 3     %O3 isotopes = 666 668 686 667 676
  fprintf(1,'The isotopes are  666,668,686,667,676 \n');
  %%%%%%%%fprintf(1,'The isotopes are  666,668,686 \n')  H92;
  which_isotope=input('enter which isotope (0 for all)');

elseif gasid == 4        %N2O isotopes 446 456 546 448 447
  fprintf(1,'The isotopes are   446,456,546,448,447 \n');
  which_isotope=input('enter which isotope (0 for all)');

elseif gasid == 5      %CO isotopes 26 36 28 27 38 37    
  fprintf(1,'The isotopes are   26,36,28,27,38,37 \n');
  %%%%fprintf(1,'The isotopes are   26,36,28,27,38 \n')   H92;
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 6 %CH4 isoptopes 211 311 212
  fprintf(1,'The isotopes are    211,311,212 \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 7  %O2 isotopes 66 68 67
  fprintf(1,'The isotopes are     66,68,67 \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 8        %NO isotopes = 46 56 48
  fprintf(1,'The isotopes are   46,56,48 \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 9         %SO2 isotopes - 626 646
  fprintf(1,'The isotopes are  626,646 \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 10    %NO2 isotope   646
  which_isotope=0;

elseif gasid == 11       %NH3 isotope 4111 5111
  fprintf(1,'The isotopes are  4111,5111 \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 12       %HNO3 isotope 146
  which_isotope=0;

elseif gasid == 13            %OH isotope 61 81 62
  fprintf(1,'The isotopes are   61,81,62 \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 14        %HF isotope = 19
  which_isotope=0;

elseif gasid == 15        %HCl isotopes 15 17
  fprintf(1,'The isotopes are   15,17 \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 16     %HBr isotopes 19 11
  fprintf(1,'The isotopes are   19,11 \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 17         %HI iso 17
  which_isotope=0;

elseif gasid == 18         %ClO  iso 56 76
  fprintf(1,'The isotopes are  56 76 \n');
  which_isotope=input('enter which isotope (0 for all)'); 

%%%note BD_ISO claims there are 5 isotopes, as does the actual database
%%%but BD_ABUN does not give any evidence of isotope 623!!!!!!
%%%neither does BD_QT!!! strange!!!!!!!!!!!
elseif gasid == 19   %OCS 622 624 632 822
  fprintf(1,'The isotopes are  622,624,632,822  \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 20         %H2CO iso 126 136 128
  fprintf(1,'The isotopes are   126 136 128   \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 21     %isotopes HOCl 165 167
  fprintf(1,'The isotopes are  165,167  \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 22     %isotopes N2 44
  which_isotope=0;

elseif gasid == 23          %isotopes HCN 124 134 125
  fprintf(1,'The isotopes are   124,134,125 \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 24      %isotopes  CH3Cl  --   215 217
  fprintf(1,'The isotopes are  215,217 \n');
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 25      %H2O2 iso = 1661 
  which_isotope=0;

elseif gasid == 26      % C2H2 iso = 1221 1231
  fprintf(1,'The isotopes are  1221 1231 \n');
  which_isotope=input('enter which isotope (0 for all)'); 
 
elseif gasid == 27     %C2H6 iso = 1221
  which_isotope=0;

elseif gasid == 28   % PH3  --  1111
  which_isotope=0;

elseif gasid == 29  %COF2  --   269
  which_isotope=0;

elseif gasid == 30 %SF6  --    29
  which_isotope=0;

elseif gasid == 31      %H2S 121 141 131
%original 1 isotope
  fprintf(1,'The isotopes are   121,141,131 \n');
  %%%%%%%fprintf(1,'The isotopes are   121 \n');  H92
  which_isotope=input('enter which isotope (0 for all)'); 

elseif gasid == 32   %HCOOH  --   126
  which_isotope=0;

  end










