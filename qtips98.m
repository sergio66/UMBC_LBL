function [a,b,c,d,g] = qtips(gasid,liso);

% function [a,b,c,d,g] = qtips(gasid);
%
% returns nuclear degeneracy factors g and Gamache's internal partition
% sum coefficients a, c, c, and d.  This is taken from Dave Edwards 
% qtips.f.
%

%C***********************************************************************
%C
%C  PROGRAM        QTIPS   BLOCK DATA
%C
%C  PURPOSE        TOTAL INTERNAL PARTITION SUMS
%C
%C  VERSION        3.0   D.P. EDWARDS   11/01/90 
%C
%C  DESCRIPTION    COEFFICIENTS FOR CALCULATING THE TOTAL INTERNAL 
%C                 PARTITION SUMS
%C  
%C                 THIS IS A COPY OF THE ROUTINE QSUMS
%C                 FROM THE PROGRAM TIPS BY R.R.GAMACHE WHICH CALCULATES
%C                 THE TOTAL INTERNAL PARTITION FUNCTIONS IN THE 
%C                 HITRAN SELECT ROUTINE
%C                 ...LAST MODIFIED MAY 24, 1991
%C
%C***********************************************************************
%C
%C...STATE INDEPENDENT DEGENERACY FACTORS: GJ 
%C...(INCLUDES GENERAL NUCLEAR FACTOR (P(2I+1)), (2S+1), AND (2-DL0)
%
%      DATA GJ/
%

% this is from /salsify/packages/HITRAN/HITRAN96/SOFTWARE/GENERIC/BD_QT.FOR
% H98 : more uptodate version to to log onto asl.umbc.edu and go to 
% /asl/data/hitran/HITRAN98/SOFTWARE/GENERIC/BD_QT.FOR
% with irange == 1 (70 <= T <= 500)

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
%% H96if gasid == 19; g = [ 1  1  2  1 ];	 end	% OCS ***iso missing!**
if gasid == 19; g = [ 1  1  2  1 ];		 end	% OCS  
if gasid == 20; g = [ 1  2  1 ];		 end	% H2CO
if gasid == 21; g = [ 8  8 ];			 end	% HOCL
if gasid == 22; g = [ 1 ];			 end	% N2 
if gasid == 23; g = [ 6  12  4 ];		 end	% HCN
if gasid == 24; g = [ 4  4 ];			 end	% CH3CL
%% H96if gasid == 25; g = [ 1    ];		 end  	% H2O2 ***
if gasid == 25; g = [ 4	];			 end  	% H2O2
if gasid == 26; g = [ 1  8 ];			 end	% C2H2
if gasid == 27; g = [ 64 ];			 end	% C2H6
if gasid == 28; g = [ 2	];			 end	% PH3
if gasid == 29; g = [ 1	];			 end 	% COF2
if gasid == 30; g = [ 1	];			 end	% SF6
%% H96if gasid == 31; g = [ 1 1 4];        	 end	% H2S ***
if gasid == 31; g = [ 1 1 1];	        	 end	% H2S
%%H96if gasid == 32; g = [ 4	];		 end	% HCOOH ***
%%if gasid == 32; g = [ 1 1 1 1 1 1 1 1 1];	 end	%HCOOH *isos missing!**
if gasid == 32; g = [ 4	];		         end	% HCOOH 

g = g';

if (length(g) ~= liso)
  fprintf(1,'there is inconsistency in number of isotopes \n');
  fprintf(1,'in qtips.m %3i and that in mass.dat %3i \n',length(g),liso);
  error('please check!!!!!!')
  end

%C...TOTAL INTERNAL PARTITION SUMS FOR 70 - 405 K RANGE

if gasid == 1  		% H20 isotopes 161  181  171  162
	abcd = [ ...
  -.37688E+01, .26168E+00, .13497E-02,-.66013E-06 ;
  -.38381E+01, .26466E+00, .13555E-02,-.65372E-06 ; 
  -.22842E+02, .15840E+01, .81575E-02,-.39650E-05 ; 
  -.20481E+02, .13017E+01, .66225E-02,-.30447E-05];
 
fprintf(1,'The isotopes are  161,181,171,162 \n');
%which_isotope=input('enter which isotope (0 for all)');
end

if gasid == 2      %CO2 sotopes 626 636 628 627 638 637 828 728
	abcd = [ ...
   -.21995E+01, .96751E+00,    -.80827E-03, .28040E-05 ;
   -.38840E+01, .19263E+01,    -.16058E-02, .58202E-05 ; 
   -.47289E+01, .20527E+01,    -.17421E-02, .60748E-05 ; 
   -.27475E+02, .11973E+02,    -.10110E-01, .35187E-04 ; 
   -.84191E+01, .41186E+01,    -.34961E-02, .12750E-04 ; 
   -.48468E+02, .23838E+02,    -.20089E-01, .73067E-04 ; 
   -.22278E+01, .10840E+01,    -.89718E-03, .32143E-05 ; 
   -.29547E+02, .12714E+02,    -.10913E-01, .38169E-04]; 
fprintf(1,'The isotopes are   626,636,628,627,638,637,828,728 \n');
%which_isotope=input('enter which isotope (0 for all)');
end

if gasid == 3     %O3 isotopes = 666 668 686 667 676
%originally only had 3 isotopes  ---H96 had 5 .. then H98 is messed up
%and inconsistent -- bd_qt.for files claims there are 3 isotopes, others claim
%there are 5 isotopes. so i've updated the data for the first 3 isotopes
	abcd = [ ...
-.13459E+03, .62255E+01,       .14811E-01, .18608E-04 ;
-.24722E+03, .12331E+02,       .38336E-01, .26446E-04 ;
-.12359E+03, .60957E+01,       .18239E-01, .13939E-04 ;
       -.20540E+04, .85998E+02,  .12667E+00, .33026E-03; 
       -.10148E+04, .42494E+02,  .62586E-01, .16319E-03]; 
fprintf(1,'The isotopes are  666,668,686,667,676 \n');
%%%%%%%%fprintf(1,'The isotopes are  666,668,686 \n')  H92;
%which_isotope=input('enter which isotope (0 for all)');
 end

if gasid == 4        %N2O isotopes 446 456 546 448 447
	abcd = [ ...
-.95291E+01, .15719E+02,       -.12063E-01, .53781E-04 ;
 .48994E+01, .10211E+02,       -.62964E-02, .33355E-04 ; 
-.28797E+01, .10763E+02,       -.78058E-02, .36321E-04 ; 
 .25668E+02, .15803E+02,       -.67882E-02, .44093E-04 ; 
 .18836E+03, .91152E+02,       -.31071E-01, .23789E-03]; 
fprintf(1,'The isotopes are   446,456,546,448,447 \n');
%which_isotope=input('enter which isotope (0 for all)');
end

if gasid == 5      %CO isotopes 26 36 28 27 38 37    
%originally had 5 isotopes  --- H96 has 6 isotopes
%H98 has 5 ... what a mess!!!!! so last line is from H96 database
	abcd = [ ...
 .31591E+00, .36205E+00,  -.22603E-05, .61215E-08 ;
 .62120E+00, .75758E+00,  -.59190E-05, .15232E-07 ; 
 .30985E+00, .38025E+00,  -.29998E-05, .76646E-08 ; 
 .18757E+01, .22289E+01,  -.15793E-04, .41607E-07 ; 
 .60693E+00, .79754E+00,  -.78021E-05, .19200E-07; 
           .32731E+01, .46577E+01,   -.69833E-04, .18853E-06];
fprintf(1,'The isotopes are   26,36,28,27,38,37 \n');
%%%%fprintf(1,'The isotopes are   26,36,28,27,38 \n')   H92;
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 6 %CH4 isoptopes 211 311 212
	abcd = [ ...
-.17475E+02, .95375E+00,     .39758E-02, -.81837E-06 ;
-.27757E+02, .17264E+01,     .93304E-02, -.48181E-05 ; 
-.89810E+03, .44451E+02,     .17474E+00, -.22469E-04];
fprintf(1,'The isotopes are    211,311,212 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
 end

if gasid == 7  %O2 isotopes 66 68 67
	abcd = [ ...
-.10000E+01, .00000E+00,      .00000E+00, .00000E+00 ;
-.10000E+01, .00000E+00,      .00000E+00, .00000E+00 ; 
-.10000E+01, .00000E+00,      .00000E+00, .00000E+00];
fprintf(1,'The isotopes are     66,68,67 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 8        %NO isotopes = 46 56 48
	abcd = [ ...
-.44213E+02, .72098E+01,   .21853E-01,-.23036E-04 ;
-.15289E+02, .33260E+01,   .10040E-01,-.10562E-04 ;
-.46938E+02, .76070E+01,   .23010E-01,-.24207E-04];
fprintf(1,'The isotopes are   46,56,48 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
 end

if gasid == 9         %SO2 isotopes - 626 646
	abcd = [ ...
  -.17187E+03, .94104E+01,    .34620E-01, .25199E-04 ;
  -.17263E+03, .94528E+01,    .34777E-01, .25262E-04]; 
fprintf(1,'The isotopes are  626,646 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 10    %NO2 isotope   646
	abcd = [ ...
  -.44875E+03, .22359E+02,   .78905E-01, .21910E-04];
end

if gasid == 11       %NH3 isotope 4111 5111
	abcd = [ ...
  -.48197E+02, .27739E+01,   .11492E-01,-.18209E-05 ;
  -.32700E+02, .18444E+01,   .77001E-02,-.12388E-05];
fprintf(1,'The isotopes are  4111,5111 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 12       %HNO3 isotope 146
	abcd = [ ...
  -.74208E+04, .34984E+03,        .89051E-01, .39356E-02];
end

if gasid == 13            %OH isotope 61 81 62
	abcd = [ ...
  .19128E+02, .28443E+00,   .97670E-03,-.10688E-05 ;
  .19035E+02, .28770E+00,   .97945E-03,-.10718E-05 ; 
  .36233E+02, .11952E+01,   .38603E-02,-.64868E-05];
fprintf(1,'The isotopes are   61,81,62 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 14        %HF isotope = 19
	abcd = [ ...
 .15649E+01, .13318E+00,       .80622E-05,-.83354E-08];
end

if gasid == 15        %HCl isotopes 15 17
	abcd = [ ...
.28877E+01, .53077E+00,       .99904E-05,-.70856E-08 ;
.28873E+01, .53157E+00,       .99796E-05,-.70647E-08]; 
fprintf(1,'The isotopes are   15,17 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 16     %HBr isotopes 19 11
	abcd = [ ...
  .28329E+01, .66462E+00,       .83420E-05,-.30996E-08 ;
  .28329E+01, .66483E+00,       .83457E-05,-.31074E-08];
fprintf(1,'The isotopes are   19,11 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 17         %HI iso 17
	abcd = [ ...
  .41379E+01, .12977E+01,          .61598E-05, .10382E-07];
end

if gasid == 18         %ClO  iso 56 76
	abcd = [ ...
  .38740E+03, .28000E+02,      .48063E-01, .10208E-04 ;
  .39320E+03, .28483E+02,      .48795E-01, .10827E-04];
fprintf(1,'The isotopes are  56 76 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 19   %OCS 622 624 632 822
	abcd = [ ...
 .18600E+02, .31185E+01,        .30405E-03, .85400E-05 ;
 .19065E+02, .31965E+01,        .31228E-03, .87535E-05 ; 
 .42369E+02, .61394E+01,        .13090E-02, .16856E-04 ;
 .21643E+02, .32816E+01,        .57748E-03, .90034E-05];
fprintf(1,'The isotopes are  622,624,632,822  \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 20         %H2CO iso 126 136 128
	abcd = [ ...
 -.89326E+02, .46062E+01,    .19019E-01,-.33930E-05 ;
 -.18321E+03, .94446E+01,    .39010E-01,-.69664E-05 ;
 -.89326E+02, .46062E+01,    .19019E-01,-.33930E-05];
fprintf(1,'The isotopes are   126 136 128   \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 21     %isotopes HOCl 165 167
	abcd = [ ...
-.62547E+03, .31546E+02,       .11132E+00, .32438E-04 ;
-.60170E+03, .31312E+02,       .11841E+00, .23717E-04];
fprintf(1,'The isotopes are  165,167  \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 22     %isotopes N2 44
	abcd = [ ...
   .36774E+00, .39331E+00,     -.91410E-06, .34386E-08];

end

if gasid == 23          %isotopes HCN 124 134 125
	abcd = [ ...
 -.97107E+00, .29506E+01,    -.16077E-02, .61148E-05 ;
 -.16460E+01, .60490E+01,    -.32724E-02, .12632E-04 ;
 -.40184E+00, .20202E+01,    -.10855E-02, .42504E-05];
fprintf(1,'The isotopes are   124,134,125 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 24      %isotopes  CH3Cl  --   215 217
	abcd = [ ...
  -.89695E+03, .40155E+02,       .82775E-01, .13400E-03 ;
  -.91113E+03, .40791E+02,       .84091E-01, .13611E-03]; 
fprintf(1,'The isotopes are  215,217 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 25      %H2O2 iso = 1661 
	abcd = [ ...
-.95255E+03, .49483E+02,          .21249E+00,-.35489E-04];
end

if gasid == 26      % C2H2 iso = 1221 1231
	abcd = [ ...
  .25863E+01, .11921E+01,              -.79281E-03, .46225E-05 ;
  .20722E+02, .95361E+01,              -.63398E-02, .36976E-04];
fprintf(1,'The isotopes are  1221 1231 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end
 
if gasid == 27     %C2H6 iso = 1221
	abcd = [ ...
      -.10000E+01  .00000E+00   .00000E+00  .00000E+00];
end

if gasid == 28   % PH3  --  1111
	abcd = [ ...
     -.11388E+03, .69602E+01,          .17396E-01, .65088E-05];
end

if gasid == 29  %COF2  --   269
	abcd = [ ...
      -.10000E+01  .00000E+00   .00000E+00  .00000E+00];
end

if gasid == 30 %SF6  --    29
	abcd = [ ...
      -.10000E+01, .00000E+00,        .00000E+00, .00000E+00];
end

if gasid == 31      %H2S 121 141 131
%original 1 isotope
	abcd = [ ...
      -.10000E+01, .00000E+00,        .00000E+00, .00000E+00 ;
      -.10000E+01, .00000E+00,        .00000E+00, .00000E+00 ;
      -.10000E+01, .00000E+00,        .00000E+00, .00000E+00];
fprintf(1,'The isotopes are   121,141,131 \n');
%%%%%%%fprintf(1,'The isotopes are   121 \n');  H92
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 32   %HCOOH  --   126
	abcd = [ ...
      -.10000E+01, .00000E+00,        .00000E+00, .00000E+00];
end

a = abcd(:,1);
b = abcd(:,2);
c = abcd(:,3);
d = abcd(:,4);

