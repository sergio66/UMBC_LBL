function [a,b,c,d,g] = qtips(gasid,liso);

% function [a,b,c,d,g] = qtips(gasid,liso);
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

% which_isotope enables the user to choose ALL isotopes (0) or one of
% the particular 1..N isotopes of gasID

% this is from HITRAN92 qtips.m
% with irange == 1 (70 <= T <= 500)

if gasid ==  1; g = [ 1  1  6  6 ]; 		 end  	% H2O
if gasid ==  2; g = [ 1  2  1  6  2  12  1  6 ] ;end   	% CO2
if gasid ==  3; g = [ 1  1  1]; 		 end 	% O3
if gasid ==  4; g = [ 9  6  6  9 54 ]; 		 end	% N2O
if gasid ==  5; g = [ 1  2  1  6 2];             end	% CO
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
if gasid == 31; g = [ 1 ];			 end	% H2S
if gasid == 32; g = [ 4	];			 end	% HCOOH

g = g';

if (length(g) ~= liso)
  fprintf(1,'there is inconsistency in number of isotopes \n');
  fprintf(1,'in qtips.m %3i and that in mass.dat %3i \n',length(g),liso);
  error('please check!!!!!!')
  end

%C...TOTAL INTERNAL PARTITION SUMS FOR 70 - 405 K RANGE

which_isotope=0;
if gasid == 1  		% H20 isotopes 161  181  171  162
	abcd = [ ...
        -.44405E+01 .27678E+00  .12536E-02  -.48938E-06
        -.43624E+01 .27647E+00  .12802E-02  -.52046E-06
        -.25767E+02 .16458E+01  .76905E-02  -.31668E-05
        -.23916E+02 .13793E+01  .61246E-02  -.21530E-05];
fprintf(1,'The isotopes are  161,181,171,162 \n');
%which_isotope=input('enter which isotope (0 for all)');
end

if gasid == 2      %CO2 sotopes 626 636 628 627 638 637 828 728
	abcd = [ ...
       -.13617E+01, .94899E+00,  -.69259E-03, .25974E-05
       -.20631E+01, .18873E+01,  -.13669E-02, .54032E-05 
       -.29175E+01, .20114E+01,  -.14786E-02, .55941E-05 
       -.16558E+02, .11733E+02,  -.85844E-02, .32379E-04 
       -.44685E+01, .40330E+01,  -.29590E-02, .11770E-04
       -.26263E+02, .23350E+02,  -.17032E-01, .67532E-04 
       -.14811E+01, .10667E+01,  -.78758E-03, .30133E-05 
       -.17600E+02, .12445E+02,  -.91837E-02, .34915E-04];
fprintf(1,'The isotopes are   626,636,628,627,638,637,828,728 \n');
%which_isotope=input('enter which isotope (0 for all)');
end

if gasid == 3     %O3 isotopes = 666 668 686 667 676
%originally only had 3 isotopes
	abcd = [ ...
       -.16443E+03, .69047E+01,  .10396E-01, .26669E-04
       -.35222E+03, .14796E+02,  .21475E-01, .59891E-04 
       -.17466E+03, .72912E+01,  .10093E-01, .29991E-04]; 
fprintf(1,'The isotopes are  666,668,686 \n');
%which_isotope=input('enter which isotope (0 for all)');
 end

if gasid == 4        %N2O isotopes 446 456 546 448 447
	abcd = [ ...
        .24892E+02, .14979E+02,    -.76213E-02, .46310E-04
        .36318E+02, .95497E+01,    -.23943E-02, .26842E-04 
        .24241E+02, .10179E+02,    -.43002E-02, .30425E-04 
        .67708E+02, .14878E+02,    -.10730E-02, .34254E-04 
        .50069E+03, .84526E+02,     .83494E-02, .17154E-03]; 
fprintf(1,'The isotopes are   446,456,546,448,447 \n');
%which_isotope=input('enter which isotope (0 for all)');
end

if gasid == 5      %CO isotopes 26 36 28 27 38 37    
%originally had 5 isotopes
	abcd = [ ...
           .27758E+00, .36290E+00,   -.74669E-05, .14896E-07
           .53142E+00, .75953E+00,   -.17810E-04, .35160E-07 
           .26593E+00, .38126E+00,   -.92083E-05, .18086E-07
           .16376E+01, .22343E+01,   -.49025E-04, .97389E-07
           .51216E+00, .79978E+00,   -.21784E-04, .42749E-07];
fprintf(1,'The isotopes are   26,36,28,27,38 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 6 %CH4 isoptopes 211 322 212
	abcd = [ ...
         -.26479E+02, .11557E+01,     .26831E-02, .15117E-05
         -.52956E+02, .23113E+01,     .53659E-02, .30232E-05 
         -.21577E+03, .93318E+01,     .21779E-01, .12183E-04]; 
fprintf(1,'The isotopes are    211,311,212 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
 end

if gasid == 7  %O2 isotopes 66 68 67
	abcd = [ ...
            .35923E+00, .73534E+00,   -.64870E-04, .13073E-06
           -.40039E+01, .15595E+01,   -.15357E-03, .30969E-06 
           -.23325E+02, .90981E+01,   -.84435E-03, .17062E-05]; 
fprintf(1,'The isotopes are     66,68,67 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 8        %NO isotopes = 46 56 48
	abcd = [ ...
           -.75888E+02, .79048E+01,    .17555E-01, -.15606E-04
           -.29980E+02, .36479E+01,    .80522E-02, -.71296E-05 
           -.80558E+02, .83447E+01,    .18448E-01, -.16323E-04]; 
fprintf(1,'The isotopes are   46,56,48 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
 end

if gasid == 9         %SO2 isotopes - 626 646
	abcd = [ ...
           -.24056E+03, .11101E+02,   .22164E-01, .52334E-04
           -.24167E+03, .11151E+02,   .22270E-01, .52550E-04]; 
fprintf(1,'The isotopes are  626,646 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 10    %NO2 isotope   646
	abcd = [ ...
         -.53042E+03, .24216E+02,        .66856E-01, .43823E-04];
end

if gasid == 11       %NH3 isotope 4111 5111
	abcd = [ ...
            -.42037E+02, .25976E+01,  .13073E-01, -.62230E-05
            -.28609E+02, .17272E+01,  .87529E-02, -.41714E-05]; 
fprintf(1,'The isotopes are  4111,5111 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 12       %HNO3 isotope 146
	abcd = [ ...
           -.10000E+01, .00000E+00,   .00000E+00, .00000E+00];
end

if gasid == 13            %OH isotope 61 81 62
	abcd = [ ...
         .17478E+02, .31954E+00,  .76581E-03,-.71337E-06 
         .17354E+02, .32350E+00,  .76446E-03,-.70932E-06 
         .30717E+02, .13135E+01,  .31430E-02,-.28371E-05];
fprintf(1,'The isotopes are   61,81,62 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 14        %HF isotope = 19
	abcd = [ ...
            .15486E+01, .13350E+00,   .59154E-05,-.46889E-08]
end

if gasid == 15        %HCl isotopes 15 17
	abcd = [ ...
          .28627E+01, .53122E+00,     .67464E-05,-.16730E-08
          .28617E+01, .53203E+00,     .66553E-05,-.15168E-08]; 
fprintf(1,'The isotopes are   15,17 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 16     %HBr isotopes 19 11
	abcd = [ ...
           .27963E+01, .66532E+00,     .34255E-05, .52274E-08
           .27953E+01, .66554E+00,     .32931E-05, .54823E-08]; 
fprintf(1,'The isotopes are   19,11 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 17         %HI iso 17
	abcd = [ ...
         .40170E+01, .13003E+01,    -.11409E-04, .40026E-07];
end

if gasid == 18         %ClO  iso 57 76
	abcd = [ ...
        .36387E+03, .28367E+02,   .46556E-01, .12058E-04
        .37039E+03, .28834E+02,   .47392E-01, .12522E-04];
fprintf(1,'The isotopes are  57 76 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 19   %OCS 622 624 632 822
	abcd = [ ...
           -.93697E+00, .36090E+01,    -.34552E-02, .17462E-04
           -.11536E+01, .37028E+01,    -.35582E-02, .17922E-04 
           -.61015E+00, .72200E+01,    -.70044E-02, .36708E-04 
           -.21569E+00, .38332E+01,    -.36783E-02, .19177E-04]; 
fprintf(1,'The isotopes are  622,624,632,822  \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 20         %H2CO iso 126 136 128
	abcd = [ ...
        -.11760E+03, .46885E+01,     .15088E-01, .35367E-05
        -.24126E+03, .96134E+01,     .30938E-01, .72579E-05 
        -.11999E+03, .52912E+01,     .14686E-01, .43505E-05]; 
fprintf(1,'The isotopes are   126 136 128   \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 21     %isotopes HOCl 165 167
	abcd = [ ...
           -.73640E+03, .34149E+02,     .93554E-01, .67409E-04
           -.74923E+03, .34747E+02,     .95251E-01, .68523E-04];
fprintf(1,'The isotopes are  165,167  \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 22     %isotopes N2 44
	abcd = [ ...
          .13684E+01, .15756E+01,      -.18511E-04, .38960E-07];
end

if gasid == 23          %isotopes HCN 124 134 125
	abcd = [ ...
            -.13992E+01, .29619E+01,    -.17464E-02, .65937E-05
            -.25869E+01, .60744E+01,    -.35719E-02, .13654E-04
            -.11408E+01, .20353E+01,    -.12159E-02, .46375E-05]; 
fprintf(1,'The isotopes are   124,134,125 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 24      %isotopes  CH3Cl  --   215 217
	abcd = [ ...
        -.91416E+03, .34081E+02,      .75461E-02, .17933E-03
        -.92868E+03, .34621E+02,      .76674E-02, .18217E-03]; 
fprintf(1,'The isotopes are  215,217 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end

if gasid == 25      %H2O2 iso = 1661 
	abcd = [ ...
            -.36499E+03, .13712E+02,  .38658E-01, .23052E-04];
end

if gasid == 26      % C2H2 iso = 1221 1231
	abcd = [ ...
           -.83088E+01, .14484E+01,      -.25946E-02, .84612E-05
           -.66736E+02, .11592E+02,      -.20779E-01, .67719E-04]; 
fprintf(1,'The isotopes are  1221 1231 \n');
%which_isotope=input('enter which isotope (0 for all)'); 
end
 
if gasid == 27     %C2H6 iso = 1221
	abcd = [ ...
      -.10000E+01  .00000E+00   .00000E+00  .00000E+00];
end

if gasid == 28   % PH3  --  1111
	abcd = [ ...
        -.15068E+03, .64718E+01,        .12588E-01, .14759E-04];

end

if gasid == 29  %COF2  --   269
	abcd = [ ...
     -.54180E+04, .18868E+03,     -.33139E+00, .18650E-02];
end

if gasid == 30 %SF6  --    29
	abcd = [ ...
      -.10000E+01, .00000E+00,        .00000E+00, .00000E+00];
end

if gasid == 31      %H2S 121 141 131
%original 1 isotope
	abcd = [ ...
      -.15521E+02, .83130E+00,   .33656E-02,-.85691E-06]; 
end

if gasid == 32   %HCOOH  --   126
	abcd = [ ...
       -.29550E+04, .10349E+03,    -.13146E+00, .87787E-03];
end

a = abcd(:,1);
b = abcd(:,2);
c = abcd(:,3);
d = abcd(:,4);










