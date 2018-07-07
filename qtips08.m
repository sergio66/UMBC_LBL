function [a,b,c,d,g] = qtips(gasid,liso);
 
% this was made by find_qnewABCD.m
% returns nuclear degeneracy factors g and Gamache's internal partition
% sum coefficients a, c, c, and d.  
%C  PURPOSE        TOTAL INTERNAL PARTITION SUMS
 
if gasid == 1; g = [1   1   6   6   6  36]; end
if gasid == 2; g = [1   2   1   6   2  12   1   6]; end
if gasid == 3; g = [1  1  1  6  6]; end
if gasid == 4; g = [9   6   6   9  54]; end
if gasid == 5; g = [1   2   1   6   2  12]; end
if gasid == 6; g = [1  2  3]; end
if gasid == 7; g = [1  1  6]; end
if gasid == 8; g = [3  2  3]; end
if gasid == 9; g = [1  1]; end
if gasid == 10; g = [3]; end
if gasid == 11; g = [3  2]; end
if gasid == 12; g = [6]; end
if gasid == 13; g = [2  2  3]; end
if gasid == 14; g = [4]; end
if gasid == 15; g = [8  8]; end
if gasid == 16; g = [8  8]; end
if gasid == 17; g = [12]; end
if gasid == 18; g = [4  4]; end
if gasid == 19; g = [1  1  2  4  1]; end
if gasid == 20; g = [1  2  1]; end
if gasid == 21; g = [8  8]; end
if gasid == 22; g = [1]; end
if gasid == 23; g = [6  12   4]; end
if gasid == 24; g = [4  4]; end
if gasid == 25; g = [1]; end
if gasid == 26; g = [1  8]; end
if gasid == 27; g = [1]; end
if gasid == 28; g = [2]; end
if gasid == 29; g = [1]; end
if gasid == 30; g = [1]; end
if gasid == 31; g = [1  1  4]; end
if gasid == 32; g = [4]; end
 
g = g';
if (length(g) ~= liso)
  fprintf(1,'there is inconsistency in number of isotopes \n'); 
  fprintf(1,'in qtips.m %3i and that in mass.dat %3i \n',length(g),liso);
  error('please check!!!!!!')
  end
 
% fitted/actual at 296K =  
%1      1.0001      1.0004      1.0001      1.0001      1.0001 
if gasid == 1 
    abcd = [ ... 
  -5.02804e+00 , 2.82766e-01, 1.25216e-03 , -5.29467e-07 ; 
  -6.17558e+00 , 2.87130e-01, 1.25090e-03 , -5.12233e-07 ; 
  -3.33768e+01 , 1.65449e+00, 7.76692e-03 , -3.41403e-06 ; 
  -3.36760e+01 , 1.42266e+00, 6.05796e-03 , -2.27960e-06 ; 
  -2.81761e+01 , 1.44304e+00, 6.11330e-03 , -2.28780e-06 ; 
  -1.67456e+02 , 8.59802e+00, 3.65157e-02 , -1.37320e-05 ]; 
  disp('The isotopes are 161  181  171  162  182  172'); 
  end 
 
% fitted/actual at 296K =  
%0.99999     0.99999           1     0.99999     0.99998     0.99999           1     0.99999 
if gasid == 2 
    abcd = [ ... 
  -1.54343e+00 , 9.55877e-01, -7.41606e-04 , 2.71906e-06 ; 
  -2.32001e+00 , 1.89937e+00, -1.45945e-03 , 5.64430e-06 ; 
  -3.21131e+00 , 2.02555e+00, -1.58345e-03 , 5.85593e-06 ; 
  -1.88899e+01 , 1.18200e+01, -9.20705e-03 , 3.39138e-05 ; 
  -4.76161e+00 , 4.02455e+00, -3.11698e-03 , 1.21716e-05 ; 
  -2.87895e+01 , 2.34977e+01, -1.81848e-02 , 7.05459e-05 ; 
  -1.64799e+00 , 1.07439e+00, -8.43865e-04 , 3.15327e-06 ; 
  -1.95541e+01 , 1.25319e+01, -9.82442e-03 , 3.65243e-05 ]; 
  disp('The isotopes are 626  636  628  627  638  637  828  728'); 
  end 
 
% fitted/actual at 296K =  
%1.0002      1.0002      1.0002      1.0002      1.0003 
if gasid == 3 
    abcd = [ ... 
  -2.07031e+02 , 7.49905e+00, 8.53171e-03 , 2.79313e-05 ; 
  -4.46056e+02 , 1.61254e+01, 1.71294e-02 , 6.32212e-05 ; 
  -2.20761e+02 , 7.93956e+00, 8.01524e-03 , 3.14794e-05 ; 
  -2.58192e+03 , 9.34348e+01, 1.02689e-01 , 3.57414e-04 ; 
  -1.28511e+03 , 4.63611e+01, 4.95838e-02 , 1.78462e-04 ]; 
  disp('The isotopes are 666  668  686  667  676'); 
  end 
 
% fitted/actual at 296K =  
%0.99991     0.99992     0.99992     0.99992     0.99994 
if gasid == 4 
    abcd = [ ... 
  4.17376e+00 , 1.55238e+01, -1.12867e-02 , 5.36376e-05 ; 
  7.43139e+00 , 1.02652e+01, -7.23450e-03 , 3.66141e-05 ; 
  3.52587e+00 , 1.07007e+01, -7.79720e-03 , 3.74229e-05 ; 
  4.48101e+00 , 1.64490e+01, -1.21577e-02 , 5.80733e-05 ; 
  2.62781e+01 , 9.59768e+01, -7.03315e-02 , 3.35301e-04 ]; 
  disp('The isotopes are 446  456  546  448  447'); 
  end 
 
% fitted/actual at 296K =  
%0.99996           1     0.99998           1           1           1 
if gasid == 5 
    abcd = [ ... 
  2.97746e-01 , 3.60668e-01, -1.62670e-06 , 7.83652e-09 ; 
  5.88252e-01 , 7.54674e-01, -4.08328e-06 , 1.80282e-08 ; 
  2.88204e-01 , 3.78895e-01, -2.64876e-06 , 9.98693e-09 ; 
  1.74992e+00 , 2.22073e+00, -1.28934e-05 , 5.33223e-08 ; 
  5.69446e-01 , 7.94627e-01, -6.10265e-06 , 2.23577e-08 ; 
  3.60364e+00 , 4.64949e+00, -1.49998e-05 , 9.09864e-08 ]; 
  disp('The isotopes are 26  36  28  27  38  37'); 
  end 
 
% fitted/actual at 296K =  
%1.0001      1.0001      1.0001 
if gasid == 6 
    abcd = [ ... 
  -2.98838e+01 , 1.17448e+00, 2.85006e-03 , 8.87956e-07 ; 
  -5.97145e+01 , 2.34786e+00, 5.70493e-03 , 1.76733e-06 ; 
  -2.43333e+02 , 9.48396e+00, 2.31274e-02 , 7.13789e-06 ]; 
  disp('The isotopes are 211  311  212'); 
  end 
 
% fitted/actual at 296K =  
%1           1     0.99999 
if gasid == 7 
    abcd = [ ... 
  3.99968e-01 , 7.33471e-01, -4.97189e-05 , 1.01140e-07 ; 
  -3.96665e+00 , 1.55603e+00, -1.22324e-04 , 2.47066e-07 ; 
  -2.29272e+01 , 9.07765e+00, -6.67353e-04 , 1.34947e-06 ]; 
  disp('The isotopes are 66  68  67'); 
  end 
 
% fitted/actual at 296K =  
%1.0003      1.0002      1.0002 
if gasid == 8 
    abcd = [ ... 
  -3.56407e+01 , 2.79503e+00, 5.61754e-03 , -5.45441e-06 ; 
  -2.47529e+01 , 1.93302e+00, 3.87320e-03 , -3.75261e-06 ; 
  -3.78764e+01 , 2.95131e+00, 5.90222e-03 , -5.71108e-06 ]; 
  disp('The isotopes are 46  56  48'); 
  end 
 
% fitted/actual at 296K =  
%1.0001      1.0001 
if gasid == 9 
    abcd = [ ... 
  -2.80225e+02 , 1.16078e+01, 2.09675e-02 , 5.19798e-05 ; 
  -2.81736e+02 , 1.16646e+01, 2.10502e-02 , 5.22151e-05 ]; 
  disp('The isotopes are 626  646'); 
  end 
 
% fitted/actual at 296K =  
%1.0001 
if gasid == 10 
    abcd = [ ... 
  -6.62785e+02 , 2.62072e+01, 5.93160e-02 , 4.96675e-05 ]; 
  disp('The isotopes are 646'); 
  end 
 
% fitted/actual at 296K =  
%1.0001      1.0001 
if gasid == 11 
    abcd = [ ... 
  -7.44915e+01 , 3.20634e+00, 9.30738e-03 , 1.36491e-06 ; 
  -4.98390e+01 , 2.14309e+00, 6.21456e-03 , 9.17732e-07 ]; 
  disp('The isotopes are 4111  5111'); 
  end 
 
% fitted/actual at 296K =  
%0.99978 
if gasid == 12 
    abcd = [ ... 
  -2.97279e+04 , 8.07443e+02, -2.76346e+00 , 9.52103e-03 ]; 
  disp('The isotopes are 146'); 
  end 
 
% fitted/actual at 296K =  
%1.0003      1.0003      1.0003 
if gasid == 13 
    abcd = [ ... 
  8.21220e+00 , 1.64964e-01, 3.80932e-04 , -3.86712e-07 ; 
  8.15332e+00 , 1.66832e-01, 3.81099e-04 , -3.86267e-07 ; 
  8.31643e+00 , 4.65260e-01, 1.02108e-03 , -1.00740e-06 ]; 
  disp('The isotopes are 61  81  62'); 
  end 
 
% fitted/actual at 296K =  
%1 
if gasid == 14 
    abcd = [ ... 
  1.54307e+00 , 1.33499e-01, 6.48797e-06 , -6.17308e-09 ]; 
  disp('The isotopes are 19'); 
  end 
 
% fitted/actual at 296K =  
%0.99998           1 
if gasid == 15 
    abcd = [ ... 
  2.86543e+00 , 5.31015e-01, 8.50433e-06 , -5.20203e-09 ; 
  2.85807e+00 , 5.31948e-01, 7.86017e-06 , -4.25403e-09 ]; 
  disp('The isotopes are 15  17'); 
  end 
 
% fitted/actual at 296K =  
%0.99999     0.99999 
if gasid == 16 
    abcd = [ ... 
  2.80758e+00 , 6.64916e-01, 6.56549e-06 , -7.73871e-10 ; 
  2.80971e+00 , 6.65069e-01, 6.83953e-06 , -1.19615e-09 ]; 
  disp('The isotopes are 19  11'); 
  end 
 
% fitted/actual at 296K =  
%0.99999 
if gasid == 17 
    abcd = [ ... 
  4.07211e+00 , 1.29849e+00, 1.77874e-06 , 1.60299e-08 ]; 
  disp('The isotopes are 17'); 
  end 
 
% fitted/actual at 296K =  
%0.99985     0.99984 
if gasid == 18 
    abcd = [ ... 
  9.34281e+01 , 6.82248e+00, 1.26240e-02 , 2.12756e-06 ; 
  9.50449e+01 , 6.94102e+00, 1.28163e-02 , 2.28590e-06 ]; 
  disp('The isotopes are 56  76'); 
  end 
 
% fitted/actual at 296K =  
%0.99989     0.99985     0.99985     0.99989     0.99991 
if gasid == 19 
    abcd = [ ... 
  1.66944e+00 , 3.58126e+00, -3.39442e-03 , 1.76043e-05 ; 
  1.56044e+00 , 3.67550e+00, -3.52908e-03 , 1.82388e-05 ; 
  6.98620e+00 , 7.11367e+00, -6.56329e-03 , 3.64857e-05 ; 
  6.53873e+00 , 1.45146e+01, -1.38383e-02 , 7.16866e-05 ; 
  2.62005e+00 , 3.80276e+00, -3.61056e-03 , 1.93447e-05 ]; 
  disp('The isotopes are 622  624  632  623  822'); 
  end 
 
% fitted/actual at 296K =  
%1.0001      1.0001      1.0001 
if gasid == 20 
    abcd = [ ... 
  -1.37946e+02 , 5.44463e+00, 1.49652e-02 , 2.39793e-06 ; 
  -2.83349e+02 , 1.11713e+01, 3.06619e-02 , 4.95420e-06 ; 
  -1.44926e+02 , 5.71380e+00, 1.56906e-02 , 2.53135e-06 ]; 
  disp('The isotopes are 126  136  128'); 
  end 
 
% fitted/actual at 296K =  
%1.0002      1.0002 
if gasid == 21 
    abcd = [ ... 
  -9.09778e+02 , 3.66590e+01, 8.47874e-02 , 7.35832e-05 ; 
  -9.24461e+02 , 3.72813e+01, 8.64210e-02 , 7.46990e-05 ]; 
  disp('The isotopes are 165  167'); 
  end 
 
% fitted/actual at 296K =  
%1 
if gasid == 22 
    abcd = [ ... 
  1.41961e+00 , 1.56718e+00, 2.37010e-06 , 1.79074e-08 ]; 
  disp('The isotopes are 44'); 
  end 
 
% fitted/actual at 296K =  
%1      1.0001           1 
if gasid == 23 
    abcd = [ ... 
  -3.12153e+00 , 3.00982e+00, -2.05181e-03 , 7.22265e-06 ; 
  -5.89465e+00 , 6.16964e+00, -4.18727e-03 , 1.49208e-05 ; 
  -1.87538e+00 , 2.07846e+00, -1.41071e-03 , 5.07762e-06 ]; 
  disp('The isotopes are 124  134  125'); 
  end 
 
% fitted/actual at 296K =  
%1.0002      1.0002 
if gasid == 24 
    abcd = [ ... 
  -8.76466e+03 , 2.87589e+02, 5.57016e-02 , 1.33462e-03 ; 
  -8.89699e+03 , 2.92018e+02, 5.72166e-02 , 1.35486e-03 ]; 
  disp('The isotopes are 215  217'); 
  end 
 
% fitted/actual at 296K =  
%0.99997 
if gasid == 25 
    abcd = [ ... 
  -3.51392e+02 , 1.30860e+01, 3.50311e-02 , 1.24476e-04 ]; 
  disp('The isotopes are 1661'); 
  end 
 
% fitted/actual at 296K =  
%0.99999     0.99997 
if gasid == 26 
    abcd = [ ... 
  -8.39273e+00 , 1.44579e+00, -2.54589e-03 , 8.38746e-06 ; 
  -3.36284e+01 , 5.78435e+00, -1.01886e-02 , 3.35575e-05 ]; 
  disp('The isotopes are 1221  1231'); 
  end 
 
% fitted/actual at 296K =  
%0.99953 
if gasid == 27 
    abcd = [ ... 
  -8.47500e+03 , 2.26466e+02, -6.35126e-01 , 2.61955e-03 ]; 
  disp('The isotopes are 1221'); 
  end 
 
% fitted/actual at 296K =  
%1.0002 
if gasid == 28 
    abcd = [ ... 
  -1.82689e+02 , 6.84261e+00, 1.20018e-02 , 1.36881e-05 ]; 
  disp('The isotopes are 1111'); 
  end 
 
% fitted/actual at 296K =  
%1.0001 
if gasid == 29 
    abcd = [ ... 
  -6.64111e+03 , 2.02915e+02, -3.51663e-01 , 1.82936e-03 ]; 
  disp('The isotopes are 269'); 
  end 
 
% fitted/actual at 296K =  
%0.98207 
if gasid == 30 
    abcd = [ ... 
  -2.08724e+06 , 4.26101e+04, -2.54542e+02 , 5.15565e-01 ]; 
  disp('The isotopes are 29'); 
  end 
 
% fitted/actual at 296K =  
%1.0001      1.0001      1.0001 
if gasid == 31 
    abcd = [ ... 
  -1.81295e+01 , 8.61898e-01, 3.31448e-03 , -9.35442e-07 ; 
  -1.81846e+01 , 8.64200e-01, 3.32234e-03 , -9.36832e-07 ; 
  -7.25663e+01 , 3.45129e+00, 1.32793e-02 , -3.75208e-06 ]; 
  disp('The isotopes are 121  141  131'); 
  end 
 
% fitted/actual at 296K =  
%1.0001 
if gasid == 32 
    abcd = [ ... 
  -3.52215e+03 , 1.09068e+02, -1.29304e-01 , 8.36893e-04 ]; 
  disp('The isotopes are 126'); 
  end 
 
a = abcd(:,1); 
b = abcd(:,2); 
c = abcd(:,3); 
d = abcd(:,4); 
