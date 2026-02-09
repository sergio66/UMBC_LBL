function [a,b,c,d,g] = qtips16(gasid,liso);
  
% this was made by find_qnewABCD_H16.m
% returns nuclear degeneracy factors g and Gamache's internal partition
% sum coefficients a, b, c, and d.  
%C  PURPOSE        TOTAL INTERNAL PARTITION SUMS

if gasid == 1; g = [1   1   6   6   6  36   1]; end
if gasid == 2; g = [1   2   1   6   2  12   1   6   1   2  12]; end
if gasid == 3; g = [1  1  1  6  6]; end
if gasid == 4; g = [9   6   6   9  54]; end
if gasid == 5; g = [1   2   1   6   2  12]; end
if gasid == 6; g = [1  2  3  6]; end
if gasid == 7; g = [1  1  6]; end
if gasid == 8; g = [3  2  3]; end
if gasid == 9; g = [1  1]; end
if gasid == 10; g = [3]; end
if gasid == 11; g = [3  2]; end
if gasid == 12; g = [6  4]; end
if gasid == 13; g = [2  2  3]; end
if gasid == 14; g = [4  6]; end
if gasid == 15; g = [8   8  12  12]; end
if gasid == 16; g = [8   8  12  12]; end
if gasid == 17; g = [12  18]; end
if gasid == 18; g = [4  4]; end
if gasid == 19; g = [1  1  2  4  1]; end
if gasid == 20; g = [1  2  1]; end
if gasid == 21; g = [8  8]; end
if gasid == 22; g = [1  6]; end
if gasid == 23; g = [6  12   4]; end
if gasid == 24; g = [4  4]; end
if gasid == 25; g = [1]; end
if gasid == 26; g = [1  8  6]; end
if gasid == 27; g = [1  2]; end
if gasid == 28; g = [2]; end
if gasid == 29; g = [1  2]; end
if gasid == 30; g = [1]; end
if gasid == 31; g = [1  1  4]; end
if gasid == 32; g = [4]; end
if gasid == 33; g = [2]; end
if gasid == 34; g = [1]; end
if gasid == 35; g = [12  12]; end
if gasid == 36; g = [3]; end
if gasid == 37; g = [8  8]; end
if gasid == 38; g = [1  2]; end
if gasid == 39; g = [2]; end
if gasid == 40; g = [4  4]; end
if gasid == 41; g = [3]; end
if gasid == 42; g = [1]; end
if gasid == 43; g = [1]; end
if gasid == 44; g = [6]; end
if gasid == 45; g = [1  6]; end
if gasid == 46; g = [1  1  2  4]; end
if gasid == 47; g = [1]; end
if gasid == 48; g = [1]; end
if gasid == 49; g = [1  16]; end
 
g = g';
if (length(g) ~= liso)
  fprintf(1,'there is inconsistency in number of isotopes \n'); 
  fprintf(1,'in qtips.m %3i and that in mass.dat %3i \n',length(g),liso);
  error('please check!!!!!!')
end
 
% fitted/actual at 296K =  
%1.0001      1.0001      1.0001      1.0001      1.0001      1.0001      1.0001 
if gasid == 1 
    abcd = [ ... 
  -5.03394e+00 , 2.82701e-01, 1.25234e-03 , -5.31230e-07 ; 
  -5.08465e+00 , 2.85152e-01, 1.26238e-03 , -5.34443e-07 ; 
  -3.03343e+01 , 1.70338e+00, 7.54910e-03 , -3.20274e-06 ; 
  -2.76699e+01 , 1.42157e+00, 6.06419e-03 , -2.29901e-06 ; 
  -2.80730e+01 , 1.43999e+00, 6.13628e-03 , -2.31943e-06 ; 
  -1.66834e+02 , 8.58464e+00, 3.65591e-02 , -1.35022e-05 ; 
  -3.65796e+01 , 1.74450e+00, 6.79623e-03 , -1.82559e-06 ]; 
  disp('The isotopes are 161  181  171  162  182  172  262'); 
end 
 
% fitted/actual at 296K =  
%1           1           1           1     0.99996           1           1           1     0.99997     0.99999     0.99999 
if gasid == 2 
    abcd = [ ... 
  -1.65407e+00 , 9.58162e-01, -7.65288e-04 , 2.74493e-06 ; 
  -2.44629e+00 , 1.90180e+00, -1.49188e-03 , 5.66308e-06 ; 
  -3.55661e+00 , 2.03289e+00, -1.65085e-03 , 5.94879e-06 ; 
  -2.06391e+01 , 1.18565e+01, -9.55268e-03 , 3.43467e-05 ; 
  -5.24975e+00 , 4.03437e+00, -3.21680e-03 , 1.22759e-05 ; 
  -3.04229e+01 , 2.35296e+01, -1.86113e-02 , 7.08549e-05 ; 
  -1.91991e+00 , 1.08022e+00, -8.92777e-04 , 3.23212e-06 ; 
  -2.21912e+01 , 1.25883e+01, -1.03145e-02 , 3.72636e-05 ; 
  -6.43127e+01 , 3.66937e+01, -2.98289e-02 , 1.07509e-04 ; 
  -2.82091e+00 , 2.14310e+00, -1.73806e-03 , 6.66988e-06 ; 
  -3.25966e+01 , 2.49767e+01, -2.00829e-02 , 7.68876e-05 ]; 
  disp('The isotopes are 626  636  628  627  638  637  828  827  727  838  837'); 
end 
 
% fitted/actual at 296K =  
%1.0003      1.0002      1.0002      1.0002      1.0002 
if gasid == 3 
    abcd = [ ... 
  -2.07620e+02 , 7.50781e+00, 8.49485e-03 , 2.79786e-05 ; 
  -4.46416e+02 , 1.61282e+01, 1.71304e-02 , 6.32027e-05 ; 
  -2.21301e+02 , 7.94769e+00, 7.98017e-03 , 3.15258e-05 ; 
  -2.58709e+03 , 9.35050e+01, 1.02420e-01 , 3.57735e-04 ; 
  -1.28793e+03 , 4.64025e+01, 4.94090e-02 , 1.78687e-04 ]; 
  disp('The isotopes are 666  668  686  667  676'); 
end 
 
% fitted/actual at 296K =  
%0.99993     0.99991     0.99992     0.99993     0.99991 
if gasid == 4 
    abcd = [ ... 
  4.45285e+00 , 1.55167e+01, -1.13743e-02 , 5.33541e-05 ; 
  7.66762e+00 , 1.02619e+01, -7.21842e-03 , 3.65911e-05 ; 
  3.60297e+00 , 1.07000e+01, -7.79599e-03 , 3.74234e-05 ; 
  4.71110e+00 , 1.64458e+01, -1.21440e-02 , 5.80562e-05 ; 
  2.78393e+01 , 9.59587e+01, -7.02611e-02 , 3.35219e-04 ]; 
  disp('The isotopes are 446  456  546  448  447'); 
end 
 
% fitted/actual at 296K =  
%1           1           1           1           1     0.99996 
if gasid == 5 
    abcd = [ ... 
  2.98465e-01 , 3.62245e-01, -3.46084e-06 , 7.70333e-09 ; 
  5.78642e-01 , 7.58118e-01, -8.99633e-06 , 1.92659e-08 ; 
  2.89806e-01 , 3.80502e-01, -4.52749e-06 , 9.73103e-09 ; 
  1.76607e+00 , 2.23018e+00, -2.36837e-05 , 5.18435e-08 ; 
  5.59504e-01 , 7.98152e-01, -1.12350e-05 , 2.36488e-08 ; 
  3.41460e+00 , 4.67304e+00, -6.06093e-05 , 1.28328e-07 ]; 
  disp('The isotopes are 26  36  28  27  38  37'); 
end 
 
% fitted/actual at 296K =  
%1.0001      1.0001      1.0001      1.0001 
if gasid == 6 
    abcd = [ ... 
  -3.01118e+01 , 1.17864e+00, 2.82546e-03 , 9.33510e-07 ; 
  -5.98998e+01 , 2.35060e+00, 5.69288e-03 , 1.78395e-06 ; 
  -2.73066e+02 , 1.01173e+01, 1.89515e-02 , 1.59335e-05 ; 
  -5.47726e+02 , 2.02772e+01, 3.77497e-02 , 3.23348e-05 ]; 
  disp('The isotopes are 211  311  212  312'); 
end 
 
% fitted/actual at 296K =  
%1           1           1 
if gasid == 7 
    abcd = [ ... 
  4.23903e-01 , 7.32863e-01, -4.69083e-05 , 9.61877e-08 ; 
  -9.68533e-01 , 1.55585e+00, -1.23252e-04 , 2.49462e-07 ; 
  -5.08467e+00 , 9.07803e+00, -6.78244e-04 , 1.37046e-06 ]; 
  disp('The isotopes are 66  68  67'); 
end 
 
% fitted/actual at 296K =  
%1.0003      1.0003      1.0003 
if gasid == 8 
    abcd = [ ... 
  -3.24106e+01 , 2.75380e+00, 5.78045e-03 , -5.65707e-06 ; 
  -2.27494e+01 , 1.90731e+00, 3.97502e-03 , -3.87913e-06 ; 
  -3.50272e+01 , 2.91437e+00, 6.04991e-03 , -5.89599e-06 ]; 
  disp('The isotopes are 46  56  48'); 
end 
 
% fitted/actual at 296K =  
%1.0001      1.0001 
if gasid == 9 
    abcd = [ ... 
  -2.80642e+02 , 1.16130e+01, 2.09526e-02 , 5.19878e-05 ; 
  -2.82243e+02 , 1.16713e+01, 2.10264e-02 , 5.22413e-05 ]; 
  disp('The isotopes are 626  646'); 
end 
 
% fitted/actual at 296K =  
%1.0002 
if gasid == 10 
    abcd = [ ... 
  -6.64424e+02 , 2.62297e+01, 5.92275e-02 , 4.97777e-05 ]; 
  disp('The isotopes are 646'); 
end 
 
% fitted/actual at 296K =  
%1.0002      1.0002 
if gasid == 11 
    abcd = [ ... 
  -7.62704e+01 , 3.24093e+00, 9.10102e-03 , 1.73588e-06 ; 
  -5.08208e+01 , 2.16436e+00, 6.08556e-03 , 1.17421e-06 ]; 
  disp('The isotopes are 4111  5111'); 
end 
 
% fitted/actual at 296K =  
%0.99978     0.99977 
if gasid == 12 
    abcd = [ ... 
  -2.96799e+04 , 8.06247e+02, -2.75407e+00 , 9.49373e-03 ; 
  -2.00691e+04 , 5.43813e+02, -1.88126e+00 , 6.44036e-03 ]; 
  disp('The isotopes are 146  156'); 
end 
 
% fitted/actual at 296K =  
%1.0003      1.0003      1.0002 
if gasid == 13 
    abcd = [ ... 
  8.15605e+00 , 1.65595e-01, 3.78222e-04 , -3.83160e-07 ; 
  8.11347e+00 , 1.67379e-01, 3.78849e-04 , -3.83419e-07 ; 
  8.25537e+00 , 4.66100e-01, 1.01768e-03 , -1.00315e-06 ]; 
  disp('The isotopes are 61  81  62'); 
end 
 
% fitted/actual at 296K =  
%1           1 
if gasid == 14 
    abcd = [ ... 
  1.54087e+00 , 1.33532e-01, 6.37621e-06 , -5.97366e-09 ; 
  2.15409e+00 , 3.82778e-01, 6.37029e-06 , -3.82965e-09 ]; 
  disp('The isotopes are 19  29'); 
end 
 
% fitted/actual at 296K =  
%1           1           1     0.99999 
if gasid == 15 
    abcd = [ ... 
  2.86398e+00 , 5.31031e-01, 8.40895e-06 , -5.01147e-09 ; 
  2.86178e+00 , 5.31878e-01, 8.11094e-06 , -4.50498e-09 ; 
  3.94533e+00 , 1.54952e+00, -8.91778e-06 , 3.69997e-08 ; 
  3.95147e+00 , 1.55400e+00, -8.78794e-06 , 3.71053e-08 ]; 
  disp('The isotopes are 15  17  25  27'); 
end 
 
% fitted/actual at 296K =  
%1           1           1           1 
if gasid == 16 
    abcd = [ ... 
  2.80517e+00 , 6.64938e-01, 6.45015e-06 , -7.08339e-10 ; 
  2.80403e+00 , 6.65159e-01, 6.37245e-06 , -5.68755e-10 ; 
  3.57469e+00 , 1.97315e+00, -4.78626e-05 , 1.14454e-07 ; 
  3.57139e+00 , 1.97440e+00, -4.81232e-05 , 1.14870e-07 ]; 
  disp('The isotopes are 19  11  29  21'); 
end 
 
% fitted/actual at 296K =  
%1     0.99997 
if gasid == 17 
    abcd = [ ... 
  4.07210e+00 , 1.29851e+00, 1.47828e-06 , 1.65719e-08 ; 
  4.06613e+00 , 3.88614e+00, -2.26261e-04 , 4.82979e-07 ]; 
  disp('The isotopes are 17  27'); 
end 
 
% fitted/actual at 296K =  
%0.99986     0.99985 
if gasid == 18 
    abcd = [ ... 
  9.40604e+01 , 6.81335e+00, 1.26640e-02 , 2.07360e-06 ; 
  9.57055e+01 , 6.93176e+00, 1.28554e-02 , 2.23517e-06 ]; 
  disp('The isotopes are 56  76'); 
end 
 
% fitted/actual at 296K =  
%0.9999     0.99987     0.99989     0.99989     0.99986 
if gasid == 19 
    abcd = [ ... 
  1.75081e+00 , 3.58034e+00, -3.39130e-03 , 1.76013e-05 ; 
  1.70990e+00 , 3.67334e+00, -3.51931e-03 , 1.82254e-05 ; 
  7.25619e+00 , 7.11005e+00, -6.54867e-03 , 3.64680e-05 ; 
  6.91429e+00 , 1.45105e+01, -1.38251e-02 , 7.16743e-05 ; 
  2.84651e+00 , 3.79919e+00, -3.59326e-03 , 1.93196e-05 ]; 
  disp('The isotopes are 622  624  632  623  822'); 
end 
 
% fitted/actual at 296K =  
%1.0001      1.0001      1.0001 
if gasid == 20 
    abcd = [ ... 
  -1.38788e+02 , 5.45906e+00, 1.48673e-02 , 2.51190e-06 ; 
  -2.84039e+02 , 1.11807e+01, 3.06242e-02 , 5.00161e-06 ; 
  -1.45325e+02 , 5.71928e+00, 1.56692e-02 , 2.55670e-06 ]; 
  disp('The isotopes are 126  136  128'); 
end 
 
% fitted/actual at 296K =  
%1.0002      1.0002 
if gasid == 21 
    abcd = [ ... 
  -9.11641e+02 , 3.66824e+01, 8.47052e-02 , 7.36702e-05 ; 
  -9.27873e+02 , 3.73310e+01, 8.62174e-02 , 7.49540e-05 ]; 
  disp('The isotopes are 165  167'); 
end 
 
% fitted/actual at 296K =  
%1     0.99999 
if gasid == 22 
    abcd = [ ... 
  1.42985e+00 , 1.57362e+00, -6.51080e-06 , 1.73729e-08 ; 
  1.87759e+00 , 2.17056e+00, -1.13269e-05 , 2.79607e-08 ]; 
  disp('The isotopes are 44  45'); 
end 
 
% fitted/actual at 296K =  
%1.0001           1      1.0001 
if gasid == 23 
    abcd = [ ... 
  -2.94210e+00 , 3.00541e+00, -2.05134e-03 , 7.14596e-06 ; 
  -6.25131e+00 , 6.16769e+00, -4.23964e-03 , 1.47740e-05 ; 
  -3.85596e+00 , 2.10112e+00, -1.62975e-03 , 5.39968e-06 ]; 
  disp('The isotopes are 124  134  125'); 
end 
 
% fitted/actual at 296K =  
%1.0002      1.0002 
if gasid == 24 
    abcd = [ ... 
  -4.38567e+03 , 1.43822e+02, 2.78256e-02 , 6.67246e-04 ; 
  -4.45515e+03 , 1.46093e+02, 2.83215e-02 , 6.77733e-04 ]; 
  disp('The isotopes are 215  217'); 
end 
 
% fitted/actual at 296K =  
%1 
if gasid == 25 
    abcd = [ ... 
  -3.04164e+02 , 1.21435e+01, 4.06816e-02 , 1.15423e-04 ]; 
  disp('The isotopes are 1661'); 
end 
 
% fitted/actual at 296K =  
%0.99998     0.99997     0.99885 
if gasid == 26 
    abcd = [ ... 
  -8.03675e+00 , 1.43785e+00, -2.50354e-03 , 8.26034e-06 ; 
  -3.35324e+01 , 5.78317e+00, -1.01839e-02 , 3.35514e-05 ; 
  -5.48522e+01 , 5.54015e+00, -1.22696e-02 , 4.12566e-05 ]; 
  disp('The isotopes are 1221  1231  1222'); 
end 
 
% fitted/actual at 296K =  
%0.99953     0.99953 
if gasid == 27 
    abcd = [ ... 
  -8.47511e+03 , 2.26495e+02, -6.35338e-01 , 2.62000e-03 ; 
  -4.32736e+03 , 1.15643e+02, -3.24477e-01 , 1.33806e-03 ]; 
  disp('The isotopes are 1221  1231'); 
end 
 
% fitted/actual at 296K =  
%1.0002 
if gasid == 28 
    abcd = [ ... 
  -1.83721e+02 , 6.86141e+00, 1.18939e-02 , 1.39107e-05 ]; 
  disp('The isotopes are 1111'); 
end 
 
% fitted/actual at 296K =  
%1.0002      1.0002 
if gasid == 29 
    abcd = [ ... 
  -6.63953e+03 , 2.02840e+02, -3.51121e-01 , 1.82774e-03 ; 
  -1.32805e+04 , 4.05705e+02, -7.02387e-01 , 3.65596e-03 ]; 
  disp('The isotopes are 269  369'); 
end 
 
% fitted/actual at 296K =  
%0.98238 
if gasid == 30 
    abcd = [ ... 
  -2.08780e+06 , 4.26226e+04, -2.54620e+02 , 5.15727e-01 ]; 
  disp('The isotopes are 29'); 
end 
 
% fitted/actual at 296K =  
%1.0001      1.0001      1.0001 
if gasid == 31 
    abcd = [ ... 
  -1.81080e+01 , 8.60483e-01, 3.33982e-03 , -9.01023e-07 ; 
  -1.82508e+01 , 8.65146e-01, 3.31840e-03 , -9.31770e-07 ; 
  -7.29360e+01 , 3.45688e+00, 1.32542e-02 , -3.71772e-06 ]; 
  disp('The isotopes are 121  141  131'); 
end 
 
% fitted/actual at 296K =  
%1.0001 
if gasid == 32 
    abcd = [ ... 
  -3.52672e+03 , 1.09134e+02, -1.29581e-01 , 8.37262e-04 ]; 
  disp('The isotopes are 126'); 
end 
 
% fitted/actual at 296K =  
%1.0001 
if gasid == 33 
    abcd = [ ... 
  -1.85848e+02 , 7.78511e+00, 2.54339e-02 , -1.77304e-06 ]; 
  disp('The isotopes are 166'); 
end 
 
% fitted/actual at 296K =  
%0.14878 
if gasid == 34 
    abcd = [ ... 
  1.00000e+00 , 2.10795e-17, -9.85073e-20 , 1.38945e-22 ]; 
  disp('The isotopes are 6'); 
end 
 
% fitted/actual at 296K =  
%0.99761     0.99762 
if gasid == 35 
    abcd = [ ... 
  -1.99763e+06 , 4.25674e+04, -2.57393e+02 , 6.44953e-01 ; 
  -2.04856e+06 , 4.36522e+04, -2.63951e+02 , 6.61375e-01 ]; 
  disp('The isotopes are 5646  7646'); 
end 
 
% fitted/actual at 296K =  
%0.99999 
if gasid == 36 
    abcd = [ ... 
  9.55611e-01 , 1.04998e+00, -3.95264e-06 , 1.08419e-08 ]; 
  disp('The isotopes are 46'); 
end 
 
% fitted/actual at 296K =  
%1.0002      1.0002 
if gasid == 37 
    abcd = [ ... 
  -1.28496e+03 , 5.27852e+01, 1.12495e-01 , 1.59960e-04 ; 
  -1.27819e+03 , 5.25529e+01, 1.12118e-01 , 1.59711e-04 ]; 
  disp('The isotopes are 169  161'); 
end 
 
% fitted/actual at 296K =  
%1.0002      1.0002 
if gasid == 38 
    abcd = [ ... 
  -1.06270e+03 , 3.20417e+01, -2.46305e-02 , 1.84327e-04 ; 
  -4.35012e+03 , 1.31151e+02, -1.00782e-01 , 7.54471e-04 ]; 
  disp('The isotopes are 221  231'); 
end 
 
% fitted/actual at 296K =  
%0.99977 
if gasid == 39 
    abcd = [ ... 
  -4.23201e+03 , 1.34986e+02, -6.20886e-02 , 1.55277e-03 ]; 
  disp('The isotopes are 2161'); 
end 
 
% fitted/actual at 296K =  
%1.0002      1.0002 
if gasid == 40 
    abcd = [ ... 
  -6.41565e+03 , 2.09236e+02, -5.07900e-02 , 1.23380e-03 ; 
  -6.43792e+03 , 2.10014e+02, -5.12088e-02 , 1.24043e-03 ]; 
  disp('The isotopes are 219  211'); 
end 
 
% fitted/actual at 296K =  
%0.9996 
if gasid == 41 
    abcd = [ ... 
  -9.24028e+03 , 2.58781e+02, -8.28660e-01 , 3.61997e-03 ]; 
  disp('The isotopes are 2124'); 
end 
 
% fitted/actual at 296K =  
%8.2461e-06 
if gasid == 42 
    abcd = [ ... 
  1.00000e+00 , 2.10795e-17, -9.85073e-20 , 1.38945e-22 ]; 
  disp('The isotopes are 29'); 
end 
 
% fitted/actual at 296K =  
%0.99697 
if gasid == 43 
    abcd = [ ... 
  -4.13718e+03 , 9.48475e+01, -5.55119e-01 , 1.32985e-03 ]; 
  disp('The isotopes are 2211'); 
end 
 
% fitted/actual at 296K =  
%0.99891 
if gasid == 44 
    abcd = [ ... 
  -5.64627e+03 , 1.45091e+02, -7.94112e-01 , 2.19926e-03 ]; 
  disp('The isotopes are 1224'); 
end 
 
% fitted/actual at 296K =  
%1.0001      1.0001 
if gasid == 45 
    abcd = [ ... 
  -5.47970e-01 , 3.57918e-02, -3.97537e-05 , 4.27607e-08 ; 
  2.78373e+00 , 8.66668e-02, 2.37226e-05 , -2.46052e-08 ]; 
  disp('The isotopes are 11  12'); 
end 
 
% fitted/actual at 296K =  
%1           1           1           1 
if gasid == 46 
    abcd = [ ... 
  -6.66930e-01 , 8.71899e-01, -1.26628e-04 , 2.81758e-07 ; 
  -7.11155e-01 , 8.86674e-01, -1.32669e-04 , 2.95116e-07 ; 
  -1.61113e+00 , 1.85089e+00, -2.94878e-04 , 6.59526e-07 ; 
  -2.75944e+00 , 3.51787e+00, -5.18995e-04 , 1.15469e-06 ]; 
  disp('The isotopes are 22  24  32  23'); 
end 
 
% fitted/actual at 296K =  
%0.9999 
if gasid == 47 
    abcd = [ ... 
  -7.60525e+02 , 2.28363e+01, -6.01361e-02 , 2.71932e-04 ]; 
  disp('The isotopes are 26'); 
end 
 
% fitted/actual at 296K =  
%0.99973 
if gasid == 48 
    abcd = [ ... 
  -1.54239e+03 , 5.24879e+01, -2.40376e-01 , 8.73150e-04 ]; 
  disp('The isotopes are 4224'); 
end 
 
% fitted/actual at 296K =  
%0.99968     0.99965 
if gasid == 49 
    abcd = [ ... 
  -2.34296e+05 , 5.94007e+03, -2.78350e+01 , 9.23237e-02 ; 
  -4.82169e+05 , 1.22222e+04, -5.72861e+01 , 1.89942e-01 ]; 
  disp('The isotopes are 2655  2657'); 
end 
 
a = abcd(:,1); 
b = abcd(:,2); 
c = abcd(:,3); 
d = abcd(:,4); 
