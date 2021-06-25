addpath /home/sergio/SPECTRA

clc
iX = input('do the FIR bands? (500-605,300-500,140-300,80-150,50-80,30-50,15-30) : ');
if iX > 0
  clc; disp('>>>>>>>> FIR1');   ckd_lookupBIN_FIR1
  clc; disp('>>>>>>>> FIR2');   ckd_lookupBIN_FIR2
  clc; disp('>>>>>>>> FIR3');   ckd_lookupBIN_FIR3
  clc; disp('>>>>>>>> FIR4');   ckd_lookupBIN_FIR4
  clc; disp('>>>>>>>> FIR5');   ckd_lookupBIN_FIR5
  clc; disp('>>>>>>>> FIR6');   ckd_lookupBIN_FIR6
  clc; disp('>>>>>>>> FIR7');   ckd_lookupBIN_FIR7
end

clc
iX = input('do main IR band??? ');
if iX > 0
  ckd_lookupBIN_IR
end

clc
iX = input('do NIR bands??? (2830-3580,3550-5550,5550-8200,8200-12000) : ');
if iX > 0
  clc; disp('>>>>>>>> NIR1');   ckd_lookupBIN_NIR1
  clc; disp('>>>>>>>> NIR2');   ckd_lookupBIN_NIR2
  clc; disp('>>>>>>>> NIR3');   ckd_lookupBIN_NIR3
  clc; disp('>>>>>>>> NIR4');   ckd_lookupBIN_NIR4
end

clc
iX = input('do VIS UV : ');
if iX > 0
  clc; disp('>>>>>>>> VIS1');   ckd_lookupBIN_VIS1
  clc; disp('>>>>>>>> UV 1');   ckd_lookupBIN_UV1
end