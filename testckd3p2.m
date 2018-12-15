f1 = 605; f2 = 705;
topts.CKD = 25; [w,d25] = run8watercontinuum(1,605,f2,'IPFILES/waterone',topts); plot(w,d25)
disp('ret finish CKD2.5'); pause
topts.CKD = 32; [w,d32] = run8watercontinuum(1,605,f2,'IPFILES/waterone',topts); plot(w,d32)
  disp('ret finish CKD3.2'); pause

f1 = 605; f2 = 805;
topts.CKD = 25; [w,d25] = run8watercontinuum(1,605,f2,'IPFILES/waterone',topts); plot(w,d25)
  disp('ret finish CKD2.5'); pause
topts.CKD = 32; [w,d32] = run8watercontinuum(1,605,f2,'IPFILES/waterone',topts); plot(w,d32)
  disp('ret finish CKD3.2'); pause
semilogy(w,d25,'b',w,d32,'r');   disp('ret finish CKD 2.5 3.2'); pause

f1 = 605; f2 = 2805;
topts.CKD = 25; [w,d25] = run8watercontinuum(1,605,f2,'IPFILES/waterone',topts); plot(w,d25)
  disp('ret finish CKD2.5'); pause
topts.CKD = 32; [w,d32] = run8watercontinuum(1,605,f2,'IPFILES/waterone',topts); plot(w,d32)
  disp('ret finish CKD3.2'); pause
semilogy(w,d25,'b.-',w,d32,'r');   disp('ret finish CKD 2.5 3.2'); pause
hl = legend('CKD 2.5','CKD 3.2','location','best'); grid