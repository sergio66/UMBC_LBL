function [matr]=getmatr(band,prb,freqq)
%this function get the appropriate matrix for us to use kfull/klor

%done from findratio.m 
if ((band == 720) & ( (prb == 'P') | (prb == 'p')))  
   matr=[... 
   6.0000e+02   4.3653e-01   4.3073e-01   4.2771e-01   4.2594e-01   4.2488e-01   4.2419e-01   4.2376e-01; 
   6.4000e+02   4.3700e-01   4.3162e-01   4.2912e-01   4.2796e-01   4.2763e-01   4.2779e-01   4.2840e-01; 
   7.5000e+02   4.4209e-01   4.3752e-01   4.3544e-01   4.3440e-01   4.3395e-01   4.3384e-01   4.3396e-01; 
   7.8000e+02   4.3863e-01   4.3321e-01   4.3046e-01   4.2885e-01   4.2792e-01   4.2736e-01   4.2709e-01]; 
elseif ((band == 720) & ( (prb == 'R') | (prb == 'r')))  
   matr=[... 
   6.6000e+02   4.3158e-01   4.2923e-01   4.2830e-01   4.2760e-01   4.2700e-01   4.2647e-01   4.2593e-01; 
   6.8000e+02   4.3224e-01   4.2996e-01   4.2923e-01   4.2876e-01   4.2839e-01   4.2809e-01   4.2776e-01; 
   8.0000e+02   4.3365e-01   4.3268e-01   4.3237e-01   4.3212e-01   4.3210e-01   4.3239e-01   4.3294e-01; 
   8.3000e+02   4.3267e-01   4.3115e-01   4.3036e-01   4.2958e-01   4.2891e-01   4.2837e-01   4.2786e-01]; 
else
  matr=ones(4,8);
  mn=floor(min(freqq)); mx=ceil(max(freqq));
  %%%doVmixSimple uses abs(f-15) to see if it wants to do k=klor*0.5
  %%%so i've used abs-16 to be safe
  matr(:,1)=[mn-20.0 mn-16.0 mx+16.0 mx+20.0]';
  end
