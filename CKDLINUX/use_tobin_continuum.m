%%% this is so that we can use a MexSized version of calconCKD5.f instead
%%% of relying on ckd_lookupBIN_v5_ieee_le.m

%%% what a pain in the butt!

ckd1 = '/home/sergio/KCARTADATA/General/CKDieee_le/CKDSelf1.bin'; 
[kx, fr, temp] = contread(ckd1); 

cd /carrot/s1/strow/Tobin/Tobin_radish/NISTdata2/New/ 
WN = 1300:1:1800;  
[Cs0_296,Cs0_260,Cf0_296]=makecons(WN,7); 
T0 = 296; 
raTFAC = (temp - T0)/(260.-T0); 
raSH2OT0 = Cs0_296; 
raSH2OT1 = Cs0_260; 
for iL = 1 : length(temp) 
  raSH2O(:,iL) = raSH2OT0.*(raSH2OT1./raSH2OT0).^raTFAC(iL); 
  raFH2O(:,iL) = Cf0_296; 
  end 
raSH2O = raSH2O'; 
raFH2O = raFH2O'; 
 
