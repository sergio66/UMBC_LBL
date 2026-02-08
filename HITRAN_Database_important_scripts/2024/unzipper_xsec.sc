# assumes this is a fresh copy from H20XY, so all g51 - 81 are unzipped by default ...
# if you don;t want son gases unzipped, just comment them out eg here we comment out g51 amd g81
#str='Unzipping g81"
#echo $str
#cd ../SF6_g81
#  unzip *.zip
#
#str='Unzipping g51"
#echo $str
#cd CFC-11_CCl3F_g51
#cp ../CFC-11.zip .
#  unzip *.zip

########################################################################

ls -ltR XSEC/
cd XSEC

########################################################################

str="Unzipping g51"
echo $str
cd CFC-11_CCl3F_g51
#cp ../CFC-11.zip .
  unzip *.zip

str="Unzipping g52"
echo $str
cd ../CFC-12_CCl2F2_g52
  unzip *.zip

str="Unzipping g53"
echo $str
cd ../CFC-13_CClF3_g53
  unzip *.zip

str="Unzipping g54"
echo $str
cd ../CFC-14_CF4_g54
#cp ../CF4.zip .
  unzip *.zip

str="Unzipping g55"
echo $str
cd ../HCFC-21_CHCl2F_g55
  unzip *.zip

str="Unzipping g56"
echo $str
cd ../HCFC22-CHClF2_g56
  unzip *.zip

str="Unzipping g57"
echo $str
cd ../C2Cl3F3-CFC-113_g57
  unzip *.zip

str="Unzipping g58"
echo $str
cd ../C2Cl2F4-CFC-114_g58
  unzip *.zip

str="Unzipping g59"
echo $str
cd ../C2ClF5-CFC-115_g59
  unzip *.zip

str="Unzipping g60"
echo $str
cd ../CCl4_g60
  unzip *.zip

str="Unzipping g61"
echo $str
cd ../ClONO2_g61
  unzip *.zip

str="Unzipping g62"
echo $str
cd ../N2O5_g62
cp ../N2O5.zip .
  unzip *.zip

str="Unzipping g63"
echo $str
cd ../HNO4_g63
  unzip *.zip

str="Unzipping g64"
echo $str
cd ../C2F6_g64
  unzip *.zip

str="Unzipping g65"
echo $str
cd ../CHCl2CF3-HCFC-123_g65
  unzip *.zip

str="Unzipping g66"
echo $str
cd ../CHCLFCF3-HCFC-124_g66
  unzip *.zip

str="Unzipping g67"
echo $str
cd ../CH3CCL2F-HCFC141b_g67
#### cp ../HCFC-141b.zip .
   unzip *.zip
####### /bin/cp -a ../../H2016/XSEC/CH3CCL2F-HCFC141b_g67/* .

str="Unzipping g68"
echo $str
cd ../CH3CCLF2-HCFC142b_g68
  unzip *.zip

str="Unzipping g69"
echo $str
cd ../CHCL2CF2CF3-HCFC-225ca_g69
  unzip *.zip

str="Unzipping g70"
echo $str
cd ../CCLF2CF2CHCLF-HCFC-225cb_g70
  unzip *.zip

str="Unzipping g71"
echo $str
cd ../CH2F2-HFC-32_g71
  unzip *.zip
   
str="Unzipping g72"
echo $str
cd ../CFH2CF3-HCF134a_g72
  unzip *.zip

str="Unzipping g73"
echo $str
cd ../CF3CH3-HFC-143a_g73
  unzip *.zip

str="Unzipping g74"
echo $str
cd ../CH3CHF2-HFC-152a_g74
  unzip *.zip

str="Unzipping g75"
echo $str
cd ../C6H6_g75
  unzip *.zip

str="Unzipping g76"
echo $str
cd ../CHF2CF3-HFC-125_g76
  unzip *.zip

str="Unzipping g77"
echo $str
cd ../CHF2CHF2-HFC-134_g77
  unzip *.zip

str="Unzipping g78"
echo $str
cd ../SF5CF3_g78
  unzip *.zip
   
str="Unzipping g79"
echo $str
cd ../CH3COOONO2-PAN_g79
  unzip *.zip

str="Unzipping g80"
echo $str
cd ../CH3CN_g80
  unzip *.zip

str="Unzipping g81"
echo $str
cd ../SF6_g81
  unzip *.zip

str="Done unzipping g51-81 or whichever ones were not commented"
echo $str
