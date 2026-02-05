#ls XSEC       gives the new files Chris downloaded
#-rw-rw-r-- 1 chepplew pi_strow    894708 Nov 18 10:09 N2O5.zip             g62
#-rw-rw-r-- 1 chepplew pi_strow  97303574 Nov 18 10:07 CF4.zip              g54
#-rw-rw-r-- 1 chepplew pi_strow 252739555 Nov 18 10:04 SF6.zip              g81
#-rw-rw-r-- 1 chepplew pi_strow 120517994 Nov 18 09:59 HCFC-141b.zip        g67  TIS IS IERD, does not have *TORR*
#-rw-rw-r-- 1 chepplew pi_strow 122215923 Nov 18 09:51 CFC-11.zip           g51

cd XSEC

#cp -a ../../H2016/XSEC/CFC-11_CCl3F_g51/*.xsc             CFC-11_CCl3F_g51/.
cp -a ../../H2016/XSEC/CFC-12_CCl2F2_g52/*.xsc            CFC-12_CCl2F2_g52/.
cp -a ../../H2016/XSEC/CFC-13_CClF3_g53/*.xsc             CFC-13_CClF3_g53/.
#cp -a ../../H2016/XSEC/CFC-14_CF4_g54/*.xsc               CFC-14_CF4_g54/.
cp -a ../../H2016/XSEC/HCFC-21_CHCl2F_g55/*.xsc           HCFC-21_CHCl2F_g55/.
cp -a ../../H2016/XSEC/HCFC22-CHClF2_g56/*.xsc            HCFC22-CHClF2_g56/.
cp -a ../../H2016/XSEC/C2Cl3F3-CFC-113_g57/*.xsc          C2Cl3F3-CFC-113_g57/.
cp -a ../../H2016/XSEC/C2Cl2F4-CFC-114_g58/*.xsc          C2Cl2F4-CFC-114_g58/.
cp -a ../../H2016/XSEC/C2ClF5-CFC-115_g59/*.xsc           C2ClF5-CFC-115_g59/.
cp -a ../../H2016/XSEC/CCl4_g60/*.xsc                     CCl4_g60/.
cp -a ../../H2016/XSEC/ClONO2_g61/*.xsc                   ClONO2_g61/.
#cp -a ../../H2016/XSEC/N2O5_g62/*.xsc                     N2O5_g62/.
cp -a ../../H2016/XSEC/HNO4_g63/*.xsc                     HNO4_g63/.
cp -a ../../H2016/XSEC/C2F6_g64/*.xsc                     C2F6_g64/.
cp -a ../../H2016/XSEC/CHCl2CF3-HCFC-123_g65/*.xsc        CHCl2CF3-HCFC-123_g65/.
cp -a ../../H2016/XSEC/CHCLFCF3-HCFC-124_g66/*.xsc        CHCLFCF3-HCFC-124_g66/.
#cp -a ../../H2016/XSEC/CH3CCL2F-HCFC141b_g67/*.xsc        CH3CCL2F-HCFC141b_g67/.      THIS IS WIERD, need to hand copy some files
cp -a ../../H2016/XSEC/CH3CCLF2-HCFC142b_g68/*.xsc        CH3CCLF2-HCFC142b_g68/.
cp -a ../../H2016/XSEC/CHCL2CF2CF3-HCFC-225ca_g69/*.xsc   CHCL2CF2CF3-HCFC-225ca_g69/.
cp -a ../../H2016/XSEC/CCLF2CF2CHCLF-HCFC-225cb_g70/*.xsc CCLF2CF2CHCLF-HCFC-225cb_g70/.
cp -a ../../H2016/XSEC/CH2F2-HFC-32_g71/*.xsc             CH2F2-HFC-32_g71/.
cp -a ../../H2016/XSEC/CFH2CF3-HCF134a_g72/*.xsc         CFH2CF3-HCF134a_g72/.
cp -a ../../H2016/XSEC/CF3CH3-HFC-143a_g73/*.xsc         CF3CH3-HFC-143a_g73/.
cp -a ../../H2016/XSEC/CH3CHF2-HFC-152a_g74/*.xsc        CH3CHF2-HFC-152a_g74/.
cp -a ../../H2016/XSEC/C6H6_g75/*.xsc                    C6H6_g75/.
cp -a ../../H2016/XSEC/CHF2CF3-HFC-125_g76/*.xsc         CHF2CF3-HFC-125_g76/.
cp -a ../../H2016/XSEC/CHF2CHF2-HFC-134_g77/*.xsc        CHF2CHF2-HFC-134_g77/.
cp -a ../../H2016/XSEC/SF5CF3_g78/*.xsc                  SF5CF3_g78/.
cp -a ../../H2016/XSEC/CH3COOONO2-PAN_g79/*.xsc          CH3COOONO2-PAN_g79/.
cp -a ../../H2016/XSEC/CH3CN_g80/*.xsc                   CH3CN_g80//
#cp -a ../../H2016/XSEC/SF6_g81/*.xsc                     SF6_g81/.
			
echo "now unzip the new files g51,g62, etc using unzipper_xsec.sc"
