:

echo 'cd ../CO_PR_sigpi/'
cd ../CO_PR_sigpi;  rm *~

rm full2.m          orderer.m        voigtmixRatio.m
rm trans_pop.m      tempratio*.m     co2_param.m
rm efit.m           voigtmix2.m      y1s.m
rm smooth_full_ratio.m  doratio.m    qqttiippss.m
rm removeneg.m          blender4um.m
rm IncludeMaxer.m
rm full2_usual.m
rm co2_param_JOHNS2002.m ral_jjohn_blend3.m ral_jjohn_blend.m small_leap.m
rm find_regions.m
rm fudger_and_threeregions.m
rm efit1.m
ln -s ../CO2_COMMON/efit1.m              efit1.m
ln -s ../CO2_COMMON/fudger_and_threeregions.m fudger_and_threeregions.m
ln -s ../CO2_COMMON/blender4um.m           blender4um.m
ln -s ../CO2_COMMON/removeneg.m            removeneg.m    
ln -s ../CO2_COMMON/tempratio4um.m         tempratio4um.m
ln -s ../CO2_COMMON/co2_param.m            co2_param.m
ln -s ../CO2_COMMON/doratio.m              doratio.m
ln -s ../CO2_COMMON/efit.m                 efit.m
ln -s ../CO2_COMMON/full2.m                full2.m
ln -s ../CO2_COMMON/full2_usual.m          full2_usual.m
ln -s ../CO2_COMMON/orderer.m              orderer.m
ln -s ../CO2_COMMON/qqttiippss.m           qqttiippss.m
ln -s ../CO2_COMMON/smooth_full_ratio.m    smooth_full_ratio.m 
ln -s ../CO2_COMMON/trans_pop.m            trans_pop.m
ln -s ../CO2_COMMON/voigtmix2.m            voigtmix2.m
ln -s ../CO2_COMMON/voigtmixRatio.m        voigtmixRatio.m
ln -s ../CO2_COMMON/y1s.m                  y1s.m
ln -s ../IncludeMaxer.m                    IncludeMaxer.m 
ln -s ../CO2_COMMON/co2_param_JOHNS2002.m  co2_param_JOHNS2002.m 
ln -s ../CO2_COMMON/small_leap.m           small_leap.m
ln -s ../CO2_COMMON/ral_jjohn_blend3.m     ral_jjohn_blend3.m
ln -s ../CO2_COMMON/ral_jjohn_blend.m      ral_jjohn_blend.m
ln -s ../CO2_COMMON/find_regions.m         find_regions.m

rm driver*.m
ln -s ../CO2_COMMON/MAIN4umDRIVERS/driver4um.m  driver4um.m

cd ..;  rm *~;
cd CO2_COMMON

echo 'done CO!!'

########################################################################
########################################################################
########################################################################
