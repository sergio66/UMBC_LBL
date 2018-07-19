# UMBC_LBL is a repo to compute mostly VoigtVanHuber lineshapes for the infrared. It also has code to do non-standard lineshapes for WV
and CO2 (though the CO2 linemixing is ancient, better to use eg LBLRTM for CO2 and CH4)
The code is mostly Matlab based though there are f77 mex files to speed up loops, and C based HITRAN line parameter readers
