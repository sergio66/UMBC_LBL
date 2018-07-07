function [fullmix4]=full2(freqq,f,W_plus,jq,w_selfq,w_forq,temperature,...
   trans_ampl,population_t,stuff,birn,ratio,strenqt,layeramt,ymix1,band,prb) 

if (length(intersect(band,[668 740 2093 2322])) > 0) %%orig, before Aug 2007
%%if (length(intersect(band,[668 2093 2322])) > 0) %%%testing 740 branch
  fullmix4=full2_deltpi(freqq,f,W_plus,jq,w_selfq,w_forq,temperature,...
   trans_ampl,population_t,stuff,birn,ratio,strenqt,layeramt,ymix1,band,prb);
else
  fullmix4=full2_usual(freqq,f,W_plus,jq,w_selfq,w_forq,temperature,...
   trans_ampl,population_t,stuff,birn,ratio,strenqt,layeramt,ymix1,band,prb);
  end

