function [jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq]=...
     orderer(jq,elowerq,w_forq,w_selfq,freqq,strenqt,strenq)

[jq,ii] = sort(jq);
elowerq = elowerq(ii);
w_forq  = w_forq(ii);
w_selfq = w_selfq(ii);
freqq   = freqq(ii);
strenqt = strenqt(ii);
strenq  = strenq(ii); 

%format short e
%[jq(1:9) freqq(1:9) strenqt(1:9) w_forq(1:9) w_selfq(1:9)]

%%% ind = find(jq <=  81);
%%% jq = jq(ind);freqq= freqq(ind);strenqt= strenqt(ind);elowerq= elowerq(ind);
%%% w_forq= w_forq(ind);w_selfq= w_selfq(ind);strenq= strenq(ind);

