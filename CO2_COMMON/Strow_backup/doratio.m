function [ratio]=doratio(population_t,trans_ampl,W_plus);
%this function finds the ratio kmix/klor at v-v0 >>>> v0
%this is from Larrabee Strow's paper

global quiet

rd=population_t.*trans_ampl;  
len=length(rd);  
ratio=W_plus.*(ones(len,1)*rd').*(trans_ampl*ones(1,len));  
ratio=triu(ratio,1);  
denom=diag(W_plus).*rd.*trans_ampl;  
denom=sum(denom);  
ratio=1+2*sum(sum(ratio))/denom;

if (ratio <= 0.0) 
  ratio=1.0e-10;         %if the sum rule give -ve numbers
  end

%%ratio = 0.25;
fprintf(1,'ratio of kfull/klor = %9.6f \n',ratio)
