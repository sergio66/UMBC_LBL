function [grr]=efitgrad(bds,jr,elowerr);    
  
grr=zeros(length(jr),3);  
  
grr(:,1)=1.0;  
grr(:,2)=(jr.*(jr+1));  
grr(:,3)=-grr(:,2).^2;  
 
grr=grr'; 
