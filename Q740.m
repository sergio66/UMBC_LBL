load Q740.dat;
d=Q740;
%%%%temp=d(:,1); p=(d(:,2)+d(:,3))/1013; pp=d(:,2)/1013; 
temp=d(:,1); p=(d(:,2)+d(:,3))/760; pp=d(:,2)/760; 

MGC=8.314674269981136  ;  
q=d(:,4); 
q=q*101325.*pp/1e9/MGC./temp;

gid=2*ones(length(q),1);

neww=[gid p pp temp q];

fid=fopen('IPFILES/q740','w');
fprintf(fid,'%3i %8.5e %8.5e %8.5f %8.5e \n',neww');
fclose(fid)
