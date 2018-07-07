%this is a script to compare HITRAN databases

gasID = input('Enter gas : ');

low  = 500;
high = 3000;
strengthM = 0.0;

%% fnamePRE='/salsify/scratch4/h96.by.gas/g';        %H96 -- old
%% fnamePRE='/salsify/scratch4/h98.by.gas/g';        %H98 -- KCARTA database 
%% fnamePRE='/salsify/scratch4/h2k.by.gas/g';        %H98 -- KCARTA database 

HITRAN        = '/asl/data/hitran/h2k.by.gas';
%%%%%HITRAN        = '/taro/s1/sergio/sergio.h2k';
if (HITRAN(length(HITRAN)) == '/')
     fnamePRE = [HITRAN 'g' ];
     else
fnamePRE = [HITRAN '/g'];
  end
fnamePOST='.dat';
fnameIN=int2str(gasID);
fname=[fnamePRE fnameIN fnamePOST];
[line2000]=hitread(low,high,strengthM,gasID,fname);

HITRAN        = '/asl/data/hitran/h98.by.gas';
if (HITRAN(length(HITRAN)) == '/')
     fnamePRE = [HITRAN 'g' ];
     else
fnamePRE = [HITRAN '/g'];
  end
fnamePOST='.dat';
fnameIN=int2str(gasID);
fname=[fnamePRE fnameIN fnamePOST];
[line1998]=hitread(low,high,strengthM,gasID,fname);

HITRAN        = '/asl/data/hitran/h96.by.gas';
if (HITRAN(length(HITRAN)) == '/')
     fnamePRE = [HITRAN 'g' ];
     else
fnamePRE = [HITRAN '/g'];
  end
fnamePOST='.dat';
fnameIN=int2str(gasID);
fname=[fnamePRE fnameIN fnamePOST];
[line1996]=hitread(low,high,strengthM,gasID,fname);

if line2000.linct == line1998.linct
  fprintf(1,'same number of lines ... plotting ratios ... \n');
  plot(line2000.wnum./line1998.wnum);       title('wnum'); pause
  plot(line2000.iso./line1998.iso);         title('iso'); pause
  plot(line2000.stren./line1998.stren);     title('stren'); pause
  plot(line2000.tprob./line1998.tprob);     title('tprob'); pause
  plot(line2000.abroad./line1998.abroad);   title('abroad'); pause
  plot(line2000.sbroad./line1998.sbroad);   title('sbroad'); pause
  plot(line2000.els./line1998.els);         title('els'); pause
  plot(line2000.abcoef./line1998.abcoef);   title('abscof'); pause
  plot(line2000.tsp./line1998.tsp);         title('tsp'); pause
  plot(line2000.iusgq./line1998.iusgq);     title('iusgq'); pause
  plot(line2000.ilsgq./line1998.ilsgq);     title('ilsgq'); pause
else
  fprintf(1,'different number of lines ... \n');
  [C,IA,IB] = intersect(line1998.wnum,line2000.wnum) ;
semilogy(line2000.wnum(IB),line2000.stren(IB)./line1998.stren(IA),'r'); pause
semilogy(line2000.wnum(IB),line2000.abroad(IB)./line1998.abroad(IA),'r'); pause
semilogy(line2000.wnum(IB),line2000.sbroad(IB)./line1998.sbroad(IA),'r'); pause
  end
