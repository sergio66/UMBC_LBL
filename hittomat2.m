function [line]=hittomat2(low,high,vers,strength,gasID)
%function [line]=hittomat2(low,high,vers,strength,gasID)
%this is a function version of  script file hittomat.m
%low   = low freq bound
%high  =  high freq bound
%vers  = line version number
%str   = minimum line strength
%gasID = gasID

WNUMLO=[];
WNUMHI=[];
IDG=[];
ICH=[];
STR=[];

WNUMLO(1)=low;
WNUMHI(1)=high;
IDG(1)=gasID;
ICH(1)=vers;
STR(1)=strength;

%ISUM=input('enter the number of gases to be loaded:');
ISUM=1;

%for i=1:ISUM
%  IDG(i)=input('Enter gas id:');
%  WNUMLO(i)=input('Enter minimum frequency:');
%  WNUMHI(i)=input('Enter maximum frequency:');
%  ICH(i)=input('Enter line version #:');
%  STR(i)=input('Enter minimum strength:');
%  end

%INFILE=input('Enter name of hitran linebase file:(hitlin92.bin)','s');
%if isempty(INFILE)
%  INFILE='hitlin92.bin';
%  end
%INFILE='Hitlin/hitlin92.bin';
%INFILE='Hitlin/hitlin_fewlines.bin';
INFILE='Hitlin/hitlin96.bin';


[line.LINCT,line.WNUMHX,line.ZLSTAT,line.ZIGAS,line.ZISO,line.ZWNUM, ...
 line.ZSTREN,line.ZTPROB,line.ZABROAD,line.ZSBROAD,line.ZELS,....
 line.ZABCOEF,line.ZTSP,line.ZIUSGQ,line.ZILSGQ,line.ZUSLQ,...
 line.ZBSLQ,line.ZAI,line.ZREF,line.GASID]=  ...
      mexhitd(WNUMLO,WNUMHI,INFILE,IDG,ICH,STR,ISUM);

line.ZAI=char(line.ZAI)';
line.ZBSLQ=char(line.ZBSLQ)';
line.ZREF=char(line.ZREF)';
line.ZUSLQ=char(line.ZUSLQ)';

