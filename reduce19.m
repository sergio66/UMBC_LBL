function lineNEW=reduce19(lineORIG);

%% function lineNEW=reduce19(lineORIG);
%% since GAS1D == 19 is problematic (H98 database claims there are 5 isotopes
%% while supplied BD_*.FOR code says there are 4 isotopes, let us remove the
%% fifth isotope from lineORIG

fprintf(1,'gasID 19 is not correctly entered into H98 database .. pruning \n');

index=find(lineORIG.iso ~= 5);
lineNEW.igas=lineORIG.igas;

lineNEW.iso=lineORIG.iso(index);

lineNEW.wnum=lineORIG.wnum(index);
lineNEW.stren=lineORIG.stren(index);
lineNEW.tprob=lineORIG.tprob(index);

lineNEW.abroad=lineORIG.abroad(index);
lineNEW.sbroad=lineORIG.sbroad(index);

lineNEW.els=lineORIG.els(index);
lineNEW.abcoef=lineORIG.abcoef(index);
lineNEW.tsp=lineORIG.tsp(index);

lineNEW.iusgq=lineORIG.iusgq(index);
lineNEW.ilsgq=lineORIG.ilsgq(index);

lineNEW.uslq=lineORIG.uslq(index);
lineNEW.bslq=lineORIG.bslq(index);

lineNEW.ai=lineORIG.ai(index);
lineNEW.ref=lineORIG.ref(index);

lineNEW.linct=length(index);
lineNEW.gasid=lineORIG.gasid(index);

