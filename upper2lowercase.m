function line2=upper2lowercase(line,gasID)
%function line2=upper2lowercase(line) changes CASE of fields

if (line.LINCT > 0)
  line2.igas    = line.GASID(1);
  line2.iso     = line.ZISO;
  line2.wnum    = line.ZWNUM;
  line2.stren   = line.ZSTREN;
  line2.tprob   = line.ZTPROB;
  line2.abroad  = line.ZABROAD;
  line2.sbroad  = line.ZSBROAD;
  line2.els     = line.ZELS;
  line2.abcoef  = line.ZABCOEF;
  line2.tsp     = line.ZTSP;
  line2.iusgq   = line.ZIUSGQ;
  line2.ilsgq   = line.ZILSGQ;
  line2.uslq    = line.ZUSLQ;
  line2.bslq    = line.ZBSLQ;
  line2.ai      = line.ZAI;
  line2.ref     = line.ZREF;
  line2.linct   = line.LINCT;
  line2.gasid   = line.GASID;
else
  line2.igas    = gasID;
  line2.iso     = line.ZISO;
  line2.wnum    = line.ZWNUM;
  line2.stren   = line.ZSTREN;
  line2.tprob   = line.ZTPROB;
  line2.abroad  = line.ZABROAD;
  line2.sbroad  = line.ZSBROAD;
  line2.els     = line.ZELS;
  line2.abcoef  = line.ZABCOEF;
  line2.tsp     = line.ZTSP;
  line2.iusgq   = line.ZIUSGQ;
  line2.ilsgq   = line.ZILSGQ;
  line2.uslq    = line.ZUSLQ;
  line2.bslq    = line.ZBSLQ;
  line2.ai      = line.ZAI;
  line2.ref     = line.ZREF;
  line2.linct   = line.LINCT;
  line2.gasid   = line.GASID;
  end
clear line
