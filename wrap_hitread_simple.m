function line = wrap_hitread_simple(gid,start,stop,strengthM,vers);

if vers == 1992
  HITRAN = '/asl/data/hitran/h92.by.gas';
elseif vers == 1996
  HITRAN = '/asl/data/hitran/h96.by.gas';
elseif vers == 1998
  HITRAN = '/asl/data/hitran/h98.by.gas';
elseif vers == 2000
  HITRAN = '/asl/data/hitran/h2k.by.gas';
elseif vers == 2004
  HITRAN = '/asl/data/hitran/h04.by.gas';
elseif vers == 2008
  HITRAN = '/asl/data/hitran/h08.by.gas';
  end

if (HITRAN(length(HITRAN)) == '/') 
  fnamePRE = [HITRAN 'g' ]; 
else 
  fnamePRE = [HITRAN '/g']; 
  end 
fnamePOST = '.dat'; 
fnameIN   = int2str(gid); 
fname     = [fnamePRE fnameIN fnamePOST]; 

line = hitread_simple(start,stop,strengthM,gid,fname);