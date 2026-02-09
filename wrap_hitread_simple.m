function line = wrap_hitread_simple(gid,start,stop,strengthM,vers);

hitranpath
if vers == 1992
  HITRANX = '/asl/data/hitran/h92.by.gas';
  HITRANX = [HITRAN '/h92.by.gas'];
elseif vers == 1996
  HITRANX = '/asl/data/hitran/h96.by.gas';
  HITRANX = [HITRAN '/h96.by.gas'];  
elseif vers == 1998
  HITRANX = '/asl/data/hitran/h98.by.gas';
  HITRANX = [HITRAN '/h98.by.gas'];  
elseif vers == 2000
  HITRANX = '/asl/data/hitran/h2k.by.gas';
  HITRANX = [HITRAN '/h2k.by.gas'];  
elseif vers == 2004
  HITRANX = '/asl/data/hitran/h04.by.gas';
  HITRANX = [HITRAN '/h04.by.gas'];  
elseif vers == 2008
  HITRANX = '/asl/data/hitran/h08.by.gas';
  HITRANX = [HITRAN '/h08.by.gas'];  
elseif vers == 2012
  HITRANX = '/asl/data/hitran/h12.by.gas';
  HITRANX = [HITRAN '/h12.by.gas'];  
elseif vers == 2016
  HITRANX = '/asl/data/hitran/h16.by.gas';
  HITRANX = [HITRAN '/h16.by.gas'];  
elseif vers == 2020
  HITRANX = '/asl/data/hitran/h20.by.gas';
  HITRANX = [HITRAN '/h20.by.gas'];  
elseif vers == 2024
  HITRANX = '/asl/data/hitran/h24.by.gas';
  HITRANX = [HITRAN '/h24.by.gas'];  
end

if (HITRANX(length(HITRANX)) == '/') 
  fnamePRE = [HITRANX 'g' ]; 
else 
  fnamePRE = [HITRANX '/g']; 
end 
fnamePOST = '.dat'; 
fnameIN   = int2str(gid); 
fname     = [fnamePRE fnameIN fnamePOST]; 

line = hitread_simple(start,stop,strengthM,gid,fname);
