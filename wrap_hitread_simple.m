function line = wrap_hitread_simple(gid,start,stop,strengthM,vers);

if vers == 1992
  XHITRAN = '/asl/data/hitran/h92.by.gas';
  XHITRAN = [hitranpath '/h92.by.gas'];
elseif vers == 1996
  XHITRAN = '/asl/data/hitran/h96.by.gas';
  XHITRAN = [hitranpath '/h96.by.gas'];  
elseif vers == 1998
  XHITRAN = '/asl/data/hitran/h98.by.gas';
  XHITRAN = [hitranpath '/h98.by.gas'];  
elseif vers == 2000
  XHITRAN = '/asl/data/hitran/h2k.by.gas';
  XHITRAN = [hitranpath '/h2k.by.gas'];  
elseif vers == 2004
  XHITRAN = '/asl/data/hitran/h04.by.gas';
  XHITRAN = [hitranpath '/h04.by.gas'];  
elseif vers == 2008
  XHITRAN = '/asl/data/hitran/h08.by.gas';
  XHITRAN = [hitranpath '/h08.by.gas'];  
elseif vers == 2012
  XHITRAN = '/asl/data/hitran/h12.by.gas';
  XHITRAN = [hitranpath '/h12.by.gas'];  
elseif vers == 2016
  XHITRAN = '/asl/data/hitran/h16.by.gas';
  XHITRAN = [hitranpath '/h16.by.gas'];  
elseif vers == 2020
  XHITRAN = '/asl/data/hitran/h20.by.gas';
  XHITRAN = [hitranpath '/h20.by.gas'];  
elseif vers == 2024
  XHITRAN = '/asl/data/hitran/h24.by.gas';
  XHITRAN = [hitranpath '/h24.by.gas'];  
end

if (XHITRAN(length(XHITRAN)) == '/') 
  fnamePRE = [XHITRAN 'g' ]; 
else 
  fnamePRE = [XHITRAN '/g']; 
end 
fnamePOST = '.dat'; 
fnameIN   = int2str(gid); 
fname     = [fnamePRE fnameIN fnamePOST]; 

line = hitread_simple(start,stop,strengthM,gid,fname);
