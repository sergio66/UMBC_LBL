function gname = gid_to_gname(gid);

%% see pg 31 of https://www.cfa.harvard.edu/hitran/  H2012 documentation
%{
%% these are in LBLRTM XSEC, cannot find these in UMBC XSEC
HNO3
NO2
ACET
CCL2FCH3
CHCL2C2F5
SO2
ISOP
CHF3
%}

%{
%% these are in UMBC XSEC, cannot find in LBLRTM XSEC so remapping some IDs
64 : C2F6
66 : CHCLFCF3 -------- HCFC 124 ------> CHCl2C2F5
67 : CH3CCL2F -------- HCFC-141b -----> CCl2FCH3
68 : CH3CCLF2 -------- HCFC-142b -----> CH3CClF2
69 : CHCL2CF2CF3 ----- HCFC-225ca ----> C3HCl2F5
70 : CCLF2CF2CHCLF --- HCFC-225cb ----> C3HCl2F5
75 : C6F6 (benzene)????
77 : CHF2CHF2
78 : SF5CF3
80 : CH3CN
81 : SF6
%}

%% see /home/sergio/KCARTA/DOC/gasids_H2012

if gid == 51
  gname = ' CCL3F';
  gname = ' F11';
elseif gid == 52
  gname = ' CCL2F2';
  gname = ' F12';  
elseif gid == 53
  gname = ' CCLF3';
elseif gid == 54
  gname = ' CF4';
elseif gid == 55
  gname = ' CHCL2F';
elseif gid == 56
  gname = ' CHCLF2';
elseif gid == 57
  gname = ' C2CL3F3';
elseif gid == 58
  gname = ' C2CL2F4';
elseif gid == 59
  gname = ' C2CLF5';
elseif gid == 60
  gname = ' CCL4';
elseif gid == 61
  gname = ' CLONO2';
elseif gid == 62
  gname = ' N2O5';
elseif gid == 63
  gname = ' HNO4';
%%% 64-74 requested by Evan Fishbein, 28 June 2010
elseif gid == 64
  gname = ' C2F6'; %          hexafluoroethane
elseif gid == 65
  gname = ' CHCl2CF3'; %      HCFC-123  
elseif gid ==  66
  gname = ' CHCLFCF3';  %      HCFC-124
  gname = ' CHCl2C2F5'; %      HCFC-124
elseif gid ==  67
  gname = ' CH3CCL2F'; %      HCFC-141b
  gname = ' CCl2FCH3'; %      HCFC-141b
elseif gid == 68
  gname = ' CH3CCLF2'; %      HCFC-142b
  gname = ' CH3CClF2'; %      HCFC-142b  
elseif gid == 69
  gname = ' CHCL2CF2CF3'; %   HCFC-225ca
  gname = ' C3HCl2F5';    %   HCFC-225ca  
elseif gid == 70
  gname = ' CCLF2CF2CHCLF'; %  HCFC-225cb
  gname = ' C3HCl2F5';      %  HCFC-225cb  
elseif gid ==  71
  gname = ' CH2F2'; %          HFC-32
elseif gid ==  72
  gname = ' CFH2CF3';    %    HFC-134a
elseif gid == 73
  gname = ' CF3CH3'; % HFC-143a
elseif gid == 74
  gname = ' CH3CHF2'; %       HFC-152a
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
%%% added from new HITRAN gases
elseif gid == 75
  gname = '  C6F6'; %          benzene
elseif gid ==  76
  gname = ' CHF2CF3'; %       HFC-125
elseif gid == 77
  gname = ' CHF2CHF2'; %      HFC-134   
elseif gid == 78
  gname = ' SF5CF3';
elseif gid == 79
  gname = 'CH3C(O)OONO2'; %  PAN % (Peroxy Acetyl Nitrate)
  gname = ' PAN';
elseif gid == 80
  gname = 'CH3CN'; % [also see 41]
elseif gid == 81
  gname = ' SF6';  % [also see 30]
end

if length(gname) < 10
  len = length(gname);
  blank = [];
  for ii = 1 : 10-len
    blank = [blank ' '];
  end
  gname = [gname blank];
end
