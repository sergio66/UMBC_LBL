function gstr = gid2mol(gid,year)

% function gstr = gid2mol(gid)
%
% gid2mol() translates an IR cross-section gas ID to the 
% molecular name from the HITRAN headers
%
% inputs
%   gid - IR cross-section gas ID
%  year - HITRAN year (default 2012)
%
% output
%   gstr - molecular name, from HITRAN record header

if nargin == 1
  idd = 1998;   %% these are the H1998 ones
  idd = 2008;   %% these are the H2008 ones
  idd = 2012;   %% these are the H2012 ones
  idd = 2016;   %% these are the H2016 ones
  idd = 2020;   %% these are the H2020 ones
else
  idd = year;
end

if idd <= 1998
  switch gid
    % /asl/data/hitran/xsec98.ok/
    case 51, gstr = 'CCl3F';	% F11
    case 52, gstr = 'CCl2F2';	% F12
    case 53, gstr = 'CClF3';	% F13
    case 54, gstr = 'CF4';	% F14
    case 55, gstr = 'CHCl2F';	% F21
    case 56, gstr = 'CHClF2';	% F22
    case 57, gstr = 'C2Cl3F3';	% F113
    case 58, gstr = 'C2Cl2F4';	% F114
    case 59, gstr = 'C2ClF5';	% F115
    case 60, gstr = 'CCl4';	% carbon tetrachloride
    case 61, gstr = 'ClONO2';	% chlorine nitrate
    case 62, gstr = 'N2O5';	% dinitrogen pentoxide
    case 63, gstr = 'HNO4';	% pernitric acid
    % gases added from HITRAN 98
    case 64, gstr = 'SF6';	% sulphur hexafluoride (?)
    case 65, gstr = 'CH2F2';	% HFC-32
    case 66, gstr = 'CHF2CHF2';	% HFC-134
    case 67, gstr = 'CF3CH3';	% HFC-143a
    otherwise, error('unknown gas id')
  end
elseif idd >= 1998 & idd < 2012
  % /asl/data/hitran/HITRAN08_SERGIO/Xsec/
  switch gid
    case 51, gstr = 'CFC-11_IR00';    % 51  (F11)
    case 52, gstr = 'CFC-12_IR00';    % 52  (F12)
    case 53, gstr = 'CFC-13_IR01';    % 53  (F13)
    case 54, gstr = 'CFC-14_IR01';    % 54  (F14)    also gid 42
    case 55, gstr = 'HCFC-21_IR00';   % 55  (F21)
    case 56, gstr = 'HCFC-22_IR01';   % 56  (F22)
    case 57, gstr = 'CFC-113_IR00';   % 57  (F113)
    case 58, gstr = 'CFC-114_IR00';   % 58  (F114)
    case 59, gstr = 'CFC-115_IR00';   % 59  (F115)
    case 60, gstr = 'CCl4_IR00';     % 60
    case 61, gstr = 'ClONO2_IR04';   % 61            also gid 35
    case 62, gstr = 'N2O5_IR04';     % 62
    case 63, gstr = 'HNO4_IR04';     % 63

    %%%these are new in H2008
    case 64, gstr = 'C2F6_IR01';       % 64
    case 65, gstr = 'HCFC-123_IR00';   % 65
    case 66, gstr = 'HCFC-124_IR00';   % 66
    case 67, gstr = 'HCFC-141b_IR00';  % 67
    case 68, gstr = 'HCFC-142b_IR00';  % 68
    case 69, gstr = 'HCFC-225ca_IR00'; % 69
    case 70, gstr = 'HCFC-225cb_IR00'; % 70
    case 71, gstr = 'HFC-32_IR00';     % 71
    case 72, gstr = 'HFC-134a_IR00';   % 72
    % case 73, gstr = 'HFC-143a_IR00';   % 73
    case 73, gstr = 'HFC-143a_IR08';   % 73
    case 74, gstr = 'HFC-152a_IR00';   % 74
    case 75, gstr = 'C6H6_IR08';       % 75
    case 76, gstr = 'HFC-125_IR08';    % 76
    case 77, gstr = 'HFC-134_IR00';    % 77
    case 78, gstr = 'SF5CF3_IR04';     % 78
    case 79, gstr = 'PAN_IR04';        % 79
    case 80, gstr = 'CH3CN_IR05';      % 80       also gid 41
    case 81, gstr = 'SF6_IR00';        % 81       also gid 30
  end

elseif idd >= 2012
  % /asl/data/hitran/HITRAN2012/IR-XSect/Uncompressed-files
  switch gid
    case 51, gstr = 'CFC-11_IR00';    % 51  (F11)
    case 52, gstr = 'CFC-12_IR00';    % 52  (F12)
    case 53, gstr = 'CFC-13_IR01';    % 53  (F13)
    case 54, gstr = 'CFC-14_IR01';    % 54  (F14)    also gid 42
    case 55, gstr = 'HCFC-21_IR00';   % 55  (F21)
    case 56, gstr = 'HCFC-22_IR01';   % 56  (F22)
    case 57, gstr = 'CFC-113_IR00';   % 57  (F113)
    case 58, gstr = 'CFC-114_IR00';   % 58  (F114)
    case 59, gstr = 'CFC-115_IR00';   % 59  (F115)
    case 60, gstr = 'CCl4_IR00';     % 60
    case 61, gstr = 'ClONO2_IR04';   % 61            also gid 35
    case 62, gstr = 'N2O5_IR04';     % 62
    case 63, gstr = 'HNO4_IR04';     % 63

    %%%these are new in H2008
    case 64, gstr = 'C2F6_IR01';       % 64
    case 65, gstr = 'HCFC-123_IR00';   % 65
    case 66, gstr = 'HCFC-124_IR00';   % 66
    case 67, gstr = 'HCFC-141b_IR11';  % 67
    case 68, gstr = 'HCFC-142b_IR11';  % 68
    case 69, gstr = 'HCFC-225ca_IR00'; % 69
    case 70, gstr = 'HCFC-225cb_IR00'; % 70
    case 71, gstr = 'HFC-32_IR00';     % 71
    case 72, gstr = 'HFC-134a_IR00';   % 72
    case 73, gstr = 'HFC-143a_IR00';   % 73
    case 74, gstr = 'HFC-152a_IR00';   % 74
    case 75, gstr = 'C6H6_IR08';       % 75
    case 76, gstr = 'HFC-125_IR08';    % 76
    case 77, gstr = 'HFC-134_IR00';    % 77
    case 78, gstr = 'SF5CF3_IR04';     % 78
    case 79, gstr = 'PAN_IR11';        % 79
    case 80, gstr = 'CH3CN_IR12';      % 80       also gid 41
    case 81, gstr = 'SF6_IR00';        % 81       also gid 30

    %% these are new in H2012, probably not used in kCARTA since no profile!
    case 82, gstr = 'C2H6_IR10';      % 82 x
    case 83, gstr = 'CH3OH_IR12';     % 83 x
    case 84, gstr = 'CH3CHO_IR11';    % 84 x
    case 85, gstr = 'C3H8_IR10';      % 85
    case 86, gstr = 'CH3COCH3_IR11';  % 86
    case 87, gstr = 'BrONO2_IR12';    % 87
    case 88, gstr = 'ClOOCl_IR12';    % 88
  end

  if gid >= 51 & gid <= 81 & idd >= 2015    
    ystr = num2str(idd);
    gstr = [gstr(1:end-2) ystr];
  end
    
end

