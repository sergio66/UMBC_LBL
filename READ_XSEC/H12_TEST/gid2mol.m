function gstr = gid2mol(gid)

% function gstr = gid2mol(gid)
%
% gid2mol() translates an IR cross-section gas ID to the 
% molecular name from the HITRAN headers
%
% inputs
%   gid  - IR cross-section gas ID
%
% output
%   gstr - molecular name, from HITRAN record header

idd = 1998;   %% these are the old ones
idd = 2008;   %% these are the new ones

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
elseif idd > 1998
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
    %%%these are new
    case 64, gstr = 'C2F6_IR01';       % 64
    case 65, gstr = 'HCFC-123_IR00';   % 65
    case 66, gstr = 'HCFC-124_IR00';   % 66
    case 67, gstr = 'HCFC-141b_IR00';  % 67
    case 68, gstr = 'HCFC-142b_IR00';  % 68
    case 69, gstr = 'HCFC-225ca_IR00'; % 69
    case 70, gstr = 'HCFC-225cb_IR00'; % 70
    case 71, gstr = 'HFC-32_IR00';     % 71
    case 72, gstr = 'HFC-134a_IR00';   % 72
    %case 73, gstr = 'HFC-143a_IR00';   % 73
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
  end
