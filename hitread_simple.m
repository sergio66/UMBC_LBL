function [line,hitran_version] = hitread_no_symbolic_links(...
                                   start,stop,strengthM,gasID,filename);

%% same as hitread.m except it does NOT update symbolic links


% function line=hitreadNEW(start,stop,strengthM,gasID);
% this calls read_hitran, and adds on an extra field : LINCT which is
% the number of lines read in

% it differentiates and calls the right MEXED HITRAN readers

current_dir = pwd;

blah = findstr('/',filename);
%%look at the fourth occurence
hitran_version = filename(blah(4)+1:blah(4)+3);

old_reader = {'h92','h96','h98','h2k'};
new_reader = {'h04','h08','h12','h16'};

iOld = -1;
iNew = -1;
for ii = 1 : length(old_reader)
  if strcmp(old_reader(ii),hitran_version)
    iOld = +1;
    end
  end
for ii = 1 : length(new_reader)
  if strcmp(new_reader(ii),hitran_version)
    iNew = +1;
    end
  end
if iOld == -1 & iNew == -1
  error('whoops : have NOT been able to figure what HITRAN version you want!')
  end

if iOld > 0
  fprintf(1,'using old MEX reader for HITRAN version %s \n',hitran_version)
  %cd /asl/matlab2012/aslutil/
  cd /asl/matlib/aslutil/  
  line = read_hitran(start,stop,strengthM,gasID,filename);
else
  fprintf(1,'using new MEX reader for HITRAN version %s \n',hitran_version)
  %cd /asl/matlab2012/read_hitr06
  cd /home/sergio/SPECTRA/read_hitr06
  line = read_hitran(start,stop,strengthM,gasID,filename);

end

%iVers = input('Enter HITRAN version 92 96 98 2k 04 08 12 16: ');
%database = ['/asl/data/hitran/h' iVers '.by.gas'];
%if iOld > 0
%  fprintf(1,'using old MEX reader for HITRAN version %s \n',hitran_version)
%  cd /asl/matlab/aslutil/
%  line = read_hitran(start,stop,strengthM,gasID,filename,database);
%else
%  fprintf(1,'using new MEX reader for HITRAN version %s \n',hitran_version)
%  cd /asl/matlab/read_hitr06
%  line = read_hitran(start,stop,strengthM,gasID,filename,database);
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
line.linct=length(line.wnum);

line.iso=line.iso';
line.wnum=line.wnum';
line.stren=line.stren';
line.tprob=line.tprob';

line.abroad=line.abroad';
line.sbroad=line.sbroad';

line.els=line.els';
line.abcoef=line.abcoef';
line.tsp=line.tsp';

line.iusgq=line.iusgq';
line.ilsgq=line.ilsgq';

line.gasid=line.igas';

if line.linct > 0
  line.igas=line.gasid(1);
  end

cder = ['cd ' current_dir];
eval(cder)

