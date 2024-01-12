function [line,hitran_version,hlist_qtips] = ...
            hitread(start,stop,strengthM,gasID,filename);

%%    hlist_qtips tells you which HXX versions use the old qtips, and which
%%    use the new Lagrange stuff from the H04 database
% function line=hitreadNEW(start,stop,strengthM,gasID,filename);
% this calls read_hitran, and adds on an extra field : LINCT which is
% the number of lines read in

current_dir = pwd;

%{
cd ../SAMPLE_HITREAD
lineMatlab = read_hitran2(start,stop,strengthM,gasID,filename); line = lineMatlab;
%cd ../read_hitr06/
which read_hitran
lineC = read_hitran(start,stop,strengthM,gasID,filename);  line = lineC;
cd ../SAMPLE_HITREAD
%}

lineMatlab = read_hitran2(start,stop,strengthM,gasID,filename); line = lineMatlab;
%which read_hitran
lineC = read_hitran(start,stop,strengthM,gasID,filename);  line = lineC;

%sum(lineC.wnum-lineMatlab.wnum)
%sum(lineC.stren-lineMatlab.stren)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
line.linct = length(line.wnum);

line.iso   = line.iso';
line.wnum  = line.wnum';
line.stren = line.stren';
line.tprob = line.tprob';

line.abroad = line.abroad';
line.sbroad = line.sbroad';

line.els    = line.els';
line.abcoef = line.abcoef';
line.tsp    = line.tsp';

line.iusgq = line.iusgq';
line.ilsgq = line.ilsgq';

line.gasid = line.igas';
