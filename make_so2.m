topts.HITRAN = [hitranpath '/h98.by.gas/'];
[w,d8_98] = run8(9,1380,1405,'IPFILES/std_n2o_lay1_3',topts);

topts.HITRAN = [hitranpath '/h04.by.gas/'];
[w,d8_h04] = run8(9,1380,1405,'IPFILES/std_n2o_lay1_3',topts);

topts.HITRAN = [hitranpath '/h08.by.gas/'];
[w,d8_h08] = run8(9,1380,1405,'IPFILES/std_n2o_lay1_3',topts);

% run 7.m
% topts.HITRAN = [hitranpath '/h98.by.gas/'];
% [w,d7_98] = run7(9,1380,1405,'IPFILES/std_n2o_lay1_3',topts);
