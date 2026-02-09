error('uou need to fix the asl/data/hitran/ to eg umbc/xfs3/strow/asl/rta/hitran/, see hitranpath.m')

topts.HITRAN = '/asl/data/hitran/h98.by.gas';
[w,d8_98] = run8(9,1380,1405,'IPFILES/std_n2o_lay1_3',topts);

topts.HITRAN = '/asl/data/hitran/h04.by.gas';
[w,d8_h04] = run8(9,1380,1405,'IPFILES/std_n2o_lay1_3',topts);

topts.HITRAN = '/asl/data/hitran/h08.by.gas';
[w,d8_h08] = run8(9,1380,1405,'IPFILES/std_n2o_lay1_3',topts);

topts.HITRAN = '/asl/data/hitran/h98.by.gas';
[w,d7_98] = run7(9,1380,1405,'IPFILES/std_n2o_lay1_3',topts);
