topts.HITRAN = '/asl/data/hitran/h98.by.gas';
[w,d8_98] = run8(4,1580,1605,'IPFILES/std_n2o_lay1_3',topts);

topts.HITRAN = '/asl/data/hitran/h04.by.gas';
[w,d8_h04] = run8(4,1580,1605,'IPFILES/std_n2o_lay1_3',topts);

topts.HITRAN = '/asl/data/hitran/h08.by.gas';
[w,d8_h08] = run8(4,1580,1605,'IPFILES/std_n2o_lay1_3',topts);

topts.HITRAN = '/asl/data/hitran/h98.by.gas';
[w,d7_98] = run7(4,1580,1605,'IPFILES/std_n2o_lay1_3',topts);
