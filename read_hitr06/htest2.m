
% compare matlab and mexized C HITRAN readers

more off

hdat = '/asl/data/hitran/h04.by.gas';
v1 = 600;
v2 = 2800;

gid = 12

fprintf(1, 'reading data for gas %d...\n', gid)
tic; L1 = read_hitran2(v1, v2, 0, gid, hdat); toc
tic; L2 = read_hitran(v1, v2, 0, gid, hdat); toc

isequal(L1, L2)

isequal(L1.igas, L2.igas)
isequal(L1.iso, L2.iso)
isequal(L1.wnum, L2.wnum)
isequal(L1.stren, L2.stren)
isequal(L1.tprob, L2.tprob)
isequal(L1.abroad, L2.abroad)
isequal(L1.sbroad, L2.sbroad)
isequal(L1.els, L2.els)
isequal(L1.abcoef, L2.abcoef)
isequal(L1.tsp, L2.tsp)
isequal(L1.iusgq, L2.iusgq)
isequal(L1.ilsgq, L2.ilsgq)
isequal(L1.uslq, L2.uslq)
isequal(L1.bslq, L2.bslq)
isequal(L1.ai, L2.ai)
isequal(L1.ref, L2.ref)
isequal(L1.flag, L2.flag)
isequal(L1.swus, L2.swus)
isequal(L1.swls, L2.swls)

