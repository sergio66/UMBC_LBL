function [fr,stren] = only_pr(band);

%%% this ensures we only pick out line params for PR bands .. essentially for 
%%% the NLTE kCARTA code, else we double add on the Q lines (they have ALREADY
%%% been included in  caaLTEWeakBack_StdOptDepth

fnameIN = ['hit' num2str(band) '.mat'];
loader = ['load ' fnameIN];
eval(loader);

fprintf(1,'lower, upper = %3i %3i \n',v_lower(1),v_upper(1))

pqr = j_lower(:,5); pr = find(pqr ~= 'Q');
pqr = j_lower(:,5); qq = find(pqr == 'Q');
if length(qq) > 0
  figure(1)
  semilogy(freq(pr),stren(pr),'.',freq(qq),stren(qq),'r.')
  fprintf('len0, len1 = %4i %4i \n',length(pr)+length(qq),length(pr))
  accuracy = accuracy(pr,:);
  dipole = dipole(pr);
  elower = elower(pr);
  freq = freq(pr);
  gas_id = gas_id;
  iso = iso(pr);
  j_lower = j_lower(pr,1:9);
  j_upper = j_upper(pr,1:9);
  line_status = line_status(pr);
  p_shift = p_shift(pr);
  reference = reference(pr,1:6);
  stren = stren(pr);
  v_lower = v_lower(pr);
  v_upper = v_upper(pr); 
  w = w(pr);
  w_s = w_s(pr);
  w_temp = w_temp(pr);
  
  clear fnameIN loader pqr pr qq
  fnameOUT = ['hit' num2str(band) '_pr.mat'];
  saver = ['save ' fnameOUT ' accuracy dipole elower freq gas_id iso '];
  saver = [saver      ' j_lower j_upper line_status p_shift reference '];
  saver = [saver      ' stren v_lower v_upper w w_s w_temp'];
  fprintf(1,'%s \n',saver);
  eval(saver);

  whos

  end