function [PQRline]=readDAVEhitlin(band,version,strengthM,fname)
%function [PQRline]=readDAVEhitlin(band,version,strengthM,fname)
%is an alternative to hittomat, in that it only pulls out a particular band

global quiet  stretch

gasID=2;

%assume isotop ======= 1 unless stated
isotope=1;

%find all PQR lines for isotope 1 (or 2 or 3) : these are Q_sigpi
if (band==618)
  v_l=2;v_u=3;
elseif (band==648)
  v_l=1;v_u=2;isotope=2;
elseif (band==662)
  v_l=1;v_u=2;isotope=3;
elseif (band == 667)
  v_l=1;  v_u=2; 
elseif (band == 720)
  v_l=2;  v_u=5; 
elseif (band==791)
  v_l=3;v_u=8;
elseif (band==2080)
  v_l=1;v_u=8;

%find all PQR lines for isotope 1 : these are Q_deltpi]
elseif (band==668)
  v_l=2;v_u=4;
elseif (band==740)
  v_l=4;v_u=8;
elseif (band==2093)
  v_l=2;v_u=14;
else
  error('cannot find band!!');
  end

%call hittomat to get the lines
start=band-stretch;
stop=band+stretch;

%line=hittomat2(start,stop,version,strengthM,gasID)
line=hitread(start,stop,strengthM,gasID,fname)

%get the branch info P Q or R
pqr=line.bslq(:,5);

indyes=find((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
    (line.iso == isotope) & ...
    ((pqr' == 'P') | (pqr' == 'Q') | (pqr' == 'R'))  );
ii=length(indyes);
ind=indyes;
%fprintf(1,'\n number of PQR lines selected = %5i \n',ii);

indno=find(~ ((line.ilsgq == v_l) & (line.iusgq == v_u) & ...
    (line.iso == isotope) & ...
    ((pqr' == 'P') | (pqr' == 'Q') | (pqr' == 'R')))  );
ii=length(indno);
%fprintf(1,'\n number of PQR lines not selected = %5i \n',ii);

PQRline.linct  = length(indyes);
PQRline.wnumhx = max(line.wnum(ind));
PQRline.lstat = line.lstat(ind);
PQRline.igas  = line.igas;

PQRline.iso   = line.iso(ind);
PQRline.wnum  = line.wnum(ind);
PQRline.stren = line.stren(ind);
PQRline.tprob = line.tprob(ind);    %transition probablility
PQRline.abroad= line.abroad(ind);  %air broadenend widths
PQRline.sbroad= line.sbroad(ind);  %self broadened widths
PQRline.els   = line.els(ind);
PQRline.abcoef= line.abcoef(ind);  %tempr correct for air broadened widths
PQRline.tsp   = line.tsp(ind);
PQRline.iusgq = line.iusgq(ind);
PQRline.ilsgq = line.ilsgq(ind);
PQRline.uslq  = line.uslq(ind,1:9);
PQRline.bslq  = line.bslq(ind,1:9);
PQRline.ai    = line.ai(ind,1:3);
PQRline.ref   = line.ref(ind,1:6);
PQRline.gasid  = gasID*ones(1,length(indyes));

semilogy(PQRline.wnum,PQRline.stren,'+');
