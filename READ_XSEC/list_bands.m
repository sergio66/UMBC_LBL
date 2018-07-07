function [bands,xs] = list_bands(gid,HITRAN)

if nargin == 1
  HITRAN = 2016;
end

xs = read_xsec(gid,[],HITRAN);

[nrec, nband] = size(xs);
  
for b = 1:nband
  
  v1 = max([xs(:,b).v1]);
  v2 = min([xs(:,b).v2]);
  
  fprintf(1, '%4d %4d %10.3f %10.3f\n', gid, b, v1, v2);
  bands.v1(b) = v1;
  bands.v2(b) = v2;  

end % band loop


