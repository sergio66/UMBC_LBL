function [bands,xs] = list_bands(gid,HITRANyear)

if nargin == 1
  HITRANyear = 2016;
  HITRANyear = 2020;
  HITRANyear = 2024;    
end

xs = read_xsec(gid,[],HITRANyear);

[nrec, nband] = size(xs);
if nband > 0
  for b = 1:nband
  
    v1 = max([xs(:,b).v1]);
    v2 = min([xs(:,b).v2]);
  
    fprintf(1, '%4d %4d %10.3f %10.3f\n', gid, b, v1, v2);
    bands.v1(b) = v1;
    bands.v2(b) = v2;  

  end % band loop
else
  bands.v1 = 0;
  bands.v2 = 0;
end
