%initialize Q fcns (partition fcns) if possible

if length(intersect(hitran_version,hlist_qtips)) == 1
  disp('old polynomial qtips will work ...');  
  [A,B,C,D,G] = qtips(gasID,liso);
else
  disp('need new version of qtips ...');
  A = zeros(liso,1);    %for polynom
  B = zeros(liso,1);    %for polynom
  C = zeros(liso,1);    %for polynom
  D = zeros(liso,1);    %for polynom
  G = zeros(liso,1);    %nuclear degeneracy factor
  end