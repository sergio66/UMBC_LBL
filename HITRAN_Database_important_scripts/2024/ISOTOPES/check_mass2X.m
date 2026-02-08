ngas = 61;

x = load('molparam.txt');

fid = fopen('molparam.txt','r');
% Read all lines & collect in cell array
txt = textscan(fid,'%s','delimiter','\n'); 
fclose(fid);

moo = find(x(:,3) == 0 & x(:,4) == 0 & x(:,5) == 0);
fprintf(1,'expected %2i molecules, found %2i places where columns 3,4,5 were 0 \n',ngas,length(moo));
sum_iso = sum(x(moo,2));

exp_count = 1;                      %% first line is  % Molecule # Iso Abundance     Q(296K)      gj    Molar Mass(g)
exp_count = exp_count + 2*ngas;  %% these would be for eah molecule      eg %   H2O (1) and 1  7  0 0 0
exp_count = exp_count + sum_iso;    %% these would be the actual iso info

[mm,nn] = size(x);
mm = cellfun(@numel, txt);
fprintf(1,'based on %2i molecules in file, we expect %3i lines in files  \n',ngas,exp_count);
fprintf(1,'and we read in %3i lines \n',mm);

if mm == exp_count
  disp('YAY  length(txt) == exp_count')
  iYorN = +1;  
else
  disp('BOO  length(txt) != exp_count')
  iYorN = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iYorN < 0
  disp(' ')
  disp('found a problem : so here is relevant info')
  disp(' ')
  
  %% find how well info joves and jives, or not
  moo = find(x(:,3) == 0 & x(:,4) == 0 & x(:,5) == 0);
  disp('quickly checking file --- if code crashes, you can figure out which gasID has a problem')
  ii = 1;    % Molecule # Iso Abundance     Q(296K)      gj    Molar Mass(g)
  for iMol = 1 : ngas
    fprintf(1,'gasID = %2i \n',iMol);
    iWoo = moo(iMol);
    
    ii = ii + 1;  fprintf(1,'%s \n',txt{1}{ii})
    ii = ii + 1;  fprintf(1,'%s \n',txt{1}{ii})
    strx = txt{1}{ii};
    ii5  = sscanf(strx,'%i');  %% gasiD num_iso(gasID) 0 0 0
    fprintf(1,'number of isotopes = %2i from %s is %2i \n',x(iWoo,2),txt{1}{ii},ii5(2))
    if x(iWoo,2) ~= ii5(2)
      fprintf(1,'gasID = %2i ... check previous gas \n',iMol)
      error('there is a problem, probably for gasID provious to this')
    end
    ii = ii + x(iWoo,2);
    disp('quickly checking file --- if code crashes, you can figure out previous gasID has a problem')  
    disp(' ')
  end
end
