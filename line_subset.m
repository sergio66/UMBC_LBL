function x = line_subset(line,index)

%% assumes you have read in bunch of lines using hitread, and now want to 
%% subset one of them eg 
%%  filename12 = '/asl/data/hitran/h12.by.gas/g3.dat';
%%  line12 = hitread(1042.5,1043.5,0,3,filename12,-1);
%%  woo12 = find(abs(line12.wnum-1042.91276) < 0.001);
%%  line_subset(line12,woo12)

x.igas  = line.igas;
x.linct = 1;

fn = fieldnames(line);
for ii = 1 : length(fn)
  linefield = fn{ii};
  if ~strcmp(linefield,'igas') & ~strcmp(linefield,'linct')
    str = ['boo = line.' linefield ';'];
    eval(str);
    [mm,nn] = size(boo);
    if mm == 1 & nn == line.linct
      str = ['x.' linefield ' = line.' linefield '(' num2str(index) ');'];
    elseif mm > 1 & nn == line.linct
      str = ['x.' linefield ' = line.' linefield '(:,' num2str(index) ');'];
    elseif mm == line.linct & nn > 1
      str = ['x.' linefield ' = line.' linefield '(' num2str(index) ',:);'];
    end
    eval(str)
  end
end