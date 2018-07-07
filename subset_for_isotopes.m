function [lineNEW,which_isotope] = subset_for_isotopes(lineORIG,which_isotope0)

%% suppose isotope list = 1 2 3 4 5 6
%% a) if which_isotope == 0   use all isotopes  1 2 3 4 5 6
%% b) if unique(lineORIG.iso) = 1 2 3 4 5 6
%%       which_isotope = [-1 5 6]
%%    this means keep isotopes 1 2 3 4, throw out 5 6
%% c) if unique(lineORIG.iso) = 1 2 3 4 5 6
%%       which_isotope = [1 2 3 4 ]
%% this means keep isotopes 1 2 3 4, throw out 5 6

%%disp('qtips04.m : for water the isotopes are 161  181  171  162  182  172');

%% see qtips04.m
fprintf(1,'there are %2i unique isotopes for gas %2i \n',length(unique(lineORIG.iso)),lineORIG.igas)
unique(lineORIG.iso)
disp(' ')

lineNEW       = lineORIG;
which_isotope = which_isotope0;

if which_isotope ~= 0
  fprintf(1,'  user has specified isotope list : \n')
  which_isotope
  fprintf(1,'  while this is the list of HITRAN isotopes \n')
  unique(lineORIG.iso)
  figure(1); 
  scatter(lineORIG.wnum,log10(lineORIG.stren),10,lineORIG.iso); colorbar
  title('colorbar');

  if length(intersect(which_isotope,-1)) == 1
    which_isotope = sort(which_isotope);
    if which_isotope(1) ~= -1
      disp(' this option expects which_isotope = [-1 iso1 iso2 .. isoN]');
      disp('    where iso1 iso2 .. isoN = list of isotopes to EXCLUDE');
      error('try again!')
    end
    %%this means keep all isotopes except those which are in the list
    %% eg if unique(lineORIG.iso) = 1 2 3 4 5 6
    %%       which_isotope = [-1 5 6]
    %% this means keep isotopes 1 2 3 4, throw out 5 6
    dada = 1 : length(lineORIG.iso);
    baba = [];
    %% this means user wants to REMOVE this isotope(s)
    for jj = 2 : length(which_isotope)
      woo = find(lineORIG.iso == which_isotope(jj));
      if length(woo) > 0
        baba = [baba woo];
      end
    end
    zaza = setdiff(dada,baba);
    %[length(dada) length(baba) length(zaza)]
    figure(2);
    plot(lineORIG.wnum,log10(lineORIG.stren),'b.',...
        lineORIG.wnum(baba),log10(lineORIG.stren(baba)),'ro')
    title('blue = all lines   red = lines to be removed');
    optvarx = fieldnames(lineORIG);
    for i = 1 : length(optvarx)
      str = ['daname = lineORIG.' optvarx{i} ';']; eval(str);
      [mmm,nnn] = size(daname);      
      if length(daname) > 1
        if nnn == 1 & mmm > 1
          str = ['lineNEW.' optvarx{i} ' = lineORIG.' optvarx{i} '(zaza);'];
	elseif mmm == 1 & nnn > 1
          str = ['lineNEW.' optvarx{i} ' = lineORIG.' optvarx{i} '(zaza);'];
        elseif mmm > 1 & nnn == lineORIG.linct
          str = ['lineNEW.' optvarx{i} ' = lineORIG.' optvarx{i} '(:,zaza);'];
        elseif nnn > 1 & mmm == lineORIG.linct
          str = ['lineNEW.' optvarx{i} ' = lineORIG.' optvarx{i} '(zaza,:);'];
	end
        eval(str);
      else
        str = ['lineNEW.' optvarx{i} ' = lineORIG.' optvarx{i} ';'];
        eval(str)
      end
    end
    lineNEW.linct = length(lineNEW.iso);
    disp('In subset_for_isotopes.m')
    fprintf(1,'  all    isotopes have total %6i lines \n',lineORIG.linct)
    fprintf(1,'  chosen isotopes have total %6i lines \n',lineNEW.linct)

    %% now do the "removal" of isotopes
    zaza = setdiff(unique(lineORIG.iso),which_isotope(2:length(which_isotope)));
    which_isotope = zaza;

  elseif length(intersect(which_isotope,-1)) == 0
    which_isotope = sort(which_isotope);
    if min(which_isotope) < 1
      disp(' this option expects which_isotope = [iso1 iso2 .. isoN]');
      disp('    where iso1 iso2 .. isoN = list of isotopes to INCLUDE');
      error('try again!')
    end
    %%this means keep only isotopes in the list
    %% eg if unique(lineORIG.iso) = 1 2 3 4 5 6
    %%       which_isotope = [1 2 3 4 ]
    %% this means keep isotopes 1 2 3 4, throw out 5 6
    dada = 1 : length(lineORIG.iso);
    baba = [];
    %% this means user wants to KEEP these isotope(s)
    for jj = 1 : length(which_isotope)
      woo = find(lineORIG.iso == which_isotope(jj));
      if length(woo) > 0
        baba = [baba woo];
      end
    end
    zaza = baba;
    [length(dada) length(baba) length(zaza)]
    figure(2);
    plot(lineORIG.wnum,log10(lineORIG.stren),'b.',...
        lineORIG.wnum(baba),log10(lineORIG.stren(baba)),'ro')
    title('blue = all lines   red = lines to be removed');
    optvarx = fieldnames(lineORIG);
    for i = 1 : length(optvarx)
      str = ['daname = lineORIG.' optvarx{i} ';']; eval(str);
      [mmm,nnn] = size(daname);
      if length(daname) > 1
        if nnn == 1 & mmm > 1
          str = ['lineNEW.' optvarx{i} ' = lineORIG.' optvarx{i} '(zaza);'];
	elseif mmm == 1 & nnn > 1
          str = ['lineNEW.' optvarx{i} ' = lineORIG.' optvarx{i} '(zaza);'];
        elseif mmm > 1 & nnn == lineORIG.linct
          str = ['lineNEW.' optvarx{i} ' = lineORIG.' optvarx{i} '(:,zaza);'];
        elseif nnn > 1 & mmm == lineORIG.linct
          str = ['lineNEW.' optvarx{i} ' = lineORIG.' optvarx{i} '(zaza,:);'];
	end
	%fprintf(1,'%s \n',str)
        eval(str);
      else
        str = ['lineNEW.' optvarx{i} ' = lineORIG.' optvarx{i} ';'];
        eval(str)
      end
    end
    %% which_isotope remains unchanged
    lineNEW.linct = length(lineNEW.iso);

    disp('In subset_for_isotopes.m')
    fprintf(1,'  all    isotopes have total %6i lines \n',lineORIG.linct)
    fprintf(1,'  chosen isotopes have total %6i lines \n',lineNEW.linct)

  else
    error('cannot understand your which_isotope list')
  end
end
