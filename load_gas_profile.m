%% parses in profile name
%% could be of form eg IPFILES/gas1.dat or simply IFILES/gas1

iMethod = -1; %% old,   does not want extensions on profname eg IPFILES/std_co2
iMethod = +1; %% new,   allows        extensions on profname eg IPFILES/std_co2.txt
iMethod = +2; %% newer, the -ascii flag allows extensions other than .mat

%a = ['load ' profname];
%eval(a);
str = ['a = load(''' profname ''',''-ascii'');'];
%% fprintf(1,'load profname str = %s \n',str)
eval(str); 

%find last occurence of '/' in profname
ii      = findstr(profname,'/');
ii      = ii(length(ii))+1;
profile = profname(ii:length(profname));


if iMethod == 2
  press       = a(:,2);
  partpress   = a(:,3);
  temperature = a(:,4);
  GasAmt      = a(:,5);
  
elseif iMethod == -1
  %find last occurence of '/' in profname
  ii      = findstr(profname,'/');
  ii      = ii(length(ii))+1;
  profile = profname(ii:length(profname));

  press       = eval([profile '(:,2)']);
  partpress   = eval([profile '(:,3)']);
  temperature = eval([profile '(:,4)']);
  GasAmt      = eval([profile '(:,5)']);

elseif iMethod == +1

  %% fix if there in an extension to the profile, by Tilak
  [ rootPath, VarName, FileExt   ] = fileparts( profname ) ;
  press       = eval( [ VarName '(:,2)' ] ) ;
  partpress   = eval( [ VarName '(:,3)' ] ) ;
  temperature = eval( [ VarName '(:,4)' ] ) ;
  GasAmt      = eval( [ VarName '(:,5)' ] ) ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf
subplot(2,2,1); plot(1:length(GasAmt),GasAmt); 
subplot(2,2,2); plot(1:length(GasAmt),temperature);
subplot(2,2,3); plot(1:length(GasAmt),press);
subplot(2,2,4); plot(1:length(GasAmt),partpress); %pause(1);
title('Gas Profile');

pause(1)
clf
