function [elower,jall] = efitter(jr,band,elowerr,elower,prb);  

% Creating the initial parameters:
%   jr are the hittomat j levels.
%   elowerr are the corresponding lower state energy levels.  Note that the 
%      lower state energies given on the tape are TOTAL energies, not just
%      rotational energies.
%   bds are the starting values for the fit.

global quiet v12p1

%disp('fitting for missing energy levels')
bds = [670; 	% Evib  for even levels 
     0.39;	% B
     1.5e-7];	% D

bds = [630;       % Evib 
     0.39;              % B 
     1.5e-7];   % D 
 
if ((prb == 'R')|(prb == 'r')) 
  bds(1) = 1000; 
  end 

if (band == 2321)
  bds(1) = 200;
  end
 

%%%%%%%%%%%%%%%%%%%% EVEN
index = find(mod(jr,2) == 0); %even j's 
jrR = jr(index);
elowerrR = elowerr(index);

if v12p1 > 0
  clear options
  options(1) = 0; 
  options(5) = 0; options(7) = 1; options(14) = 200;
  ebd  =  leastsq('efit',bds,options,'efitgrad',jrR,elowerrR);
else
  clear options; options  =  optimset('MaxIter',1000);
  %ebd  =  lsqnonlin(@efit,bds,[],[],options,jrR,elowerrR);
  ebd  =  bds;   
  params  =  [jrR elowerrR];   
  ebd  =  fsolve(@(ebd) efit1(ebd,params),bds,options);   
  end
evib = ebd(1);      
bs = ebd(2);        
ds = ebd(3);        
erot = (bs*jrR.*(jrR+1) -ds*(jrR.*(jrR+1)).^2);  
y = erot+evib-elowerrR;  
  
pflag = quiet;  
if pflag>0  
        clf;hold off;  
        plot(jrR,y)  
        xlabel('f');ylabel('Final Diff');pause(0.01)
        end 
ebdeven = ebd;

%%%%%%%%%%%%%%%%%%%% ODD 

bds = [670; 	% Evib  for odd levels 
     0.39;	% B
     1.5e-7];	% D

if (band == 2321)
  bds(1) = 200;
  end

index = find(mod(jr,2) == 1); %odd j's 
jrR = jr(index);
elowerrR = elowerr(index);

if v12p1 > 0
  options(1) = 0; 
  options(5) = 0; options(7) = 1; options(14) = 200;
  ebd  =  leastsq('efit',bds,options,'efitgrad',jrR,elowerrR);
else
  clear options; options = optimset('MaxIter',1000);
  %ebd = lsqnonlin(@efit,bds,[],[],options,jrR,elowerrR);
  ebd = bds;   
  params = [jrR elowerrR];   
  ebd = fsolve(@(ebd) efit1(ebd,params),bds,options);   
  end
evib = ebd(1);      
bs = ebd(2);        
ds = ebd(3);        
erot = (bs*jrR.*(jrR+1) -ds*(jrR.*(jrR+1)).^2);  
y = erot+evib-elowerrR;  
  
pflag = quiet;  
if pflag > 0  
        clf;hold off;  
        plot(jrR,y)  
        xlabel('f');ylabel('Final Diff');pause(0.01)
        end 
ebdodd = ebd;

%%%%%%%%%%%%% NOW PUT IN ENERGYIES FOR MISSING LEVELS

% Coefficients which minimize the difference between the functional form of the
% energy and the given values are returned in the variable ebd.
% Details of this step are shown in efit.m

% The energy levels for ALL J are calculated:

jall = [1:max(jr)]'; 

%do odds and evens separately

%%%%%%%%%%%%%%%%%%%% ODD 

index = find(mod(jall,2) == 1); %odd j's 
ebd = ebdodd;
elower(index) =  (ebd(2)*jall(index).*(jall(index)+1) - ...
                ebd(3)*(jall(index).*(jall(index)+1)).^2);
%now all values which came from tape are re-inserted
index = find(mod(jr,2) == 1); %odd j's 
if ((prb  ==  'p')|(prb == 'P'))
  index_r = index+1;
else
  index_r = index;
  end

elower(index_r) = elowerr(index)-ebd(1);

%%%%%%%%%%%%%%%%%%%% EVEN

%do odds and evens separately
index = find(mod(jall,2)  ==  0); %even j's 
ebd = ebdeven;
elower(index) =  (ebd(2)*jall(index).*(jall(index)+1) - ...
                ebd(3)*(jall(index).*(jall(index)+1)).^2);
index = find(mod(jr,2)  ==  0); %even j's 
if ((prb == 'p')|(prb == 'P'))
  index_r = index+1;
else
  index_r = index;
  end

elower(index_r) = elowerr(index)-ebd(1);

%plot(jr,elowerr,'+',jall,elower); pause(0.01)

%  elower will be used as input for calculating W
%  NoTe that elower is only the ROTATIONAL energy.
% elowerrf = (ebd(2)*jr.*(jr+1) - ebd(3)*(jr.*(jr+1)).^2);

clear start2 bds
%save efittervar.mat elowerr elower jr
