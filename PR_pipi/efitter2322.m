function [elower,jall] = efitter2322(jr,band,elowerr,elower,prb);  

% Creating the initial parameters:
%   jr are the hittomat j levels.
%   elowerr are the corresponding lower state energy levels.  Note that the 
%      lower state energies given on the tape are TOTAL energies, not just
%      rotational energies.
%   bds are the starting values for the fit.

%have to worry about ee ff duplication

global quiet v12p1

%disp('fitting for missing energy levels')
bds = [670; 	% Evib  for even levels 
       0.39;	% B
       1.5e-7];	% D

%%%%%%%%%%%%%%%%%%%% EVEN
index = 1:2:length(jr); %all e's or all f's 
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
if pflag > 0  
        clf;hold off;  
        plot(jrR,y)  
        xlabel('f');ylabel('Final Diff');pause(0.1)
        pause
        end 
ebdeven = ebd;

%%%%%%%%%%%%%%%%%%%% ODD 

bds = [670; 	% Evib  for odd levels 
     0.39;	% B
     1.5e-7];	% D

index = 2:2:length(jr); %all e's or all f's 

jrR = jr(index);
elowerrR = elowerr(index);

if v12p1 > 0
  options(1) = 0; 
  options(5) = 0; options(7) = 1; options(14) = 200;
  ebd = leastsq('efit',bds,options,'efitgrad',jrR,elowerrR);
else
  %ebd = lsqnonlin(@efit,bds,[],[],options,jrR,elowerrR);
  ebd  = bds;   
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
        xlabel('f');ylabel('Final Diff');pause(0.1)
        pause
        end 
ebdodd = ebd;

%%%%%%%%%%%%% NOW PUT IN ENERGYIES FOR MISSING LEVELS

% Coefficients which minimize the difference between the functional form of the
% energy and the given values are returned in the variable ebd.
% Details of this step are shown in efit.m

% The energy levels for ALL J are calculated:

jallall = [1:max(jr)]; 
jall = zeros(1,2*length(jallall));
jall(1:2:length(jall)) = jallall;
jall(2:2:length(jall)) = jallall;
jall = jall';

%do es and fs separately

%%%%%%%%%%%%%%%%%%%% ODD 

index = 1:2:length(jall); %odd j's 
ebd = ebdodd;
elower(index) = (ebd(2)*jall(index).*(jall(index)+1) - ...
                ebd(3)*(jall(index).*(jall(index)+1)).^2);
%now all values which came from tape are re-inserted
index = 1:2:length(jr); %odd j's 
if ((prb == 'p')|(prb=='P'))
  index_r = index+1;
else
  index_r = index;
  end

elower(index_r) = elowerr(index)-ebd(1);

%%%%%%%%%%%%%%%%%%%% EVEN

%do odds and evens separately
index = 2:2:length(jall); %even j's 
ebd = ebdeven;
elower(index) =  (ebd(2)*jall(index).*(jall(index)+1) - ...
                ebd(3)*(jall(index).*(jall(index)+1)).^2);
index = 1:2:length(jr); %even j's 
if ((prb == 'p')|(prb == 'P'))
  index_r = index+1;
else
  index_r = index;
  end

elower(index_r) = elowerr(index)-ebd(1);

%plot(jr,elowerr,'+',jall,elower); pause(0.1)

%  elower will be used as input for calculating W
%  NoTe that elower is only the ROTATIONAL energy.
% elowerrf =  (ebd(2)*jr.*(jr+1) - ebd(3)*(jr.*(jr+1)).^2);

clear start2 bds
%save efittervar.mat elowerr elower jr
