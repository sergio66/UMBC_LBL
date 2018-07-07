function [lor,lormix]=klormix(freqq,f,ymix,jq,temperature,w_forq,w_selfq,...
                       strenqt,stuff);

%%%%%%%%%%%%%%%%%%%%%%%%%%% klormix.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	This file computes:
%		(1) Lorentz with 1st order mixing absorption coefficient
%		(2) Lorentz absorption coefficient
%	It also checks SUM{Si*Yi} to see if detailed balance worked.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% K_scales are used to convert absorption coefficients into the corrrect units.
% This was taken from Dr. Strow's programs and I haven't really taken the time
% yet to go through the actual derivation.

density=stuff.density; 
pressure_self=stuff.pressure_self; 
pressure_for=stuff.pressure_for; 
pressure_ref=stuff.pressure_ref; 
temperature_ref=stuff.temperature_ref; 
path_length=stuff.path_length; 
bsm=stuff.bsm; 
 
K_scale_mixing=density*pressure_self/pressure_ref*temperature_ref*... 
        path_length/temperature/pi; 
no_lines=length(freqq); no_pts=length(f); 
k=zeros(no_pts,1); 
 
frequency_shift=0; 

%%%%%%%%%%%%%%%%%% Compute Klor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('computing lorentz absorption coefficient')
lor=zeros(length(f),1);
for i=1:length(freqq)
   w_tot=(pressure_self*w_selfq(i)+pressure_for* ...
	w_forq(i))/pressure_ref;
   lor=lor+strenqt(i)*(w_tot)*(f/freqq(i))./ ...
	((f-freqq(i)).^2 +(w_tot).^2); 
end

lor=lor*K_scale_mixing*bsm;


%%%%%%%%%%%%%%%%%% Use Rosenkranz' approximation for Kmix %%%%%%%%%%%
disp('computing Lorentz with 1st order mixing absorption coefficient')
lormix=zeros(length(f),1);  
for i=1:length(freqq)
   w_tot=(pressure_self*w_selfq(i)+pressure_for* ...
	w_forq(i))/pressure_ref;
   lormix=lormix+strenqt(i)*(f/freqq(i)).*(w_tot+...
   (f-freqq(i))*ymix(i))./((f-freqq(i)).^2+(w_tot).^2);
end
lormix=lormix*K_scale_mixing*bsm;

%%%%%%%%%%%% Calculate SUM{Si*Yi} to check detailed balance %%%%%%%%%%% 
%  See Dr. Strow's SPIE paper for an explanation.
%sum=0;
%for i=1:length(freqq) 
%	sum=sum+strenqt(i)*ymix(i);
%end
%sum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




