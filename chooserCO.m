function [fhaha,y,fc,theRATIO]=chooserCO(tempr,bandtype,p,ps,qamt,...
                whichPQR,LVF,IO,birn,nbox,pqrfreq,outwave,theRATIO,homepath)
%this file calls the relevant PQR CO2 branch

%whichPQR is code to tell which bands to add in 
%                                                  -11 = P sig  pi  
%                                                  +11 = R sig  pi

%LVF  is for lorentz, voigt, full line mixing
%IO   is for no, first order      line mixing
%birn is for no, yes birnbaum (always turned off for Q_sigpi,Q_deltpi)

global quiet

clear loader efitter y1s y1ser trans_pop klormix full_mix4 
clear efit efitter orderer wfun1co2 wfunco2 wgradco2 wfunco2er 

IOt=IO;
%%%%%%%%%% these next 4 lines are commented out on March 14, 2002 %%%%%%%%%%
%if ((birn=='c') | (birn == 'C'))
%  %cousin includes line mixing!!!!!!!!!!, about 10 cm-1 away from line center
%  IOt='0';
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% 15 um %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      allows line mixing with birnbaum     %%%%
%%%         and then can blend in cousin      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (whichPQR == -11)  
  eval(['cd ' homepath 'CO_PR_sigpi']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  [fhaha,y,fc,theRATIO] = yrun_sigsig_co(tempr,bandtype,p,ps,qamt,'P',LVF,IOt,...
             birn,nbox,pqrfreq,outwave,theRATIO); 
elseif (whichPQR == 11)  
  eval(['cd ' homepath ' CO_PR_sigpi']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  [fhaha,y,fc,theRATIO] = yrun_sigsig_co(tempr,bandtype,p,ps,qamt,'R',LVF,IOt,...
             birn,nbox,pqrfreq,outwave,theRATIO);  
end 

if quiet > 0
  fprintf(1,'\n');
end

eval(['cd ' homepath]) 

