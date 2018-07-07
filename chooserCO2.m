function [fhaha,y,fc,theRATIO]=chooserCO2(tempr,bandtype,p,ps,qamt,...
                whichPQR,LVF,IO,birn,nbox,pqrfreq,outwave,theRATIO,homepath)
%this file calls the relevant PQR CO2 branch

%whichPQR is code to tell which bands to add in     01 = Q delt pi  
%                                                   02 = Q sig  pi
 
%                                                  -11 = P sig  sig  
%                                                  -12 = P delt delt  
%                                                  -13 = P pi   pi  
%                                                  -14 = P sig  pi  
%                                                  -15 = P delt pi  
 
%                                                  +11 = R sig  sig  
%                                                  +12 = R delt delt  
%                                                  +13 = R pi   pi  
%                                                  +14 = R sig  pi  
%                                                  +14 = R delt pi  

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
if (whichPQR == 01)  
  eval(['cd ' homepath 'Q_deltpi']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  birnT = birn;        %allow cousin or no chi fcn
  [fhaha,y,fc,theRATIO] = yrun_deltpi(tempr,bandtype,p,ps,qamt,LVF,IOt,...
                   birnT,nbox,pqrfreq,outwave,theRATIO); 
elseif (whichPQR == 02)  
  eval(['cd ' homepath 'Q_sigpi']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  birnT = birn;        %allow cousin or no chi fcn
  [fhaha,y,fc,theRATIO] = yrun_sigpi(tempr,bandtype,p,ps,qamt,LVF,IOt,...
                  birnT,nbox,pqrfreq,outwave,theRATIO); 

elseif (whichPQR == -14)  
  eval(['cd ' homepath 'PR_sigpi']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  birnT = birn;        %allow cousin or no chi fcn
  [fhaha,y,fc,theRATIO] = yrun_sigpi(tempr,bandtype,p,ps,qamt,'P',LVF,IOt,...
        birnT,nbox,pqrfreq,outwave,theRATIO); 
elseif (whichPQR == 14)  
  eval(['cd ' homepath 'PR_sigpi']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  birnT = birn;
  [fhaha,y,fc,theRATIO]=yrun_sigpi(tempr,bandtype,p,ps,qamt,'R',LVF,IOt,...
        birnT,nbox,pqrfreq,outwave,theRATIO); 

elseif (whichPQR == -15)  
  eval(['cd ' homepath 'PR_deltpi']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  birnT = birn;
  [fhaha,y,fc,theRATIO] = yrun_deltpi(tempr,bandtype,p,ps,qamt,'P',LVF,IOt,...
        birnT,nbox,pqrfreq,outwave,theRATIO); 
elseif (whichPQR == 15)  
  eval(['cd ' homepath 'PR_deltpi']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  birnT = birn;
  [fhaha,y,fc,theRATIO] = yrun_deltpi(tempr,bandtype,p,ps,qamt,'R',LVF,IOt,...
        birnT,nbox,pqrfreq,outwave,theRATIO); 
 

%%%%%%%%%%%%%%%%%%  4 um %%%%%%%%%%%%%%%%%%%%%%%%%
%%%      allows line mixing with birnbaum     %%%%
%%%         and then can blend in cousin      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif (whichPQR == -11)  
  eval(['cd ' homepath 'PR_sigsig']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  [fhaha,y,fc,theRATIO] = yrun_sigsig(tempr,bandtype,p,ps,qamt,'P',LVF,IOt,...
             birn,nbox,pqrfreq,outwave,theRATIO); 
elseif (whichPQR == 11)  
  eval(['cd ' homepath 'PR_sigsig']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  [fhaha,y,fc,theRATIO] = yrun_sigsig(tempr,bandtype,p,ps,qamt,'R',LVF,IOt,...
             birn,nbox,pqrfreq,outwave,theRATIO); 
 
elseif (whichPQR == -12)  
  eval(['cd ' homepath 'PR_deltdelt'])
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  [fhaha,y,fc,theRATIO] = ...
         yrun_deltdelt(tempr,bandtype,p,ps,qamt,'P',LVF,IOt,...
                             birn,nbox,pqrfreq,outwave,theRATIO); 
elseif (whichPQR == 12)  
  eval(['cd ' homepath 'PR_deltdelt']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  [fhaha,y,fc,theRATIO] = ...
      yrun_deltdelt(tempr,bandtype,p,ps,qamt,'R',LVF,IOt,...
                             birn,nbox,pqrfreq,outwave,theRATIO); 
 
elseif (whichPQR == -13)  
  eval(['cd ' homepath 'PR_pipi']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  [fhaha,y,fc,theRATIO] = yrun_pipi(tempr,bandtype,p,ps,qamt,'P',LVF,IOt,...
           birn,nbox,pqrfreq,outwave,theRATIO); 

elseif (whichPQR == 13)  
  eval(['cd ' homepath 'PR_pipi']) 
  fprintf(1,'  >> doing linemixing in %s \n',pwd)
  [fhaha,y,fc,theRATIO] = yrun_pipi(tempr,bandtype,p,ps,qamt,'R',LVF,IOt,...
           birn,nbox,pqrfreq,outwave,theRATIO); 
 
end 

if quiet > 0
  fprintf(1,'\n');
end

eval(['cd ' homepath]) 

