%this script reads in lines and plots them!!!!!!!

run6do=-1;

usetoth=1;
fmin=400;
fmax=800;

fmin=200;
fmax=3000;

fmin=1400;
fmax=1700;


%%   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];
gid=[1 2 3 4 5 6 7 8 9 10 11 12 15 16 18 19 20 21 22 23 25 26 27];

for pp=1:23
  gasID=gid(pp);

  fnamePRE='/asl/data/hitran/h98.by.gas/g';
  if ((gasID == 1)&(fmin >= 600.103510) & (fmax <= 1650.832720) & ...
       (usetoth == 1))
    fnamePRE='/salsify/users/sergio/KCARTA/SPECTRA/TOTH/newtoth';
    end

  fnamePOST='.dat';
  fnameIN=int2str(gasID);
  fname=[fnamePRE fnameIN fnamePOST];
  [lineORIG]=hitread(fmin,fmax,1e-28,gasID,fname);

  if gasID == 19
    lineORIG=reduce19(lineORIG);
    end

  profname=['load RefProf/refgas' fnameIN];
  eval(profname);
  %find last occurence of '/' in profname
  jj=findstr(profname,'/');jj=jj(length(jj))+1;
  profile=profname(jj:length(profname));

  press=eval([profile '(:,2)']);
  partpress=eval([profile '(:,3)']);
  temperature=eval([profile '(:,4)']);
  GasAmt=eval([profile '(:,5)']);

  maxQ=max(GasAmt);  index=find(GasAmt == maxQ);
  sumQ=sum(GasAmt);

  press=press(index);
  partpress=partpress(index);
  temperature=temperature(index); temperature=296;
  GasAmt=GasAmt(index);

  tempp=[gasID press partpress temperature GasAmt];
  save tempfile tempp -ascii -double

  Q(pp).wnum = lineORIG.wnum;
  Q(pp).q    = lineORIG.stren*maxQ;
  Q(pp).qsum = lineORIG.stren*sumQ;

  if run6do > 0         %runs kind slowly
    figure(1)
    [outfreq,outvect]=run6(gasID,fmin,fmax,0.1,0.1,0.1,1,1,2,25,1,...
                         1e-28,1e-28,'L',-1,'../SPECTRA/tempfile');
  else
    load mass.dat        %get the mass of the isotopes   
    [x,y]=size(mass); 
    %first 32 are the numbers of the isotopes  <iGasID  NumIsotopes  0.0> 
    isotope_num=mass(1:32,2:2); 
    %mext bunch are masses of the isotopes  <iGasID  MassIsotopes Abundance> 
    themass=mass(33:x,2:2); 
    %now pick up the gas masses we need 
    dummy1=sum(isotope_num(1:gasID-1)); 
    dummy2=sum(isotope_num(1:gasID)); 
    dummy=dummy1+1:dummy2; 
    mass_iso=themass(dummy); 
    liso=length(mass_iso); 

    %initialize Q fcns (partition fcns) 
    [A,B,C,D,G] = qtips(gasID,liso); 

    outfreq=fmin:0.5:fmax;
    qfcn=q(A,B,C,D,G,lineORIG,temperature);
    pwr=lineORIG.abcoef;
    for_brd=lineORIG.abroad;
    self_brd=lineORIG.sbroad;
    freq=lineORIG.wnum+press*lineORIG.tsp;  %freq shift
    energy=lineORIG.els;
    s0=lineORIG.stren;
    brd=broad(press,partpress,1.0,for_brd,self_brd,pwr,...
            temperature,gasID);
    strength=find_stren(qfcn,freq,temperature,energy,s0,GasAmt);
    outvect=loop(lineORIG.iso,mass_iso,brd,strength,freq,outfreq,...
                   temperature,lineORIG.linct,length(outfreq),-1);
    end
pause

%  str=['gasID = ' num2str(gasID) ' strength * maxQ'];
%  subplot(211); plot(Q(pp).wnum,Q(pp).q);
%  str=['gasID = ' num2str(gasID) ' lorentz profile',12];
%  subplot(212); plot(outfreq,outvect); title(str,'FontSize',12); pause(0.2)

  figure(2);
  str=['gasID = ' num2str(gasID)];
  if (mod(pp,4) == 1)
    h1=subplot(411); plot(outfreq,outvect); set(h1,'FontSize',12);
    axis([fmin fmax 0 max(outvect)]); title(str,'FontSize',12); pause(0.1)
  elseif (mod(pp,4) == 2)
    h2=subplot(412); plot(outfreq,outvect); set(h2,'FontSize',12);
    axis([fmin fmax 0 max(outvect)]); title(str,'FontSize',12); pause(0.1)
  elseif (mod(pp,4) == 3)
    h3=subplot(413); plot(outfreq,outvect); set(h3,'FontSize',12);
    axis([fmin fmax 0 max(outvect)]); title(str,'FontSize',12); pause(0.1)
  elseif (mod(pp,4) == 0)
    h4=subplot(414); plot(outfreq,outvect); set(h4,'FontSize',12);
    axis([fmin fmax 0 max(outvect)]); title(str,'FontSize',12); pause(0.1)

    adjust41(h1,h2,h3,h4,'even');
    print -dps -append spectralatm.ps
    end

  maxq=max(Q(pp).q);
  maxqsum=max(Q(pp).qsum);

  fprintf(1,'gid= %3i : qmax= %8.5e stren*qmx= %8.5e \n',gasID,maxQ,maxq);

  max1(pp)=maxq;
  max2(pp)=maxqsum;

  end

clf
semilogy(gid,max1,'+',gid,max2,'r+'); grid
print -dps -append spectralatm.ps

clear str p* r* GasAmt f* gasID pp jj maxQ temperature usetoth lineORIG
clear maxq maxqsum sumQ

%now sort

%   refgas1
%   refgas2
%   refgas3
%   refgas4
%   refgas5
%   refgas6
%   refgas7
%   refgas8
%   refgas9
%   refgas10
%   refgas11
%   refgas12
%   refgas15
%   refgas16
%   refgas18
%   refgas19
%   refgas20
%   refgas21
%   refgas22
%   refgas23
%   refgas25
%   refgas26
%   refgas27
%   refgas51
%   refgas52
%   refgas53
%   refgas54
%   refgas55
%   refgas56
%   refgas57
%   refgas58
%   refgas59
%   refgas60
%   refgas61
%   refgas62
%   refgas63
