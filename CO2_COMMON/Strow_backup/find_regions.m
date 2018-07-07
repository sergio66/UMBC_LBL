global quiet p2311_21_jmax p2350_jmax r2350_jmax pr2351_jmax

matr=getmatr(band,prb,freqr);    %%%this is really dummy, as everything
                                 %%%important is reset anyways!!!
if (prb == 'p')
  prb='P';
elseif (prb == 'r')
  prb='R';
  end

orig=-1;      
if (orig > 0) 
  %%%% this is the code prior to March 1999 where we either had line mixing 
  %%%% w/wo birnbaum, or cousin
  if ((INratio <= 0.05) & (INratio >= 0.005))
    %things could get dicey so change bounds
    matr(1,1)=matr(1,1)-40;
    matr(2,1)=matr(2,1)-40;
    matr(3,1)=matr(3,1)+40;
    matr(4,1)=matr(4,1)+40;
    end
else 
  %%%this is the new code, where we turn over to Cousin after 15 wavenumbers
  %%% and close to bands, we have line mixing

  x=band;

  %%%after July 2002
  jmax = 60;
  if (x==2351)
    [gg,hh]=sort(jr); im=find(gg<=jmax); im=im(length(im));
    %%%only do low jq's for mixing upto jq=80
    matr(1,1)=freqr(hh(im))-1.0;              %%%%%%irrelevant
    matr(2,1)=freqr(hh(im));
    matr(3,1)=max(freqr)+0.025;
    matr(4,1)=max(freqr)+0.05;                %%%%%%irrelevant
    end

  %%%orig, before July 2002
  if (x==2351) 
    matr(1,1)=min(freqr)-20;                %%%%%%irrelevant
    matr(2,1)=min(freqr)-15;
    matr(3,1)=max(freqr)+15;
    matr(4,1)=max(freqr)+20;                %%%%%%irrelevant
    end

  if (x==2350) 
    if (prb=='P')  
      %%%%%%%before aug 5, 2002
      %%%%%%%has 40 cm-1 blending regions!!!!!!! on right
      %%%%%%%has 80 cm-1 blending regions!!!!!!! on left  for P branch
      matr(1,1)=min(freqr)-85;                %%%%%%irrelevant
      matr(2,1)=min(freqr)-80;
      matr(3,1)=max(freqr)+40;
      matr(4,1)=max(freqr)+45;                %%%%%%irrelevant

    elseif (prb=='R')  
      %%%%%%%has 40 cm-1 blending regions!!!!!!! on right
      %%%%%%%has 10 cm-1 blending regions!!!!!!! on left  for R branch
      %% gives bugs matr(1,1)=min(freqr)-15;                %%%%%%irrelevant
      %% gives bugs matr(2,1)=min(freqr)-10;
      matr(1,1)=min(freqr)-26;                %%%%%%irrelevant
      matr(2,1)=min(freqr)-21;
      matr(3,1)=max(freqr)+40;
      matr(4,1)=max(freqr)+45;                %%%%%%irrelevant
      end
    end

  %%% originally : (x==2310)|(x == 2320) were considered below!!! jmax === 50
  %%% after      : July 2002, let jmax = 90
  if ((x==2310)|(x == 2320))
    if (prb=='P')
      jmax = p2350_jmax;
      [gg,hh]=sort(jr); im=find(gg<=jmax); im=im(length(im));
      %%%only do low jq's for mixing upto jq=80
      matr(1,1)=freqr(hh(im))-1.0;              %%%%%%irrelevant
      matr(2,1)=freqr(hh(im));
      matr(3,1)=max(freqr)+0.025;
      matr(4,1)=max(freqr)+0.05;                %%%%%%irrelevant
    elseif (prb=='R')
      jmax = 90;
      [gg,hh]=sort(jr); im=find(gg<=jmax); im=im(length(im));
      matr(1,1)=min(freqr)-0.05;                %%%%%%irrelevant
      matr(2,1)=min(freqr)-0.025;
      matr(3,1)=freqr(hh(im));
      matr(4,1)=freqr(hh(im))+1.0;              %%%%%%irrelevant
      end
    end

  %%% originally upto July 2002
  %%% originally : if ((x==2310)|(x == 2320)|(x==2352)|(x==2353)|(x==2354))
  %%% jmax ======== 50 for these bands
  jmax = 50;
  if ((x==2352)|(x==2353)|(x==2354))
    if (prb=='P')
      [gg,hh]=sort(jr); im=find(gg<=jmax); im=im(length(im));
      %%%only do low jq's for mixing upto jq=50
      matr(1,1)=freqr(hh(im))-1.0;              %%%%%%irrelevant
      matr(2,1)=freqr(hh(im));
      matr(3,1)=max(freqr)+0.025;
      matr(4,1)=max(freqr)+0.05;                %%%%%%irrelevant
    elseif (prb=='R')
      [gg,hh]=sort(jr); im=find(gg<=jmax); im=im(length(im));
      matr(1,1)=min(freqr)-0.05;                %%%%%%irrelevant
      matr(2,1)=min(freqr)-0.025;
      matr(3,1)=freqr(hh(im));
      matr(4,1)=freqr(hh(im))+1.0;              %%%%%%irrelevant
      end
    end

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now check that matr(1,1) < matr(2,1)
if matr(1,1) > matr(2,1)
  tempppp=matr(2,1);
  matr(2,1)=matr(1,1);
  matr(1,1)=tempppp;
  end

%now check that matr(3,1) < matr(4,1)
if matr(4,1) > matr(4,1)
  tempppp=matr(4,1);
  matr(4,1)=matr(3,1);
  matr(3,1)=tempppp;
  end

fprintf(1,'blending for %s %4i = %8.6f %8.6f \n',prb,band,matr(2,1),matr(3,1))