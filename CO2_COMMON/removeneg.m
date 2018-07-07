function [voinew]=removeneg(outwave,voi);
%this function makes sures all k's >= 0

global quiet 

voinew=voi;

b=sum(voinew<0); %if any of the terms are negative, have to fix
if (b > 0)
  fprintf(1,'some of the k''s are negative .. resetting to 0 \n');
  clear b
  b=find(voinew < 0);      %find indices where k is negative
  bd=diff(b);              %find if indices are all adjacent or spaced out

  if (sum(bd==ones(size(bd))) == length(bd)) 
    %%%%all the places where the k's are negative, are adjacent
    %%%% eg k=[5 2 1 -0.1 -0.5 -1 -0.5 -0.25 -0.00125 0.025 1 5 10 20 30 40];
    %%%% eg k=[-0.1 -0.5 -1 -0.5 -0.25 -0.00125 0.025 1 5 10 20 30 40];
    %%%% eg k=[5 2 1 -0.1 -0.5 -1 -0.5 -0.25 -0.00125 ];
    if ((b(1) > 1)  & (b(length(b)) < length(voinew)))
      %%%% eg k=[5 2 1 -0.1 -0.5 -1 -0.5 0.5 1.5 10.5];
      bm=b(1)-1;
      bp=b(length(b))+1;
      %%%slope=voinew(bp)-voinew(bm); slope=slope/(bp-bm); old code
      %%%voinew(b)=slope*(b-bm)+voinew(bm);                old code
      bee=[bm:bp];
      slope=voinew(bp)-voinew(bm); slope=slope/(outwave(bp)-outwave(bm));
      voinew(bee)=slope*(outwave(bee)-outwave(bm))+voinew(bm);
    elseif ((b(1) == 1)  & (b(length(b)) < length(voinew)))
      %%%% eg k=[ -0.1 -0.5 -1 -0.5 0.5 1.5 10.5];
      bm=b(1);
      bp=b(length(b))+1;
      bee=[bm:bp];
      slope=voinew(bp)-voinew(bm); slope=slope/(outwave(bp)-outwave(bm));
      voinew(bee)=slope*(outwave(bee)-outwave(bm))+voinew(bm);
    elseif ((b(1) > 1)  & (b(length(b)) == length(voinew)))
      %%%% eg k=[5 2 1 -0.1 -0.5 -1 -0.5];
      bm=b(1)-1;
      bp=b(length(b));
      bee=[bm:bp];
      slope=voinew(bp)-voinew(bm); slope=slope/(outwave(bp)-outwave(bm));
      voinew(bee)=slope*(outwave(bee)-outwave(bm))+voinew(bm);
      end

    %%%if things are still messed, force naughty k's to 0.0
    b=find(voinew < 0);
    voinew(b) = 0.0;

    if quiet > 0
      plot(outwave,exp(-voinew),outwave(b),exp(-voinew(b)),'r+')
      title('after getting rid of negative ks');
      end

  else 
    %%%%the places where the k's are negative are not adjacent
    %break b into chunks ..........     

    %%%%%%break it into chunks where the k's are negative
    %% x=[5 2 1 -0.1 -0.5 -1 -0.5 0.5 1.0 2.0 -1.0 -2.0 -1.0 2.0 3.0 -1 -2];
    %% b= 4     5     6     7    11    12    13    16    17    %bad indices
    %% bd=1     1     1     4     1     1     3     1          %finds adjacent

    index=(bd==ones(size(bd)));
    rara=find((~index==1));
    %%rara = 4 7         %%tells us in bd, evrthing adjacent except 4,7
    %%thus indices 1-4 ok, 5-7 ok, 8-9ok

    ii=1;
    raratop(ii)=b(ii);
    rarabot(ii)=b(rara(ii));
    for ii=2:length(rara)       %there are length(rara)+1 bad regions
      raratop(ii)=b(rara(ii-1)+1);
      rarabot(ii)=b(rara(ii));
      end
    ii=length(rara)+1;
    raratop(ii)=b(rara(ii-1)+1);
    rarabot(ii)=b(length(b));
    %%%%so now raratop=[4 11 16], rarabot=[7 12 17];

    for ii=1:length(raratop)
      bm=raratop(ii); bp=rarabot(ii);
      if ((bm > 1)  & (bp < length(voinew)))
        bm=bm-1;
        bp=bp+1;
        bee=[bm:bp];
        slope=voinew(bp)-voinew(bm); slope=slope/(outwave(bp)-outwave(bm));
        voinew(bee)=slope*(outwave(bee)-outwave(bm))+voinew(bm);
      elseif ((bm == 1)  & (bp < length(voinew)))
        bm=bm;
        bp=bp+1;
        bee=[bm:bp];
        slope=voinew(bp)-voinew(bm); slope=slope/(outwave(bp)-outwave(bm));
        voinew(bee)=slope*(outwave(bee)-outwave(bm))+voinew(bm);
      elseif ((bm > 1)  & (bp == length(voinew)))
        bm=bm-1;
        bp=bp;
        bee=[bm:bp];
        slope=voinew(bp)-voinew(bm); slope=slope/(outwave(bp)-outwave(bm));
        voinew(bee)=slope*(outwave(bee)-outwave(bm))+voinew(bm);
        end

%orig stuff
%      start=rarat(ii); stop=rarat(ii+1);
%      bm=b(start)-1;
%      bp=b(stop)+1;
%      slope=voinew(bp)-voinew(bm); slope=slope/(bp-bm);
%      voinew(b)=slope*(b-bm)+voinew(bm);
      end

    %%%if things are still messed, force naughty k's to 0.0
    b=find(voinew < 0);
    voinew(b) = 0.0;

   if quiet > 0
      plot(outwave,exp(-voinew),outwave(b),exp(-voinew(b)),'r+')
      title('after getting rid of negative ks');
      end

    end
  end

%figure(1)
%%%%%plot(outwave,voinew,outwave,voi,'+',outwave,voi); grid
%plot(outwave,voinew,outwave,voi); grid; pause

%figure(2)
%plot(outwave,exp(-voinew),outwave,exp(-voi)); grid;
%title('after getting rid of negative ks'); pause
%error('hoho')
