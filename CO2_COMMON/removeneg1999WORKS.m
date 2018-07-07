function [voinew]=removeneg(outwave,voi);
%this function makes sures all k's >= 0
%%%%this is the function that was used in the stuff done June 1999 - March 2000
%%%but it might have a problem, as it does the slopes based on indices,
%%%so if resolution is differnt, it does different slopes!!!
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
      slope=voinew(bp)-voinew(bm); slope=slope/(bp-bm);
      voinew(b)=slope*(b-bm)+voinew(bm);
    elseif ((b(1) == 1)  & (b(length(b)) < length(voinew)))
      %%%% eg k=[ -0.1 -0.5 -1 -0.5 0.5 1.5 10.5];
      bm=b(1);
      bp=b(length(b))+1;
      slope=voinew(bp)-voinew(bm); slope=slope/(bp-bm);
      voinew(b)=slope*(b-bm)+voinew(bm);
    elseif ((b(1) > 1)  & (b(length(b)) == length(voinew)))
      %%%% eg k=[5 2 1 -0.1 -0.5 -1 -0.5];
      bm=b(1)-1;
      bp=b(length(b));
      slope=voinew(bp)-voinew(bm); slope=slope/(bp-bm);
      voinew(b)=slope*(b-bm)+voinew(bm);
      end

    if quiet > 0
      plot(outwave,exp(-voinew),outwave(b),exp(-voinew(b)),'r+')
      title('after getting rid of negative ks');
      end

  else 
    %%%%the places where the k's are negative are not adjacent
    %break b into chunks ..........     
    fprintf(1,'spaced! \n');

    index=(bd==ones(size(bd)));
    rara=find((~index==1));
    rarat(1)=1; rarat(2:length(rara)+1)=rara; rarat(length(rara)+2)=length(b);

    for ii=1:length(rarat)-1
      start=rarat(ii); stop=rarat(ii+1);
      bm=b(start);    bp=b(stop);
      if ((bm > 1)  & (bp < length(voinew)))
        bm=bm-1;
        bp=bp+1;
        slope=voinew(bp)-voinew(bm); slope=slope/(bp-bm);
        voinew(b)=slope*(b-bm)+voinew(bm);
      elseif ((bm == 1)  & (bp < length(voinew)))
        bm=bm;
        bp=bp+1;
        slope=voinew(bp)-voinew(bm); slope=slope/(bp-bm);
        voinew(b)=slope*(b-bm)+voinew(bm);
      elseif ((bm > 1)  & (bp == length(voinew)))
        bm=bm-1;
        bp=bp;
        slope=voinew(bp)-voinew(bm); slope=slope/(bp-bm);
        voinew(b)=slope*(b-bm)+voinew(bm);
        end

%orig stuff
%      start=rarat(ii); stop=rarat(ii+1);
%      bm=b(start)-1;
%      bp=b(stop)+1;
%      slope=voinew(bp)-voinew(bm); slope=slope/(bp-bm);
%      voinew(b)=slope*(b-bm)+voinew(bm);
      end

    b=find(voinew < 0);
    voinew(b) = 0.0;

   if quiet > 0
      plot(outwave,exp(-voinew),outwave(b),exp(-voinew(b)),'r+')
      title('after getting rid of negative ks');
      end

    end
  end

plot(outwave,voinew,outwave,voi)
title('after getting rid of negative ks'); pause
plot(outwave,exp(-voinew),outwave,exp(-voi))
title('after getting rid of negative ks'); pause
