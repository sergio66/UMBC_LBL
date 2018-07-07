
dofudge = +1;      %%%%this is when we want to fit the data
dofudge = -1;      %%%%this are the TRUE run7co2 results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% if you wanna use only RAL results,   set dofudge = +1 here, 
%%%                                          dofudge = -1 in driver4um.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% if you wanna fit beta/doc,           set dofudge = +1 here, 
%%%                                          dofudge = +1 in driver4um.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% if you wanna use RAL/JJOHNS results, set dofudge = -1 here, 
%%%                                          dofudge = -1 in driver4um.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dofudge = +1;      %%%%this is to fit the data, or only use OLD RAL results
dofudge = -1;      %%%%this is RAL results + JJOHNS blend

if (dofudge == +1)
  disp('yrun_sigsig : dofudge = +1 ==> fitting data');
  end

dofudge = -1; 

if ((band == 2350) & (dofudge == -1))
  %%%%% 2350 needs 4 sets of line mix parameters : one in MAIN entire region
  %%%%% (from above in "loader(temperature,band,path_length,ptotal,pself,prb)")
  %%%%% then one in the 2385-2390 region, then one in the 2390-2395 region,
  %%%%% then one in the 2395-2405 region
    %%% JULY 4, 2002
      %%% region == 1 ==> 2380-2391
      %%% region == 2 ==> 2391-2400
    %%% JULY 7, 2002
      %%% region == 1 ==> 2380-2391
      %%% region == 2 ==> 2391-2395
      %%% region == 3 ==> 2395-2405
  [stJJOHNS]=loader_jjohns(temperature,band,path_length,ptotal,pself,prb,1);
  stuff.beta_pureJJ(1)       = stJJOHNS.beta_pure;
  stuff.beta_forJJ(1)        = stJJOHNS.beta_for;
  stuff.frequency_shiftJJ(1) = stJJOHNS.frequency_shift;
  stuff.betaJJ(1)            = stJJOHNS.beta;
  stuff.durationJJ(1)        = stJJOHNS.duration;
  [stJJOHNS]=loader_jjohns(temperature,band,path_length,ptotal,pself,prb,2);
  stuff.beta_pureJJ(2)       = stJJOHNS.beta_pure;
  stuff.beta_forJJ(2)        = stJJOHNS.beta_for;
  stuff.frequency_shiftJJ(2) = stJJOHNS.frequency_shift;
  stuff.betaJJ(2)            = stJJOHNS.beta;
  stuff.durationJJ(2)        = stJJOHNS.duration;
  [stJJOHNS]=loader_jjohns(temperature,band,path_length,ptotal,pself,prb,3);
  stuff.beta_pureJJ(3)       = stJJOHNS.beta_pure;
  stuff.beta_forJJ(3)        = stJJOHNS.beta_for;
  stuff.frequency_shiftJJ(3) = stJJOHNS.frequency_shift;
  stuff.betaJJ(3)            = stJJOHNS.beta;
  stuff.durationJJ(3)        = stJJOHNS.duration;
elseif ((band == 2351) & (dofudge == -1))
  %%%%% 2351 needs 2 sets of line mix parameters : one in entire region,
  %%%%% (from above in "loader(temperature,band,path_length,ptotal,pself,prb)")
  %%%%% then separate one in the 2280-2290 region
  %% july 10, 2002 : use region 2 parameters .. not very good
  %%[stJJOHNS]=loader_jjohns(temperature,band,path_length,ptotal,pself,prb,2);
  %% july 11, 2002 : use region 1 parameters
  [stJJOHNS]=loader_jjohns(temperature,band,path_length,ptotal,pself,prb,1);
  for ii = 1:3
    stuff.beta_pureJJ(ii)       = stJJOHNS.beta_pure;
    stuff.beta_forJJ(ii)        = stJJOHNS.beta_for;
    stuff.frequency_shiftJJ(ii) = stJJOHNS.frequency_shift;
    stuff.betaJJ(ii)            = stJJOHNS.beta;
    stuff.durationJJ(ii)        = stJJOHNS.duration;
    end
else
  %%%%% others 1 sets of line mix parameters : one in entire region
  %%%%% (from above in "loader(temperature,band,path_length,ptotal,pself,prb)")
  for ii = 1:3
    stuff.beta_pureJJ(ii)       = stuff.beta_pure;
    stuff.beta_forJJ(ii)        = stuff.beta_for;
    stuff.frequency_shiftJJ(ii) = stuff.frequency_shift;
    stuff.betaJJ(ii)            = stuff.beta;
    stuff.durationJJ(ii)        = stuff.duration;
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
