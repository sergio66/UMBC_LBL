function line = show_unc(wv1,wv2,gid,HITRAN,HorG,index2show);

%% lineH16 = show_unc(1275,1325,6,2016,1,2);    HITRAN2016
%% lineG15 = show_unc(1275,1325,6,2015,-1,2);   GEISA 2015

%     wnum_unc_index   = str2num(unc_index(:,1));    line center
%     stren_unc_index  = str2num(unc_index(:,2));    line strength
%     abroad_unc_index = str2num(unc_index(:,3));    air broadening
%     sbroad_unc_index = str2num(unc_index(:,4));    self broadening
%     abcoef_unc_index = str2num(unc_index(:,5));    temp dependence of broadening
%     tsp_unc_index    = str2num(unc_index(:,6));    pressure shift

%% hgload /home/sergio/PAPERS/CONFERENCES/HITRAN/2018/Figs/UNC/strengthplus_airs.fig

if nargin == 3
  HITRAN = 2016;
  HorG = +1;      %% HITRAN not GEISA
  index2show = 1;
end
if nargin == 4
  HorG = +1;      %% HITRAN not GEISA
  index2show = 1;
end
if nargin == 5
  index2show = 1;
end

[iYes,line] = findlines_plot(wv1,wv2,gid,HITRAN,HorG);

unc_index = line.ai;
if index2show < 1 | index2show > 6
  error('can only show indices 1 .. 6')
elseif index2show == 1
  str = 'wnum';
  xindex = str2num(unc_index(:,1));
elseif index2show == 2
  str = 'stren';
  xindex = str2num(unc_index(:,2));  
elseif index2show == 3
  str = 'abroad';
  xindex = str2num(unc_index(:,3));  
elseif index2show == 4
  str = 'sbroad';
  xindex = str2num(unc_index(:,4));  
elseif index2show == 5
  str = 'abcoef = T dependence';
  xindex = str2num(unc_index(:,5));    
elseif index2show == 6
  str = 'press shift';
  xindex = str2num(unc_index(:,6));    
end

figure(1); scatter(line.wnum,xindex,20,line.iso,'filled'); grid on
  xlabel('wavnumber cm-1'); ylabel([str ' unc index']); title('colorbar = isotope'); colormap jet; colorbar

figure(2); scatter(line.wnum,line.stren,20,xindex,'filled'); grid on
  xlabel('wavnumber cm-1'); ylabel('stren'); title(['colorbar = ' str ' unc index']); colormap jet; colorbar
  set(gca,'yscale','log')