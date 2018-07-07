function bands = plot_bands(gid);

bands = list_bands(gid);

for ii = 1 : length(bands.v1)
  line([bands.v1(ii) bands.v2(ii)],[0 0],'color','k','linewidth',2);
  plot([bands.v1(ii) bands.v2(ii)],[0 0],'ko','markersize',5)
end
