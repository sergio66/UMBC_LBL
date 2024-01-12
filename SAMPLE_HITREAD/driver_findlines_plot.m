wv1 = 2000;  %% start wavenumber
wv2 = 2200;  %% end wavenumber
gid = 5;     %% CO

line160 = [];

fnamePRE = '/home/sergio/SPECTRA/SAMPLE_HITREAD/g';

fnamePOST    = '.dat';
fnameIN      = int2str(gid);
hitlin_fname = [fnamePRE fnameIN fnamePOST];
fprintf(1,'in findlines_plot.m we are opening %s for lineparams for gid %2i \n',hitlin_fname,gid)

start = wv1;
stop  = wv2;

line160 = hitread(start,stop,0,gid,hitlin_fname);
line100 = translate2oldHITparams(line160);

if line100.linct > 0
  semilogy(line160.wnum,line160.stren,'bo',line100.wnum,line100.stren,'r.'); title(['gid = ' num2str(gid)]);
    xlabel('Wavenumber cm^{-1}'); ylabel('log_{10}(Stren)'); hl = legend('160','100','location','best');
end
