clear all
addpath /asl/matlib/h4tools
addpath /home/sergio/git/SPECTRA

frtp = '/umbc/xfs2/strow/asl/s1/sergio/home/git/KCARTA/WORK/feb2002_raw_op_day_airs.rad.constemiss.rtp';
  f1 = 2230; f2 = 2235; iProf = 1;
  topts.tempname = '/asl/s1/sergio/JUNK/junk_wc.txt';

[h,ha,p,pa] = rtpread(frtp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /umbc/xfs2/strow/asl/s1/sergio/home/git/SPECTRA
rtpname = '/umbc/xfs2/strow/asl/s1/sergio/home/git/KCARTA/WORK/feb2002_raw_op_day_airs.rad.constemiss.rtp';
[w,dwater] = run8_rtp(1,f1,f2,rtpname,topts,iProf);
[w,do3]    = run8_rtp(3,f1,f2,rtpname,topts,iProf);
[w,dco2]   = run8_rtp(2,f1,f2,rtpname,topts,iProf);

semilogy(w,sum(dwater,1),w,sum(do3,1),w,sum(dco2,1));
legend('WV','O3','CO2','location','best','fontsize',12);
