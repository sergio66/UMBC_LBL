f1 = 2300; f2 = 2600;
f1 = 2375; f2 = 2425;
f1 = 2375; f2 = 2400;

gid = input('Enter gid : 2/4/8/10 (default 2) : ');
if length(gid) == 0
  gid = 2;
end

if gid == 2
  f1 = 2300; f2 = 2600;
  f1 = 2375; f2 = 2425;
  f1 = 2375; f2 = 2400;
  f1 = 2375; f2 = 2425;
elseif gid == 4
  f1 = 2400; f2 = 2600;
end

[w0,d0] = run8(gid,f1,f2,'IPFILES/co2four_hitemp');
semilogy(w0,d0);

if gid == 2
  [wlm,dlm] = run8co2(gid,f1,f2,'IPFILES/co2four_hitemp');
end

%disp('ret to cont'); pause

topts.HITRAN = '/asl/data/hitran/HITEMP/h16.by.gas';
[wx,dx] = run8(gid,f1,f2,'IPFILES/co2four_hitemp',topts);

%{
if gid == 2
  saver = ['save testHITEMP_g' num2str(gid) '.mat d0 dx dlm gid f1 f2 topts w0 wx wlm'];
else
  saver = ['save testHITEMP_g' num2str(gid) '.mat d0 dx gid f1 f2 topts w0 wx'];
end
eval(saver)
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE
if gid == 2
  figure(1); semilogy(w0,d0,'--',wx,dx,'-');
  figure(1); 
    semilogy(w0,d0(1,:),'b',w0,d0(2,:),'g',w0,d0(3,:),'k',w0,d0(4,:),'r',...
           w0,dx(1,:),'b.-',w0,dx(2,:),'g.-',w0,dx(3,:),'k.-',w0,dx(4,:),'r.-')
    hl = legend('250','500','750','1000','location','best'); ylabel('OD')

  figure(1); 
    plot(w0,exp(-d0(1,:)),'b',w0,exp(-d0(2,:)),'g',w0,exp(-d0(3,:)),'k',w0,exp(-d0(4,:)),'r',...
         w0,exp(-dx(1,:)),'b.-',w0,exp(-dx(2,:)),'g.-',w0,exp(-dx(3,:)),'k.-',w0,exp(-dx(4,:)),'r.-')
    hl = legend('250','500','750','1000','location','best'); ylabel('OD')
    grid

  [fc,qc0] = quickconvolve(w0,d0,1,1); qc0 = qc0';
  [fc,qcx] = quickconvolve(wx,dx,1,1); qcx = qcx';
  [fc,qclm] = quickconvolve(wlm,dlm,1,1); qclm = qclm';

  figure(2); clf
    plot(fc,exp(-qc0(1,:)),'b',fc,exp(-qc0(2,:)),'g',fc,exp(-qc0(3,:)),'k',fc,exp(-qc0(4,:)),'r',...
         fc,exp(-qcx(1,:)),'b.-',fc,exp(-qcx(2,:)),'g.-',fc,exp(-qcx(3,:)),'k.-',fc,exp(-qcx(4,:)),'r.-',...
         fc,exp(-qclm(1,:)),'bx-',fc,exp(-qclm(2,:)),'gx-',fc,exp(-qclm(3,:)),'kx-',fc,exp(-qclm(4,:)),'rx-')
    hl = legend('250','500','750','1000','location','best'); ylabel('trans')

  figure(2); plot(w0,dx./d0); hl = legend('250','500','750','1000'); ylabel('hitemp/usual');
  figure(2); 
    subplot(221); ii = 1; semilogy(w0,d0(ii,:),'b',w0,d0(ii,:),'r'); title('250 K');
%      ax = axis; ax(1) = 2500; ax(2) = 2510; axis(ax)
    subplot(222); ii = 2; semilogy(w0,d0(ii,:),'b',w0,dx(ii,:),'r'); title('500 K');
%      ax = axis; ax(1) = 2500; ax(2) = 2510; axis(ax)
    subplot(223); ii = 3; semilogy(w0,d0(ii,:),'b',w0,dx(ii,:),'r'); title('750 K');
%      ax = axis; ax(1) = 2500; ax(2) = 2510; axis(ax)
    subplot(224); ii = 4; semilogy(w0,d0(ii,:),'b',w0,dx(ii,:),'r'); title('1000 K');
%      ax = axis; ax(1) = 2500; ax(2) = 2510; axis(ax)

elseif gid == 4
  figure(1); semilogy(w0,d0,'--',wx,dx,'-');
  figure(1); 
    semilogy(w0,d0(1,:),'b',w0,d0(2,:),'g',w0,d0(3,:),'k',w0,d0(4,:),'r',...
           w0,dx(1,:),'b.-',w0,dx(2,:),'g.-',w0,dx(3,:),'k.-',w0,dx(4,:),'r.-')
    hl = legend('250','500','750','1000','location','best'); ylabel('OD')
    axis([2450 2550 1e-4 1e-1])
    axis([2500 2520 1e-4 1e-1])
    axis([2500 2510 1e-4 1e-1])

  figure(2); plot(w0,dx./d0); hl = legend('250','500','750','1000'); ylabel('hitemp/usual');
  figure(2); 
    subplot(221); ii = 1; semilogy(w0,d0(ii,:),'b',w0,d0(ii,:),'r'); title('250 K');
      ax = axis; ax(1) = 2500; ax(2) = 2510; axis(ax)
    subplot(222); ii = 2; semilogy(w0,d0(ii,:),'b',w0,dx(ii,:),'r'); title('500 K');
      ax = axis; ax(1) = 2500; ax(2) = 2510; axis(ax)
    subplot(223); ii = 3; semilogy(w0,d0(ii,:),'b',w0,dx(ii,:),'r'); title('750 K');
      ax = axis; ax(1) = 2500; ax(2) = 2510; axis(ax)
    subplot(224); ii = 4; semilogy(w0,d0(ii,:),'b',w0,dx(ii,:),'r'); title('1000 K');
      ax = axis; ax(1) = 2500; ax(2) = 2510; axis(ax)
end

