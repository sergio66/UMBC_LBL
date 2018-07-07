function [a,b]=quicklorentz(temp,freqq,strenqt,w_forq,w_selfq,stuff);

density=stuff.density; 
pressure_self=stuff.pressure_self; 
pressure_for=stuff.pressure_for; 
pressure_ref=stuff.pressure_ref; 
temperature_ref=stuff.temperature_ref; 
path_length=stuff.path_length; 
bsm=stuff.bsm; 
K_scale_mixing=density*pressure_self/pressure_ref*temperature_ref*... 
        path_length/temp/pi; 
ff=min(freqq)-50:0.01:max(freqq)+50;
lor=zeros(size(ff));
for i=1:length(freqq)
  w_tot=(pressure_self*w_selfq(i)+pressure_for*w_forq(i))/pressure_ref;
  lor=lor+strenqt(i)*(w_tot)./((ff-freqq(i)).^2 +(w_tot).^2); 
  end
lor=lor*K_scale_mixing;

vvv=find((ff <= min(freqq)-10)|(ff >= max(freqq)+10));
lor1=lor; lor1(vvv)=lor1(vvv)*0.5;
lor=exp(-lor); lor1=exp(-lor1);

figure(2)
semilogy(ff,lor,ff,lor1,ff(vvv),lor1(vvv)); grid
axis([min(ff) max(ff) 1e-6 1]);
figure(1)

vvv=find((lor1 >= 1e-4) & (ff <= min(freqq)-5));
a=ff(vvv(length(vvv)));
vvv=find((lor1 >= 1e-4) & (ff >= max(freqq)+5));
b=ff(vvv(1));
