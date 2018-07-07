pooh0 = fit_run7(X0,topts,labdata,sim);
poohF = fit_run7(X,topts,labdata,sim);
subplot(211); plot(f,t_rawdata); grid
subplot(212); plot(f,pooh0,f,poohF,'r'); grid
filesave = ['save file' NI '.mat X labdata sim'];
eval([filesave])

