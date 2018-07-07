function [freqs,num_lines,strens]=run1(gasID,fmin,fmax);

%[line]=hittomat2(low,high,vers,strengthM,gasID);

if (gasID == 1)
  vers=12;
else
  vers=10;
  end

[line]=hittomat2(fmin,fmax,vers,1e-30,gasID);

freqs=fmin:25:fmax-25;
len=length(freqs);
for ii=1:len
  %now sort the lines into near and far wing lines
  [total,number]=sortwings_stren(line,freqs(ii),freqs(ii)+25,0);
  num_lines(ii)=total;
  strens(ii,:)=number;
  fprintf(1,'\n freq = %7.1f, num lines = %7i',freqs(ii),num_lines(ii));
  end

%plot(freqs,num_lines,':',freqs,sum(strens'))
plot(freqs,num_lines,':',freqs,strens)

fprintf(1,'\n ');

%strengths(1)=0        1e-26
%strengths(2)=1e-26    1e-25
%...
%strengths(10)=1e-15   1e-14
%strengths(11)=1e-14   1

linestrens=[-30 -26 -25 -24 -23 -22 -21 -20 -19 -18 -17 -16]
