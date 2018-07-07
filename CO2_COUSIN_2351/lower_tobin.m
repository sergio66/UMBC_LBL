function diff=lower(f_in,red_num)

n_in=length(f_in);
n_out=fix(n_in/red_num);

%for n=1:n_out
%  f_out(n)=f_in(n*red_num);
%  end
%diff=f_out';

%fprintf(1,'n_in,n_out = %6i %6i \n',n_in,n_out);
nout=1:n_out;
f_out=f_in(nout*red_num);
diff=f_out;
