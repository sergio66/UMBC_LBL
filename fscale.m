
% fscale -- divide up octaves into integer scale steps
%
% fscale takes a start frequency f and a "chunk" count n and
% divides the octave f to 2*f into n intervals with integer
% boundaries, with the boundaries as close as possible to the
% more regular f*2^(j/n) spacing.  (When such a division is
% not possible, a warning message is printed.)  fscale plots
% the regular and integer spacings, and prints tables of the
% chunk boundaries and resulting chunk dv.
%
% H Motteler
% 22 Jan 01

f1 = 64;    % default start frequency
n1 = 16;    % default chunks per octave
noct1 = 6;  % default number of octaves/tables
k = 10000;  % tabulation points per chunk

f1 = 48;
noct1 = 1;


f = input(sprintf('start frequency [%g] > ', f1));
if isempty(f), f = f1; end

n = input(sprintf('chunks per octave [%g] > ', n1));
if isempty(n), n = n1; end

noct = input(sprintf('number of octaves [%g] > ', noct1));
if isempty(noct), noct = noct1; end

% integer grid over the octave
a = f : 2*f;  

% 2^(j/n) "scale" grid over the octave
b = f *  2.^((0:n)./n);

% discretized interval distances
db = round(diff(b));

if ~isempty(find(db == 0))
  fprintf(1, 'WARNING: interval distance rounds to 0\n');
end

% rebuild the scale with integer steps
c = f + [0,cumsum(db)];

% warn if it doesn't fit
if c(length(c)) ~= 2*f
  fprintf(1, 'WARNING: last discretized point %g != %g\n', c(length(c)), 2*f);
end

% warn if last dv too big
dv1 = (c(2) - c(1)) / k;
dv2 = (c(n+1) - c(n)) / k;
if dv2 > 2*dv1 
  fprintf(1, 'WARNING: last chunks too big\n');
end

% plot the optimal and discretized scales
bpts = ones(length(b),1) * .6;
cpts = ones(length(c),1) * .5;
apts = ones(length(a),1) * .4;

plot(b, bpts, 'g+', c, cpts, 'b+', a, apts, 'r+')
axis([f,2*f,0,1]);
legend('ideal 2\^(j/n)', 'discretized scale', 'integer grid');
title(sprintf('%d step scale', n))
grid

% print out chunk tables, assuming k-point chunks
fprintf(1, '\n       octave chunk tabulation\n\n');
fprintf(1, ' chunk  size   start    end      dv \n');

for j = 1:noct
  p = 2^(j-1);
  fprintf (1, ' ------------- octave %d -------------\n', j);
  for i = 1 : n
    dv = (c(i+1) - c(i)) / k;
    fprintf (1, ' %4d  %4d   %6.1f  %6.1f   %g\n', ...
             i, db(i)*p, c(i)*p, c(i+1)*p, dv*p);
  end
end

