fid = fopen('gN_ir_list.txt');
A = fscanf(fid, '%s', [1 inf]);
fclose(fid);

n = length(A)/9;
B = reshape(A,9,n);
B = B';
for ii = 1 : 33
  fprintf(1,'%s \n',B(ii,:))
  if mod(ii,11) == 0
    disp(' ')
  end
end

p = randperm(n);
plot(sort(p)-(1:n))

iRemoveSome = input('Do you want to remove any gas? -1 for no, +N for yes : ');
if iRemoveSome > 0
  iRemoveSome = input('Enter range of gaseses to remove : ');
  iRemoveSome = iRemoveSome(1) : iRemoveSome(2);
end

iCnt = 0;
fid2 = fopen('gN_ir_list_permute.txt','w');
for ii = 1 : n
  boo = B(p(ii),:);
  if iRemoveSome < 0
    iCnt = iCnt + 1;
    fprintf(fid2,'%s \n',B(p(ii),:));
  else
    booGas = str2num(boo(1:2));
    if length(intersect(booGas,iRemoveSome)) == 0
      iCnt = iCnt + 1;
      fprintf(fid2,'%s \n',B(p(ii),:));
    end
  end
end
fclose(fid2);

if iRemoveSome > 0
  fprintf(1,'originally had %6i finally had %6i \n',n,iCnt);
end
