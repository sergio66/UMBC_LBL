:
#copy file over and then take out first 3 lines

echo "cp ../../DATA/RefProf/refgasX refYgasX" > temp
echo "sed '1,3d' refYgasX > refgasX" > temp1

ll=1
llmax=63
while test $ll -le $llmax
do
  sed -e "s/X/$ll/g" temp  > infile
  sed -e "s/X/$ll/g" temp1 > infile1
  chmod +x infile
  chmod +x infile1
  infile
  infile1
  ll=`expr $ll + 1`
done
rm refYgas* temp temp1 infile infile1
