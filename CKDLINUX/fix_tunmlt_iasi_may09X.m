ll = load('tunmlt_iasi_may09.txt');
ff = ll(:,2); con = ll(:,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iA = find(ll(:,2) == 2662); iB = find(ll(:,2) == 2664);
slope = (con(iB)-con(iA))/(ff(iB)-ff(iA));

conNew = con;
woof = find(ff > 2663 & ff < 2666);
conNew(woof) = slope*(ff(woof)-ff(iB)) + con(iB);

plot(ll(:,2),ll(:,5),'bo-',ff,conNew,'r*-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iA = find(ll(:,2) == 2190); iB = find(ll(:,2) == 2200);
slope = (con(iB)-con(iA))/(ff(iB)-ff(iA));

woof = find(ff > 2163 & ff < 2185);
conNew(woof) = slope*(ff(woof)-ff(iB)) + con(iB);

plot(ll(:,2),ll(:,5),'bo-',ff,conNew,'r*-'); grid on

ll(:,5) = conNew;
fid = fopen('tunmlt_iasi_may09X.dat','w');
fprintf(fid,'%5i %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f \n',ll')
fclose(fid)
