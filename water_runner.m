f0   = input('Enter start chunk : ');
fend = input('Enter end   chunk : ');
fe = f0 + 25;
fnameIN = 'IPFILES/newcamex_water_run6';
fnameIN = 'IPFILES/camex3_nast_run6';
fnameIN = 'IPFILES/wintex1_water_run6';

while f0 <= (fend-25)
  fprintf(1,'doing chunk %5i \n',f0);
  [fr,k00]=run6water(1,f0,fe,0.0005,0.1,0.5,1,1,2,25,5,0.0,0.0,'V',...
                     -1,1,1,+2,0,fnameIN);    %%%need local lineshape
  [fr,k00]=run7water(1,f0,fe,fnameIN);
  fnameOUT = '/taro/s1/sergio/HITRAN2000WATER/CAMEX1/water';
  fnameOUT = '/taro/s1/sergio/HITRAN2000WATER/CAMEX3/water';
  fnameOUT = '/taro/s1/sergio/HITRAN2000WATER/WINTEX/water';
  fnameOUT = [fnameOUT num2str(f0) '.mat'];
  saver = ['save ' fnameOUT ' fr k00'];
  eval([saver]);
  f0 = f0 + 25;
  fe = fe + 25;
  end
