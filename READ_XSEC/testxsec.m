
gid = 64;
dv = .0045;
v1 = 930;
v2 = v1 + 20 - dv;

% [k1, fr1] = calc_xsec(gid, v1, v2, dv, 280, 800, 1);
% [k2, fr2] = calc_xsecV1(gid, v1, v2, dv, 280, 800, 2);
[k1, fr1] = calc_xsec(gid, v1, v2, dv, 220, 600, 1);
[k2, fr2] = calc_xsecV1(gid, v1, v2, dv, 220, 600, 2);

figure(3)
plot(fr1,k1,fr2,k2)

legend('new', 'old');

sum(k1 ~= k2)

