clear all;
close all;

fid = fopen('f2564.9_s1.92_g0_10s.bin');
fopen(fid);
s = fread(fid, inf, 'uint8');
fclose(fid);
s = s((9*1.92e6*2+1):end);

fid = fopen('f2564.9_s1.92_g0_1s.bin', 'w');
fopen(fid);
fwrite(fid, s, 'uint8');
fclose(fid);

