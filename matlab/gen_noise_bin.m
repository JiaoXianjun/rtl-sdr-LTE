% function gen_noise_bin
clear all;
close all;

len = 1920000*2;
fid = fopen('noise.bin', 'w');
fwrite(fid, randn(1, len)*10 + 128, 'uint8');
fclose(fid);

fid = fopen('noise.bin');
a = fread(fid, inf, 'uint8');
fclose(fid);
plot(a);
