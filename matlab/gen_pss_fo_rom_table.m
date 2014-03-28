function a = gen_pss_fo_rom_table(f_search_set)
close all; 
clear all;
f_search_set = -60e3:5e3:55e3;
[~, td_pss] = pss_gen;
pss_fo_set = pss_fo_set_gen(td_pss, f_search_set);
pss_fo_set = pss_fo_set.';

% figure;
% subplot(2,3,1); plot(real(pss_fo_set(13,:))); hold on;
% subplot(2,3,2); plot(real(pss_fo_set(end,:))); hold on;
% subplot(2,3,3); plot(real(pss_fo_set(1,:))); hold on;
% subplot(2,3,4); plot(real(pss_fo_set(30,:))); hold on;
% subplot(2,3,5); plot(real(pss_fo_set(50,:))); hold on;
% subplot(2,3,6); plot(real(pss_fo_set(70,:))); hold on;

max_val1 = max(max(abs(real(pss_fo_set))));
max_val2 = max(max(abs(imag(pss_fo_set))));
max_val = max(max_val1, max_val2);

enlarge_ratio = ((2^15)-1)/max_val;

enlarge_ratio = floor( log2(enlarge_ratio) );

enlarge_ratio = 2^enlarge_ratio;


pss_fo_set = round(real(pss_fo_set).*enlarge_ratio) + 1i.*round(imag(pss_fo_set).*enlarge_ratio);

% figure;
% subplot(2,3,1); plot(real(pss_fo_set(13,:))); hold on;
% subplot(2,3,2); plot(imag(pss_fo_set(end,:))); hold on;
% subplot(2,3,3); plot(real(pss_fo_set(1,:))); hold on;
% subplot(2,3,4); plot(imag(pss_fo_set(30,:))); hold on;
% subplot(2,3,5); plot(real(pss_fo_set(50,:))); hold on;
% subplot(2,3,6); plot(imag(pss_fo_set(70,:))); hold on;

enlarge_ratio_log2 = log2(enlarge_ratio);

[num_pss, len_pss] = size(pss_fo_set);

fid = fopen('pss_fo_rom_table.txt', 'w');

fprintf(fid, '#define ENLARGE_RATIO (2^%d)\n\n', enlarge_ratio_log2);

fprintf(fid, 'constant short2 coef[%d] = { \\ \n', num_pss*len_pss);

for i=1:(num_pss-1)
    for j=1:len_pss
        fprintf(fid, '(short2)(%d,%d), ', real(pss_fo_set(i,j)), imag(pss_fo_set(i,j)));
    end
    fprintf(fid, '\\ \n');
end
i = num_pss;
for j=1:(len_pss-1)
    fprintf(fid, '(short2)(%d,%d), ', real(pss_fo_set(i,j)), imag(pss_fo_set(i,j)));
end
j=len_pss;
fprintf(fid, '(short2)(%d,%d)};\n\n', real(pss_fo_set(i,j)), imag(pss_fo_set(i,j)));

fclose(fid);
