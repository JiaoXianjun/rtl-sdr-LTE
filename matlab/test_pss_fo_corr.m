% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% analyze PSS detection under big frequency offset
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

% according to simulation result, gold configuration is that use whole PSS (length 137), and fo step 50kHz
% simulation result:
% len pss 137 fo step 50kHz SNR 0dB 10000 cost 29.9569s fail rate 0.003
% len pss 137 fo step 50kHz SNR 3dB 10000 cost 27.4461s fail rate 0
% len pss 137 fo step 25kHz SNR 0dB 10000 cost 40.6559s fail rate 0.0025
% len pss 137 fo step 25kHz SNR 3dB 10000 cost 41.6056s fail rate 0
% len pss 128 fo step 50kHz SNR 0dB 10000 cost 23.6994s fail rate 0.0091
% len pss 128 fo step 50kHz SNR 3dB 10000 cost 23.0074s fail rate 0.0002
% len pss 128 fo step 25kHz SNR 0dB 10000 cost 41.1991s fail rate 0.0078
% len pss 128 fo step 25kHz SNR 3dB 10000 cost 41.4224s fail rate 0

% function test_pss_fo_corr
clear all; 
close all;

[fd_pss, td_pss] = pss_gen;

pss_idx = 1;
% pss = td_pss(end-127:end, pss_idx);
pss = td_pss(:, pss_idx);
len_pss = length(pss);

fo_search_set = [-200e3:50e3:200e3];
% fo_search_set = 0;
pss_set = kron(ones(1, length(fo_search_set)), pss).*exp(1i.*2.*pi.*(1./1.92e6).*(0:(length(pss)-1)).'*fo_search_set);
fo_step = (fo_search_set(2) - fo_search_set(1))/1e3;

var_pss = mean(abs(td_pss(:, pss_idx)).^2);

r = zeros(1, 2*137);
tmp_total = zeros(1, length(fo_search_set));
num_test = 10000;
num_fail = 0;
snr = 3;
sigma2 = var_pss/(10^(snr/10));
correct_pos = 147 - (length(pss)-128);
tic;
for idx = 1:num_test
    r_pss = [(randn(137,1)+1i.*randn(137,1)).*sqrt(var_pss./2); td_pss(:, pss_idx); (randn(137,1)+1i.*randn(137,1)).*sqrt(var_pss./2)];
    r_pss = r_pss + sqrt(sigma2/2).*(randn(length(r_pss),1)+1i.*randn(length(r_pss),1));

    fo = (2*rand-1)*200e3;
%     fo = (-7.5e3/1) + 0*15e3;
    r_pss_fo = r_pss.*exp(1i.*2.*pi.*fo.*(1./1.92e6).*(0:(length(r_pss)-1))');

    for i=1:(2*137)
        s = r_pss_fo(i:(i+length(pss)-1));
        
%         %  %---------method 1-----------------------
%         for j=1:length(fo_search_set)
%             tmp = conj(pss_set(:,j)).*s;
%             tmp = vec2mat(tmp, 8);
%             tmp = sum( abs(sum(tmp, 2)).^2 );
% %             r(i) = r(i) + tmp;
%             tmp_total(j) = tmp;
%         end
%         r(i) = max(tmp_total);
%         %  %---------end of method 1-----------------------
        
%         %  %---------method 2-----------------------
%         r(i) = abs(s(1:64)'*s(128:-1:65))^2;
%         %  %---------end of method 2-----------------------
        
        %  %---------method 3-----------------------
        tmp = conj(pss).*s;
        r(i) = max(abs(fft(tmp)).^2);
        %  %---------end of method 3-----------------------
    end
    plot(r); drawnow;
    % plot(138, r(138), 'rs');
    
    if mod(idx, 100) == 0
        disp(['len pss ' num2str(len_pss) ' fo step ' num2str(fo_step) 'kHz SNR ' num2str(snr) 'dB ' num2str(idx) ' cost ' num2str(toc) 's fail rate ' num2str(num_fail/idx) ]);
        tic;
    end
    
    [~, max_idx] = max(r);
    if max_idx ~= correct_pos
%         plot(r); hold on;
%         plot(correct_pos, r(correct_pos), 'rs');
        num_fail = num_fail + 1;
        disp(['idx ' num2str(idx) ' max_idx ' num2str(max_idx) ' fail rate ' num2str(num_fail/idx)]);
%         break;
    end
end

disp(['fail rate ' num2str(num_fail/idx)]);
