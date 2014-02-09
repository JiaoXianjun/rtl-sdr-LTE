% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% test_td_lte_pss.m
% Try to identify if the signal I captured in Beijing is TD-LTE or not by using PSS correlation (written by James Peroulas).
% Some scripts are borrowed from:
% https://github.com/Evrytania/LTE-Cell-Scanner
% https://github.com/Evrytania/Matlab-Library
% https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

% example: test_td_lte_pss('../scan-capture/frequency-2635-2655MHz/f2645_s1.92_g0_1s.bin')
% Those bin files are captured by rtl_sdr:
% rtl_sdr -f 900e6 -s 1.92e6 -n 1.92e6 tmp.bin
% -f 900MHz
% -s 1.92MHz sampling rate
% -n 1.92e6 samples
% result in tmp.bin with size 2*1.92e6 bytes because two bytes (I&Q) for one sample.

% See also README in root directory and ../scan-capture.

% clear all;
% close all;

function test_td_lte(rtl_sdr_bin_filename)

sampling_rate = 1.92e6; % LTE spec. baseband sampling rate of 1.4MHz bandwidth option.

len_time = 40e-3;
num_samples = len_time*sampling_rate;

fid = fopen(rtl_sdr_bin_filename);
s = fread(fid, inf, 'uint8');
fclose(fid);

if length(s) < (2*num_samples)
    disp('Please capture signal with at least 40ms length!');
    return;
end

s = raw2iq(s(1:num_samples));

f_search_set = [-100e3:25e3:100e3];
ds_comb_arm = 3;
fc = (647+1998)*1e6;

[xc_incoherent_collapsed_pow xc_incoherent_collapsed_frq n_comb_xc n_comb_sp xc_incoherent_single xc_incoherent sp_incoherent xc sp]=xcorr_pss(s.',f_search_set,ds_comb_arm,fc);

for i=1:3
    figure(i);
    for idx = 1:length(f_search_set)
        subplot(3,3,idx);
        plot(abs(xc_incoherent(i,:,idx)));
        title(['pss(' num2str(i-1) ')' ' freq offset ' num2str(f_search_set(idx)/1e3) ' kHz']);
    end
end

% plot(abs(s).^2);
