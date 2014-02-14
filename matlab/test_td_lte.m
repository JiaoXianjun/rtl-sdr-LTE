% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% test_td_lte.m
% Try to demodulate China TD-LTE signal in 1.4MHz bandwidth(1.92Msps) like Lte-Cell-Scanner (written by James Peroulas).
% Some scripts are borrowed from:
% https://github.com/Evrytania/LTE-Cell-Scanner
% https://github.com/Evrytania/Matlab-Library
% https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

% example: test_td_lte('../scan-capture/frequency-2635-2655MHz/f2645_s1.92_g0_1s.bin')
% Those bin files are captured by rtl_sdr:
% rtl_sdr -f 900e6 -s 1.92e6 -n 1.92e6 tmp.bin
% -f 900MHz
% -s 1.92MHz sampling rate
% -n 1.92e6 samples
% result in tmp.bin with size 2*1.92e6 bytes because two bytes (I&Q) for one sample.

% See also README in root directory and ../scan-capture.

% function test_td_lte(rtl_sdr_bin_filename)

clear all;
close all;
% rtl_sdr_bin_filename = '../scan-capture/frequency-2635-2655MHz-know-PPM/f2645_s1.92_g0_SamplingPPM26.2_1s.bin';
% rtl_sdr_bin_filename = '../scan-capture/frequency-2635-2655MHz-know-PPM/f2645_s1.92_g20_SamplingPPM26.2_1s.bin';
rtl_sdr_bin_filename = '../scan-capture/frequency-2635-2655MHz/f2645_s1.92_g0_1s.bin';
% rtl_sdr_bin_filename = '../scan-capture/frequency-2595-2615MHz/f2604.9_s1.92_g0_1s.bin';
% rtl_sdr_bin_filename = '../scan-capture/frequency-2575-2595MHz/f2585_s1.92_g0_1s.bin';

num_radioframe = 5; % each radio frame length 10ms. MIB period is 4 radio frame

sampling_rate = 1.92e6; % LTE spec. 30.72MHz/16.

num_sym_per_slot = 7; % normal CP
num_slot_per_subframe = 2;
num_subframe_per_radioframe = 10;

len_time_subframe = 1e-3; % 1ms. LTE spec
len_time_radioframe = len_time_subframe*num_subframe_per_radioframe;
len_time_total = len_time_radioframe*num_radioframe;

[fd_pss, td_pss] = pss_gen;

num_sample = len_time_total*sampling_rate;
r = read_rtl_sdr_bin2IQ(rtl_sdr_bin_filename, num_sample);

% LTE channle filter for 6 RB
coef = fir1(46, (0.18e6*6+16*15e3)/sampling_rate);
% freqz(coef, 1, 1024);
% channel filter
r = filter(coef, 1, r);

% ppm_val = 116.4426;
ppm_val = 26.2009;
r_timing_correct = sampling_frequency_correction(r, ppm_val);

% % because MMDS LNB is used. No relationship between timing and frequency
% % try to correct timing firstly. Then to use LTE-Cell-Scanner with
% % timing-frequency relationship avoiding.
% [position, pss_idx, r_timing_correct] = PSS_detection_correction(r, fd_pss, td_pss);

capbuf = r_timing_correct.';
% freq_start = 2e9;
% ppm = 100;
% n_extra=floor((freq_start*ppm/1e6+2.5e3)/5e3);
% f_search_set = (-n_extra*5e3) : 5e3 : (n_extra*5e3);
% f_search_set = -100e3 : 5e3 : 100e3;
f_search_set = -95e3 : 5e3 : -85e3;
DS_COMB_ARM = 2;
fc = inf;

[xc_incoherent_collapsed_pow xc_incoherent_collapsed_frq n_comb_xc n_comb_sp xc_incoherent_single xc_incoherent sp_incoherent xc sp]= ...
xcorr_pss(capbuf,f_search_set,DS_COMB_ARM,fc);

FS_LTE = 30720000;
thresh1_n_nines=12;
R_th1=chi2cdf_inv(1-(10.0^(-thresh1_n_nines)), 2*n_comb_xc*(2*DS_COMB_ARM+1));
rx_cutoff=(6*12*15e3/2+4*15e3)/(FS_LTE/16/2);
Z_th1=R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*DS_COMB_ARM+1);
peaks=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set);