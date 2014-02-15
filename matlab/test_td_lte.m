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
rtl_sdr_bin_filename = '../scan-capture/frequency-1880-1900MHz/f1890_s1.92_g0_1s.bin'; % -40kHz, 14.881ppm

% rtl_sdr_bin_filename = '../scan-capture/frequency-2555-2575MHz/f2564.9_s1.92_g0_1s.bin'; % -35kHz, 116.4426PPM

% rtl_sdr_bin_filename = '../scan-capture/frequency-2635-2655MHz-know-PPM/f2645_s1.92_g0_SamplingPPM26.2_1s.bin'; % -89531.0772Hz, 26.2009PPM
% rtl_sdr_bin_filename = '../scan-capture/frequency-2635-2655MHz-know-PPM/f2645_s1.92_g20_SamplingPPM26.2_1s.bin'; % -89531.0772Hz, 26.2009PPM
% rtl_sdr_bin_filename = '../scan-capture/frequency-2635-2655MHz/f2645_s1.92_g0_1s.bin'; % -89531.0772Hz, 26.2009PPM

% rtl_sdr_bin_filename = '../scan-capture/frequency-2575-2595MHz/f2585_s1.92_g0_1s.bin'; % -87975.8941Hz, 26.2009PPM

% rtl_sdr_bin_filename = '../scan-capture/frequency-2595-2615MHz/f2605_s1.92_g0_1s.bin'; % -65kHz, 116.4426PPM

carrier_frequency = 2605;

num_radioframe = 10; % each radio frame length 10ms. MIB period is 4 radio frame

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

% % % because MMDS LNB is used. No relationship between timing and frequency
% % % try to correct timing firstly. Then to use LTE-Cell-Scanner with
% % % timing-frequency relationship avoiding.
% [position, pss_idx, r_timing_correct] = PSS_detection_correction(r, fd_pss, td_pss);

% ppm_val = 133.9286;
% ppm_val = 116.4426;
% ppm_val = 26.2009;
ppm_val = 14.881;
r_timing_correct = sampling_frequency_correction(r, ppm_val);

capbuf = r_timing_correct.';
% freq_start = 2e9;
% ppm = 100;
% n_extra=floor((freq_start*ppm/1e6+2.5e3)/5e3);
% f_search_set = (-n_extra*5e3) : 5e3 : (n_extra*5e3);
f_search_set = -40e3;
% f_search_set = -100e3 : 5e3 : -50e3;
DS_COMB_ARM = 2;
fc = inf;

% if isempty(dir('xcorr_pss.mat'))
    [xc_incoherent_collapsed_pow xc_incoherent_collapsed_frq n_comb_xc n_comb_sp xc_incoherent_single xc_incoherent sp_incoherent xc sp]= ...
    xcorr_pss(capbuf,f_search_set,DS_COMB_ARM,fc);
%     save xcorr_pss_result.mat xc_incoherent_collapsed_pow xc_incoherent_collapsed_frq n_comb_xc n_comb_sp xc_incoherent_single xc_incoherent sp_incoherent xc sp capbuf f_search_set DS_COMB_ARM fc; 
% else
%     load xcorr_pss.mat;
% end

FS_LTE = 30720000;
thresh1_n_nines=12;
R_th1=chi2cdf_inv(1-(10.0^(-thresh1_n_nines)), 2*n_comb_xc*(2*DS_COMB_ARM+1));
rx_cutoff=(6*12*15e3/2+4*15e3)/(FS_LTE/16/2);
Z_th1=R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*DS_COMB_ARM+1);

% if isempty(dir('peak_search.mat'))
    peaks=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set);
%     save peak_search.mat peaks xc_incoherent_collapsed_pow xc_incoherent_collapsed_frq Z_th1 f_search_set;
% else
%     load peak_search.mat;
% end

THRESH2_N_SIGMA = 3;
detect_flag = zeros(1, length(peaks));
for i=1:length(peaks)
    peak = peaks(i);
    [peak sss_h1_np_est sss_h2_np_est sss_h1_nrm_est sss_h2_nrm_est sss_h1_ext_est sss_h2_ext_est]=sss_detect(peak,capbuf,THRESH2_N_SIGMA,fc);
    if isnan( peak.n_id_1 )
        continue;
    end
    peak=pss_sss_foe(peak,capbuf,fc);
    [tfg tfg_timestamp]=extract_tfg(peak,capbuf,fc);
    [tfg_comp tfg_comp_timestamp peak]=tfoec(peak,tfg,tfg_timestamp,fc);
    peak=decode_mib(peak,tfg_comp);
    if isnan( peak.n_rb_dl)
        continue;
    end
    disp('  Detected a cell!');
    disp(['    cell ID: ' num2str(peak.n_id_cell)]);
    disp(['    RX power level: ' num2str(10*log10(peak.pow))]);
    disp(['    residual frequency offset: ' num2str(peak.freq_superfine)]);
    peaks(i) = peak;
    detect_flag(i) = 1;
end

if sum(detect_flag)
    hit_idx = find(detect_flag);
    for i=1:length(hit_idx);
        disp(' ');
        disp(['Cell ' num2str(i) ' information:']);
        peak = peaks(hit_idx(i));
        disp(['              Cell ID: ' num2str(peak.n_id_cell)]);
        disp(['   Num. eNB Ant ports: ' num2str(peak.n_ports)]);
        disp(['    Carrier frequency: ' num2str(carrier_frequency) 'MHz']);
        disp(['Residual freq. offset: ' num2str(peak.freq_superfine) 'Hz']);
        disp(['       RX power level: ' num2str(10*log10(peak.pow))]);
        disp(['              CP type: ' peak.cp_type]);
        disp(['              Num. RB: ' num2str(peak.n_rb_dl)]);
        disp(['       PHICH duration: ' peak.phich_dur]);
        disp(['  PHICH resource type: ' num2str(peak.phich_res)]);
        disp([' Sampling freq. error: ' num2str(ppm_val) 'PPM']);
    end
else
    disp('No TD-LTE cells were found...');
end