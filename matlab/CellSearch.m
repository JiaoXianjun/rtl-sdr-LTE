% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% CellSearch.m
% Improved LTE-Cell-Scanner (written by James Peroulas: https://github.com/Evrytania/LTE-Cell-Scanner).
% Add 1) TD-LTE; 2) external mixer (no assumption on relationship between sampling and carrier error) support
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

% Some scripts are borrowed from:
% https://github.com/Evrytania/LTE-Cell-Scanner
% https://github.com/Evrytania/Matlab-Library
% https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

% See also README in root directory and ../scan-capture.

clear all;
close all;

% freq_set = [2564.9e6] - 1998e6;
freq_set = [2645e6] - 1998e6;
% freq_set = [2585e6] - 1998e6;
% freq_set = [2604.9e6] - 1998e6;
% freq_set = 1890e6;
use_file_flag = 1; % set to 1 to use pre-captured file; set to 0 to use live dongle IQ samples (run "rtl_tcp -p 1234 -d 0" in shell first!)
% ------------------------------------------------------------------------------------
% rtl_sdr_bin_filename = '../scan-capture/frequency-1850-1880MHz/f1860_s1.92_g0_1s_strong.bin';% 17.3611PPM; -45kHz
% rtl_sdr_bin_filename = '../scan-capture/frequency-1850-1880MHz/f1860_s1.92_g0_1s.bin'; % 17.3611PPM; -45kHz
rtl_sdr_bin_filename = '../scan-capture/frequency-1880-1900MHz/f1890_s1.92_g0_1s.bin'; % -40kHz, 14.881ppm
% rtl_sdr_bin_filename = '../scan-capture/frequency-2555-2575MHz/f2564.9_s1.92_g0_1s.bin'; % -35kHz, 116.4426PPM
% rtl_sdr_bin_filename = '../scan-capture/frequency-2635-2655MHz-know-PPM/f2645_s1.92_g0_SamplingPPM26.2_1s.bin'; % -89531.0772Hz, 26.2009PPM
% rtl_sdr_bin_filename = '../scan-capture/frequency-2635-2655MHz-know-PPM/f2645_s1.92_g20_SamplingPPM26.2_1s.bin'; % -89531.0772Hz, 26.2009PPM
% rtl_sdr_bin_filename = '../scan-capture/frequency-2635-2655MHz/f2645_s1.92_g0_1s.bin'; % -89531.0772Hz, 26.2009PPM
% rtl_sdr_bin_filename = '../scan-capture/frequency-2575-2595MHz/f2585_s1.92_g0_1s.bin'; % -87975.8941Hz, 26.2009PPM
% rtl_sdr_bin_filename = '../scan-capture/frequency-2595-2615MHz/f2605_s1.92_g0_1s.bin'; % -65kHz, 116.4426PPM
% Bin file can be captured by rtl_sdr:
% rtl_sdr -f 900e6 -s 1.92e6 -n 1.92e6 tmp.bin
% -f 900MHz
% -s 1.92MHz sampling rate
% -n 1.92e6 samples
% result in tmp.bin with size 2*1.92e6 bytes because two bytes (I&Q) for one sample.
% you'd better capture 10s signal, and cut the last 1s signal by:
% extract_part_from_rtl_sdr_bin('sig_10s.bin',9*2*1.92e6,1*2*1.92e6,'sig_1s.bin')
% This ensures signal is stable enough.
% ------------------------------------------------------------------------------------

% If external mixer (like MMDS LNB) is used, this must be set to 0!
sampling_carrier_twist = 0; % ATTENTION! If this is 1, make sure fc is aligned with rtl_sdr_bin_filename!!!
gain = 0; % when external MMDS LNB is used, consider using fixed gain instead of AGC
% gain = 49.6; %seems that a fixed high gain is better

num_radioframe = 8; % each radio frame length 10ms. MIB period is 4 radio frame
sampling_rate = 1.92e6; % LTE spec. 30.72MHz/16.
num_subframe_per_radioframe = 10;
len_time_subframe = 1e-3; % 1ms. LTE spec
num_sample_per_radioframe = num_subframe_per_radioframe*len_time_subframe*sampling_rate;
num_sample = num_radioframe*num_sample_per_radioframe;
% LTE channle filter for center 6 RB
coef = fir1(46, (0.18e6*6+16*15e3)/sampling_rate); % freqz(coef, 1, 1024);
% coef = fir1(128, (0.16e6*6)/sampling_rate); % freqz(coef, 1, 1024);

DS_COMB_ARM = 2;
FS_LTE = 30720000;
thresh1_n_nines=12;
rx_cutoff=(6*12*15e3/2+4*15e3)/(FS_LTE/16/2);
THRESH2_N_SIGMA = 3;

[~, td_pss] = pss_gen;

if use_file_flag == 0
    num_dongle = 1;
    % check if previous tce objects existed. if so clear them
    if ~isempty(who('tcp_obj'))
        for i=1:length(tcp_obj)
            fclose(tcp_obj{i});
            delete(tcp_obj{i});
        end
        clear tcp_obj;
    end
    % construct tcp objects
    tcp_obj = cell(1, num_dongle);
    for i=1:num_dongle
        tcp_obj{i} = tcpip('127.0.0.1', 1233+i); % for dongle i
    end
    % set some parameters to tcp objects, and open them.
    for i=1:num_dongle
        set(tcp_obj{i}, 'InputBufferSize', 2*num_sample);
        set(tcp_obj{i}, 'Timeout', 2);
    end
    for i=1:num_dongle
        fopen(tcp_obj{i});
    end
    % set gain
    for i=1:num_dongle
        set_gain_tcp(tcp_obj{i}, gain*10); %be careful, in rtl_sdr the 10x is done inside C program, but in rtl_tcp the 10x has to be done here.
    end
    % set sampling rate
    for i=1:num_dongle
        set_rate_tcp(tcp_obj{i}, sampling_rate);
    end
    % set different start freq to different dongle
    for i=1:num_dongle
        set_freq_tcp(tcp_obj{i}, freq_set(1));
    end
    % read and discard to flush
    for i=1:num_dongle
        fread(tcp_obj{i}, 2*num_sample, 'uint8');
    end
end

if use_file_flag == 1
    loop_size = 1;
else
    loop_size = length(freq_set);
    real_count = zeros(1, num_dongle);
    s = zeros(2*num_sample, num_dongle);
    f_search_set = -100e3:5e3:100e3; % default frequency offset searching range if sampling_carrier_twist==1
end
peaks_store = cell(1,loop_size);
detect_flag_store = cell(1,loop_size);
tdd_flags_store = cell(1,loop_size);
r_all = zeros(num_sample, loop_size);
for freq_idx = 1 : loop_size
    if use_file_flag == 1
        r = read_rtl_sdr_bin2IQ(rtl_sdr_bin_filename, num_sample);
        if sampling_carrier_twist==0
            fc = inf;
        else
            fc = 1890e6; % Be careful! This must be aligned with your captured file!
        end
    else
        fc = freq_set(freq_idx);
        while 1 % read data at current frequency until success
            for i=1:num_dongle
                set_freq_tcp(tcp_obj{i}, fc); % set different frequency to different dongle
            end
            for i=1:num_dongle
                fread(tcp_obj{i}, 2*num_sample, 'uint8'); % flush to wait for its stable
            end
            for i=1:num_dongle
                fread(tcp_obj{i}, 2*num_sample, 'uint8'); % flush to wait for its stable
            end
            for i=1:num_dongle
                [s(:, i), real_count(i)] = fread(tcp_obj{i}, 2*num_sample, 'uint8'); % read samples from multi-dongles
            end

            if sum(real_count-(2*num_sample)) ~= 0
                disp(num2str([idx 2*num_sample, real_count]));
            else
                fid = fopen(['f' num2str(fc/1e6) '_s1.92_g' num2str(gain) '_' num2str(num_sample/sampling_rate) 's.bin'], 'w');
                fwrite(fid, s(:,1), 'uint8'); % only dongle 1 is stored!
                fclose(fid);
                r = raw2iq(s(:,1));
                break;
            end
        end
    end
    r_all(:,freq_idx) = r;
end
if use_file_flag == 0
    % close TCP
    for i=1:num_dongle
        fclose(tcp_obj{i});
    end
    for i=1:num_dongle
        delete(tcp_obj{i});
    end
    clear tcp_obj;
end

plot(real(r)); drawnow;

for freq_idx = 1 : loop_size
    r = r_all(:,freq_idx);
    
    disp(['Input averaged abs: ' num2str( mean(abs([real(r); imag(r)])) )]);

    r = filter(coef, 1, r);
    
    if use_file_flag == 1
        if sampling_carrier_twist==0
            fc = inf;
        else
            fc = 1890e6; % Be careful! This must be aligned with your captured file!
        end
        disp(['Processing  at ' rtl_sdr_bin_filename]);
    else
        fc = freq_set(freq_idx);
        disp(['Processing  at ' num2str(fc/1e6) 'MHz']);
    end

    if sampling_carrier_twist == 0
        disp('sampling_ppm_f_search_set_by_pss: try ... ... ');
        [period_ppm, f_search_set] = sampling_ppm_f_search_set_by_pss(r, td_pss);
        if period_ppm == inf
            disp('No valid PSS is found at pre-proc phase! Please try again.');
            peaks = [];
            detect_flag = [];
            tdd_flags = [];
            continue;
        else
            r = sampling_period_correction(r, period_ppm);
%             r = sampling_frequency_correction(r, period_ppm);
        end
%         [period_ppm, f_search_set] = sampling_ppm_f_search_set_by_pss(r, td_pss);
%         r = sampling_period_correction(r, period_ppm);
%         [period_ppm, f_search_set] = sampling_ppm_f_search_set_by_pss(r, td_pss);
%         r = sampling_period_correction(r, period_ppm);
%         [period_ppm, f_search_set] = sampling_ppm_f_search_set_by_pss(r, td_pss);
%         r = sampling_period_correction(r, period_ppm);
    end

    capbuf = r.';
    [xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, n_comb_xc, n_comb_sp, xc_incoherent_single, xc_incoherent, sp_incoherent, xc, sp]= ...
    xcorr_pss(capbuf,f_search_set,DS_COMB_ARM,fc, sampling_carrier_twist);

    R_th1=chi2inv(1-(10.0^(-thresh1_n_nines)), 2*n_comb_xc*(2*DS_COMB_ARM+1));
%     R_th1=chi2cdf_inv(1-(10.0^(-thresh1_n_nines)), 2*n_comb_xc*(2*DS_COMB_ARM+1));
    Z_th1=R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*DS_COMB_ARM+1);

    peaks=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set);
    
    detect_flag = zeros(1, length(peaks));
    tdd_flags = zeros(1, length(peaks));
    for i=1:length(peaks)
        for tdd_flag=0:1
            [peak, sss_h1_np_est, sss_h2_np_est, sss_h1_nrm_est, sss_h2_nrm_est, sss_h1_ext_est, sss_h2_ext_est]=...
                sss_detect(peaks(i),capbuf,THRESH2_N_SIGMA,fc,sampling_carrier_twist,tdd_flag);
            if ~isnan( peak.n_id_1 )
                break;
            end
        end
        if isnan( peak.n_id_1 )
            continue;
        end
        peak=pss_sss_foe(peak,capbuf,fc,sampling_carrier_twist,tdd_flag);
        [tfg, tfg_timestamp]=extract_tfg(peak,capbuf,fc,sampling_carrier_twist);
        [tfg_comp, tfg_comp_timestamp, peak]=tfoec(peak,tfg,tfg_timestamp,fc,sampling_carrier_twist);
        peak=decode_mib(peak,tfg_comp);
        if isnan( peak.n_rb_dl)
            continue;
        end
        if tdd_flag == 1
            disp('  Detected a TDD cell!');
        else
            disp('  Detected a FDD cell!');
        end
        if use_file_flag == 1
            disp(['  at ' rtl_sdr_bin_filename]);
        else
            disp(['  at ' num2str(fc/1e6) 'MHz']);
        end
        disp(['    cell ID: ' num2str(peak.n_id_cell)]);
        disp(['    RX power level: ' num2str(10*log10(peak.pow))]);
        disp(['    residual frequency offset: ' num2str(peak.freq_superfine)]);
        peaks(i) = peak;
        detect_flag(i) = 1;
        tdd_flags(i) = tdd_flag;
    end
    peaks_store{freq_idx} = peaks;
    detect_flag_store{freq_idx} = detect_flag;
    tdd_flags_store{freq_idx} = tdd_flags;
    if sum(detect_flag)==0
        disp('No LTE cells were found...');
    end
end

% show all Cells information
disp(' ');
disp('-------------------------------Cells information summary-------------------------------');
for freq_idx = 1 : loop_size
    if use_file_flag == 1
        disp(['At ' rtl_sdr_bin_filename]);
    else
        fc = freq_set(freq_idx);
        disp(['At ' num2str(fc/1e6) 'MHz']);
    end
    peaks = peaks_store{freq_idx};
    detect_flag = detect_flag_store{freq_idx};
    tdd_flags = tdd_flags_store{freq_idx};
    if isempty(detect_flag)
        disp('No valid PSS is found at pre-proc phase! Please try again.');
    else
        if sum(detect_flag)
            hit_idx = find(detect_flag);
            for i=1:length(hit_idx);
                peak = peaks(hit_idx(i));
                tdd_flag = tdd_flags(hit_idx(i));
                if tdd_flag == 1
                    cell_mode_str = 'TDD';
                else
                    cell_mode_str = 'FDD';
                end
                disp(['Cell ' num2str(i) ' information:--------------------------------------------------------']);
                disp(['            Cell mode: ' num2str(cell_mode_str)]);
                disp(['              Cell ID: ' num2str(peak.n_id_cell)]);
                disp(['   Num. eNB Ant ports: ' num2str(peak.n_ports)]);
                disp(['    Carrier frequency: ' num2str(fc/1e6) 'MHz']);
                disp(['Residual freq. offset: ' num2str(peak.freq_superfine/1e3) 'kHz']);
                disp(['       RX power level: ' num2str(10*log10(peak.pow))]);
                disp(['              CP type: ' peak.cp_type]);
                disp(['              Num. RB: ' num2str(peak.n_rb_dl)]);
                disp(['       PHICH duration: ' peak.phich_dur]);
                disp(['  PHICH resource type: ' num2str(peak.phich_res)]);
            end
        else
            disp('No LTE cells were found...  Please try again.');
        end
    end
end