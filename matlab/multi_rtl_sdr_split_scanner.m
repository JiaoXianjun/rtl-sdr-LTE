% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% multi_rtl_sdr_split_scanner.m
% Frequency band scanning via multiple rtl-sdr dongles. Each dongle for each sub-band to speedup!
% A script of project: https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

% Assume that you have installed rtl-sdr
% (http://sdr.osmocom.org/trac/wiki/rtl-sdr) and have those native utilities run correctly already.

% For example, you have multiple dongles, please run multiple rtl_tcp in multiple shell respectively as:
% rtl_tcp -p 1234 -d 0
% rtl_tcp -p 1235 -d 1
% rtl_tcp -p 1236 -d 2
% ...

% Then run this script in MATLAB.

% ATTENTION! In some computer, every time before you run this script, maybe you need to terminate multiple rtl_tcp and re-launch them again.
% ATTENTION! Please reduce number of inspected points by reducing frequency range or increasing step size, if your computer hasn't enough memory installed. Because all signlas are stored firstly before processing.

% Change following parameters as you need:

% Number of dongles you have connected to your computer
num_dongle = 1; % more dongles, much faster.

% Beginning of the band you are interested in
start_freq = 502e6; % for test
% start_freq = 935e6; % Beginning of Primary GSM-900 Band downlink
% start_freq = (1575.42-15)*1e6; % GPS L1
% start_freq = (1207.14-30)*1e6; % COMPASS B2I
% start_freq = (1561.098-30)*1e6; % COMPASS B1I
% start_freq = 1.14e9; % Beginning of GNSS(GPS/GLONASS/COMPASS/Galileo) Band

% End of the band you are interested in
end_freq = 702e6; % for test
% end_freq = 960e6; % End of Primary GSM-900 Band downlink
% end_freq = (1575.42+30)*1e6; % GPS L1
% end_freq = (1207.14+30)*1e6; % COMPASS B2I
% end_freq = (1561.098+30)*1e6; % COMPASS B1I
% end_freq = 1.63e9; % Beginning of GNSS(GPS/GLONASS/COMPASS/Galileo) Band

freq_step = 0.5e6; % less step, higher resolution, narrower FIR bandwidth, slower speed

observe_time = 0.2; % observation time at each frequency point. ensure it can capture your signal!

RBW = freq_step; % Resolution Bandwidth each time we inspect

gain = 20; % If this is larger than 0, the fixed gain will be set to dongles

% use high sampling rate and FIR to improve estimation accuracy
sample_rate = 2.048e6; % sampling rate of dongles

coef_order = (2^(ceil(log2(sample_rate/RBW))))-1;
coef_order = min(coef_order, 127);
coef_order = max(coef_order, 31);
coef = fir1(coef_order, RBW/sample_rate);
% freqz(coef, 1, 1024);

num_samples = observe_time*sample_rate;

clf;
close all;

% construct freq set for each dongle
freq_orig = start_freq:freq_step:end_freq;
num_freq_per_sub_band = ceil(length(freq_orig)/num_dongle);
num_pad = num_freq_per_sub_band*num_dongle - length(freq_orig);
freq = [freq_orig freq_orig(end)+(1:num_pad).*freq_step];
freq = vec2mat(freq, num_freq_per_sub_band);

real_count = zeros(1, num_dongle);
s_all = uint8( zeros(2*num_samples, num_freq_per_sub_band*num_dongle) );
decimate_ratio = floor(sample_rate/(2*RBW));

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
    set(tcp_obj{i}, 'InputBufferSize', 8*2*num_samples);
    set(tcp_obj{i}, 'Timeout', 60);
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
    set_rate_tcp(tcp_obj{i}, sample_rate);
end

% set different start freq to different dongle
for i=1:num_dongle
    set_freq_tcp(tcp_obj{i}, freq(i,1));
end

% read and discard to flush
for i=1:num_dongle
    fread(tcp_obj{i}, 8*2*num_samples, 'uint8');
end

% capture samples of all frequencies firstly!
tic;
for freq_idx = 1:num_freq_per_sub_band
    while 1 % read data at current frequency until success
        for i=1:num_dongle
            set_freq_tcp(tcp_obj{i}, freq(i, freq_idx)); % set different frequency to different dongle
        end
        for i=1:num_dongle
            [s_all(:, freq_idx + (i-1)*num_freq_per_sub_band), real_count(i)] = fread(tcp_obj{i}, 2*num_samples, 'uint8'); % gather data from different dongles to s_all
        end

        if sum(real_count-(2*num_samples)) ~= 0
            disp(num2str([idx 2*num_samples, real_count]));
        else
            break;
        end
    end
end
e = toc;
ideal_time_cost = observe_time*num_freq_per_sub_band;

% close TCP
for i=1:num_dongle
    fclose(tcp_obj{i});
end
for i=1:num_dongle
    delete(tcp_obj{i});
end
clear tcp_obj;

disp('Scanning done!');
disp(['actual time cost ' num2str(e) ' ideal cost ' num2str(ideal_time_cost) ' efficiency(ideal/actual) ' num2str(ideal_time_cost/e)]);
disp(' ');
disp('Begin process ...');

% generate power spectrum
tic;
r = raw2iq( double( s_all ) ); % remove DC. complex number constructed.
r_flt = filter(coef, 1, r);% filter target band out
power_spectrum = mean(abs(r_flt(1:decimate_ratio:end, :)).^2, 1);% get averaged power
e1 = toc;
disp(['time cost ' num2str(e1) ' scan/process ' num2str(e/e1)]);
disp(['total time cost ' num2str(e1+e)]);

% plot power spectrum (converted to dB)
freq_linear = freq';
freq_linear = freq_linear(:)';
figure;
format_string = {'b.-', 'r.-', 'k.-', 'm.-'};
for i=1:num_dongle
    plot(freq(i,:).*1e-6, 10.*log10(power_spectrum( (i-1)*num_freq_per_sub_band+1: i*num_freq_per_sub_band)), format_string{i}); hold on;
end

legend_string = cell(1, num_dongle);
for i=1:num_dongle
    legend_string{i} = ['dongle ' num2str(i)];
end
legend(legend_string);

filename = ['split_scan_' num2str(start_freq) '_' num2str(end_freq) '_gain' num2str(gain) '_' num2str(num_dongle) 'dongles.mat'];
save(filename, 'power_spectrum', 'start_freq', 'end_freq', 'freq_step', 'observe_time', 'RBW', 'gain', 'sample_rate', 'coef', 'freq');
