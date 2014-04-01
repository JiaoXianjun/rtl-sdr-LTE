% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% read specific number of IQ samples from rtl-sdr captured bin file
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function s = read_rtl_sdr_bin2IQ(rtl_sdr_bin_filename, num_sample)

header_exist = false;

[fc_requested_exist, fc_programmed_exist, fs_requested_exist, fs_programmed_exist] = read_header_from_bin(rtl_sdr_bin_filename);

if fc_requested_exist~=inf
    header_exist = true;
    disp('There is already a header!');
    disp(['fc_requested ' num2str(fc_requested_exist) ' fc_programmed ' num2str(fc_programmed_exist)  ' fs_requested ' num2str(fs_requested_exist)  ' fs_programmed ' num2str(fs_programmed_exist) ]);
end

fid = fopen(rtl_sdr_bin_filename);

if fid==-1
    disp('read_rtl_sdr_bin2IQ: Can not open input file!');
    return;
end

if header_exist
    read_header_from_bin_fid(fid);
end

[s, count] = fread(fid, 2*num_sample, 'uint8');
fclose(fid);

if count < (2*num_sample)
    disp('read_rtl_sdr_bin2IQ: Please ensure number of samples is sufficient in bin file!');
    s = -1;
    return;
end

s = raw2iq(s);
