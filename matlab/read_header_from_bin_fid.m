function [fc_requested, fc_programmed, fs_requested, fs_programmed] = read_header_from_bin_fid(fid)

fc_requested = inf;
fc_programmed = inf;
fs_requested = inf;
fs_programmed = inf;

% fid = fopen(bin_filename);
% 
% if fid==-1
%     disp('read_header_from_bin: Can not open file for read!');
%     return;
% end

magic1 = fread(fid, 1, 'double');
tmp1 = fread(fid, 1, 'uint64');

magic2 = fread(fid, 1, 'double');
tmp2 = fread(fid, 1, 'uint64');

magic3 = fread(fid, 1, 'double');
tmp3 = fread(fid, 1, 'uint64');

magic4 = fread(fid, 1, 'double');
tmp4 = fread(fid, 1, 'uint64');

magic5 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

magic6 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

magic7 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

magic8 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

% fclose(fid);

fc_requested_magic = 73492.215;
fc_programmed_magic = -0.7923597;
fs_requested_magic = -189978508;
fs_programmed_magic = 93.126712;

reserve1_magic = -53243.129;
reserve2_magic = 0.0008123898;
reserve3_magic = -6.0098321;
reserve4_magic = 237.09983;

if magic1 == fc_requested_magic && ...
   magic2 == fc_programmed_magic && ...
   magic3 == fs_requested_magic && ...
   magic4 == fs_programmed_magic && ...
   magic5 == reserve1_magic && ...
   magic6 == reserve2_magic && ...
   magic7 == reserve3_magic && ...
   magic8 == reserve4_magic

    fc_requested = tmp1;
    fc_programmed = tmp2;
    fs_requested = tmp3;
    fs_programmed = tmp4;
end