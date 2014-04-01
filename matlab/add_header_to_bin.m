function add_header_to_bin(bin_filename, fc_requested, fc_programmed, fs_requested, fs_programmed)

[fc_requested_exist, fc_programmed_exist, fs_requested_exist, fs_programmed_exist] = read_header_from_bin(bin_filename);
if fc_requested_exist ~= inf
    disp('There is already a header! Just return.');
    disp(['fc_requested ' num2str(fc_requested_exist) ' fc_programmed ' num2str(fc_programmed_exist)  ' fs_requested ' num2str(fs_requested_exist)  ' fs_programmed ' num2str(fs_programmed_exist) ]);
    return;
end

fid = fopen(bin_filename, 'r');

if fid==-1
    disp('add_header_to_bin: Can not open file for read!');
    return;
end

tmp_store = fread(fid, inf, 'uint8');
fclose(fid);

fc_requested_magic = 73492.215;
fc_programmed_magic = -0.7923597;
fs_requested_magic = -189978508;
fs_programmed_magic = 93.126712;

reserve1_magic = -53243.129;
reserve2_magic = 0.0008123898;
reserve3_magic = -6.0098321;
reserve4_magic = 237.09983;

fid = fopen(bin_filename, 'w');

if fid==-1
    disp('add_header_to_bin: Can not open file for write!');
    return;
end

fwrite(fid, fc_requested_magic, 'double');
fwrite(fid, fc_requested, 'uint64');

fwrite(fid, fc_programmed_magic, 'double');
fwrite(fid, fc_programmed, 'uint64');

fwrite(fid, fs_requested_magic, 'double');
fwrite(fid, fs_requested, 'uint64');

fwrite(fid, fs_programmed_magic, 'double');
fwrite(fid, fs_programmed, 'uint64');

fwrite(fid, reserve1_magic, 'double');
fwrite(fid, 0, 'uint64');

fwrite(fid, reserve2_magic, 'double');
fwrite(fid, 0, 'uint64');

fwrite(fid, reserve3_magic, 'double');
fwrite(fid, 0, 'uint64');

fwrite(fid, reserve4_magic, 'double');
fwrite(fid, 0, 'uint64');

fwrite(fid, tmp_store, 'uint8');

fclose(fid);
