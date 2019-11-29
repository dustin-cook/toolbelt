function [ eq_data_g, eq_dt ] = fn_reformat_txt_gm_3( file_name, input_dir, output_dir, g_factor, dt_idx )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Load GM File and remove not GM signal lines
txt_data = fileread([input_dir filesep file_name]);
txt_data_cell = strsplit(txt_data,'\n');
data_raw = str2double(txt_data_cell);
eq_data = data_raw(~isnan(data_raw));
eq_data_g = g_factor*eq_data;

% Grab dt
eq_dt = str2double(strrep(txt_data_cell(dt_idx),'dt = ',''));

% Loop through data to create single line signal for running in opensees
fileID = fopen([output_dir filesep erase(file_name,'.txt') '.tcl'],'w');
for j = 1:length(eq_data_g)
    fprintf(fileID,'%d \n',eq_data_g(j));
end
fclose(fileID);

% Loop through data to create single line signal with g-factor and dt as txt doc
output_txt = [g_factor,eq_dt,eq_data_g];
fileID = fopen([output_dir filesep file_name],'w');
for j = 1:length(output_txt)
    fprintf(fileID,'%d \n',output_txt(j));
end
fclose(fileID);

end

