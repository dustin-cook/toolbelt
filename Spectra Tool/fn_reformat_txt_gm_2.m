function [ eq_data_g, eq_dt ] = fn_reformat_txt_gm_2( file_name, input_dir, output_dir, g_factor, eq_data_idx, dt_idx )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Load GM File
txt_data = fileread([input_dir filesep file_name]);
txt_data_cell = strsplit(txt_data,'\n');
data_raw = str2double(txt_data_cell);
data = data_raw(~isnan(data_raw));

% Crop data to just Earthquake and translate to data
eq_dt = str2double(strrep(txt_data_cell(dt_idx),'dt = ',''));
eq_data = data_raw(eq_data_idx:end);
eq_data = eq_data(~isnan(eq_data));
eq_data_g = eq_data/g_factor;

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

