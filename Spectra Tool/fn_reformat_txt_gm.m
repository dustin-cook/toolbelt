function [ output ] = fn_reformat_txt_gm( file_name, input_dir, output_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Load GM File
txt_data = fileread([input_dir filesep file_name]);

data_raw = str2double(strsplit(txt_data,' '));
data = data_raw(~isnan(data_raw));

% Find G Conversion
g_factor = data(1);

% Find dt
dt = data(2);

% Crop data to just Earthquake and translate to data
eq_data = data(3:end);
eq_data_g = eq_data/g_factor;

% Loop through data to create single line signal
output = [dt,eq_data_g];
fileID = fopen([output_dir filesep erase(file_name,'.txt') '.tcl'],'w');
for j = 1:length(output)
    fprintf(fileID,'%d \n',output(j));
end
fclose(fileID);


end

