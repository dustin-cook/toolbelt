clear all
close all
clc

input_dir = [pwd filesep 'inputs'];


ground_motion_set = table;
eq_id = 1;
i_files = dir([input_dir filesep '*.AT2']);
for i = 1:length(i_files)
    ground_motion_set.id(i,1) = i;
    ground_motion_set.set_id(i,1) = eq_id;
    filename = i_files(i).name;
    eq_name = strrep(filename,'.AT2','');
    output_dir = [pwd  filesep 'outputs' filesep eq_name];
    mkdir(output_dir)
    [ output_txt, pga, dt, eq_length  ] = fn_reformat_NGA_gm( filename, input_dir, output_dir );
    if rem(i,2) % odd pair
        ground_motion_set.pair(i,1) = 1;
    else % even pair
        ground_motion_set.pair(i,1) = 2;
        eq_id = eq_id + 1;
    end
    
    ground_motion_set.eq_dt(i,1) = dt;
    ground_motion_set.eq_length(i,1) = eq_length;
    ground_motion_set.eq_name{i,1} = eq_name;
    ground_motion_set.pga(i,1) = pga;
end

writetable(ground_motion_set,[pwd  filesep 'outputs' filesep 'ground_motion_set.csv'])

% gm_data_NGA = gm_data(3:end);
% 
% mygm_dir = 'C:\Users\dtc2\Desktop\Repos\Opensees\ground_motions\FEMA_far_field\1_MULHOLLAND009';
% mygm_name = '1_MULHOLLAND009.tcl';
% 
% fileID = fopen([mygm_dir filesep mygm_name],'r');
% gmTH = fscanf(fileID,'%f');
% fclose(fileID);
% 
% plot(gm_data_NGA)
% hold on
% plot(gmTH)
% close
% mean(gm_data_NGA./gmTH')
% 
% plot(gm_data_NGA./gmTH')
% 
% max(abs(gm_data_NGA))
% max(abs(gmTH))