clear 
close all
clc

%% Notes
% Roof Nodes: 6170, 6171, 6176, 6177

%% Define Inputs Directory
gm_set_name = 'ICSB_analysis';
g_factor = 1; % denominator to convert raw data to g
dt_idx = 1; % location of the dt value
eq_data_idx = 2; % location of the start of the ground motion
plot_spectra = 1;

%% Run Function for each GM in inputs folder
% Load all files in Input folder and convert them to single column signals
% with dt as the first entry
input_dir = [pwd filesep 'inputs' filesep gm_set_name];
files = [dir([input_dir filesep '*.AT2']),dir([input_dir filesep '*.tcl']),dir([input_dir filesep '*.txt'])];
for i = 1:length(files)
    gm_set.id(i,1) = i;
    % Grab GM signal
    if contains(files(i).name,'.AT2') % NGA file
        gm_name = erase(files(i).name,'.AT2');
        [ output_dir ] = fn_make_dir( [pwd filesep 'outputs' filesep gm_set_name filesep gm_name] );
        [ gm_data ] = fn_reformat_NGA_gm( files(i).name, input_dir, output_dir );
    elseif contains(files(i).name,'.txt') % Preformatted TCL signal
        gm_name = erase(files(i).name,'.txt');
        [ output_dir ] = fn_make_dir( [pwd filesep 'outputs' filesep  gm_set_name filesep gm_name] );
%         [ gm_data ] = fn_reformat_txt_gm( files(i).name, input_dir, output_dir );
        [ ag, dt ] = fn_reformat_txt_gm_2( files(i).name, input_dir, output_dir, g_factor, eq_data_idx, dt_idx );
    else
        error('Ground Motion Data Structure Not Recognized')
    end
    
    % Define GM set info
%     gm_set.set_id(i,1) = str2double(regexp(gm_name,'^\d*(?=(_))','match'));
%     if sum(gm_set.set_id == gm_set.set_id(i,1)) == 1
%         gm_set.pair(i,1) = 1;
%     elseif sum(gm_set.set_id == gm_set.set_id(i,1)) == 2
%         gm_set.pair(i,1) = 2;
%     else
%         error()
%     end
%     gm_set.eq_dt(i,1) = dt;
%     gm_set.eq_length(i,1) = length(ag);
%     gm_set.eq_name{i,1} = gm_name;
%     gm_set.pga(i,1) = max(abs(ag));

    % Calculate Spectra
    [ spectra ] = fn_single_spectra( dt, ag );

    % Save Spectra
    writetable(spectra,[output_dir filesep 'spectra.csv'])

    % Plot Spectra
    if plot_spectra
        figure
        hold on
        plot(spectra.period,spectra.psa_1,'LineWidth',1.5,'DisplayName','1% Damping') 
        plot(spectra.period,spectra.psa_2,'LineWidth',1.5,'DisplayName','2% Damping') 
        plot(spectra.period,spectra.psa_3,'LineWidth',1.5,'DisplayName','3% Damping') 
        plot(spectra.period,spectra.psa_5,'LineWidth',1.5,'DisplayName','5% Damping') 
        xlabel('Period (s)')
        ylabel('PSa (g)')
        plot_name = 'spectra';
        fn_format_and_save_plot( output_dir, plot_name, 1 )
    end
    
end

% Save ground motion set info
% gm_set_table = struct2table(gm_set);
% writetable(gm_set_table,[pwd filesep 'outputs' filesep gm_set_name filesep 'ground_motion_set.csv'])


