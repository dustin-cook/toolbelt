function [] = main_calc_spectra(gm_set_name, g_factor, dt_idx, eq_data_idx, plot_spectra)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
        [ output_dir ] = fn_make_dir( ['outputs' filesep gm_set_name filesep gm_name] );
        [ gm_data ] = fn_reformat_NGA_gm( files(i).name, input_dir, output_dir );
        ag.(['dir_' num2str(i)]) = gm_data(3:end);
        dt.(['dir_' num2str(i)]) = gm_data(2);
    elseif contains(files(i).name,'.txt') % Preformatted TCL signal
        gm_name = erase(files(i).name,'.txt');
        [ output_dir ] = fn_make_dir( ['outputs' filesep  gm_set_name filesep gm_name] );
%         [ gm_data ] = fn_reformat_txt_gm( files(i).name, input_dir, output_dir );
        [ ag, dt ] = fn_reformat_txt_gm_3( files(i).name, input_dir, output_dir, g_factor, dt_idx );
%         [ ag.(['dir_' num2str(i)]) ] = fn_reformat_txt_gm( files(i).name, input_dir, output_dir );
%         dt.(['dir_' num2str(i)]) = 0.01;
    else
        error('Ground Motion Data Structure Not Recognized')
    end
    
    % Define GM set info
%     gm_set.eq_dt(i,1) = dt.(['dir_' num2str(i)]);
%     gm_set.eq_length(i,1) = length(ag.(['dir_' num2str(i)]));
%     gm_set.eq_name{i,1} = gm_name;
%     gm_set.pga(i,1) = max(abs(ag.(['dir_' num2str(i)])));
    gm_set.eq_dt(i,1) = dt;
    gm_set.eq_length(i,1) = length(ag);
    gm_set.eq_name{i,1} = gm_name;
    gm_set.pga(i,1) = max(abs(ag));

    % Calculate Spectra
    [ spectra ] = fn_single_spectra( dt/sqrt(3), ag );

    % Save Spectra
    writetable(spectra,[output_dir filesep 'spectra.csv'])
    spectra_table.T = spectra.T;
    spectra_table.(['gm_' gm_name ]) = spectra.psa_005;
    
    % Plot Spectra
    if plot_spectra
        figure
        hold on
        plot(spectra.T,spectra.psa_001,'LineWidth',1.5,'DisplayName','1% Damping') 
        plot(spectra.T,spectra.psa_002,'LineWidth',1.5,'DisplayName','2% Damping') 
        plot(spectra.T,spectra.psa_003,'LineWidth',1.5,'DisplayName','3% Damping') 
        plot(spectra.T,spectra.psa_005,'LineWidth',1.5,'DisplayName','5% Damping') 
        xlabel('Period (s)')
        ylabel('PSa (g)')
        plot_name = ['spectra'];
        fn_format_and_save_plot( output_dir, plot_name, 1 )
    end    
end

%% Create Max Direction Spectra
% if isfield(ag,'dir_2')
%     [ spectra ] = fn_single_spectra( dt.dir_1, ag.dir_1, ag.dir_2 );
% 
%     % Save Spectra
%     writetable(spectra,['outputs' filesep  gm_set_name filesep 'spectra_max_dir.csv'])
% end
% 
%% Save ground motion set info
gm_set_table = struct2table(gm_set);
writetable(gm_set_table,[pwd filesep 'outputs' filesep gm_set_name filesep 'ground_motion_set.csv'])

% Save aggregate 5% damped spectra table
spectra_table = struct2table(spectra_table);
writetable(spectra_table,[pwd filesep 'outputs' filesep gm_set_name filesep 'all_gms_5%_spectra.csv'])
end

