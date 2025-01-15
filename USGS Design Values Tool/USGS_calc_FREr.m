%% Develpp FR and LS seismic design values
% uses data from USGS API calls loaded into matlab data structure to
% develop 2-point and multi-point design spectrum. FREr values are
% calculated using a risk targetted approah. Choose from the site set and
% code version options listed in the user inputs. All data is saved into a
% remote directory as specified.
clear all
close all
clc

%% User Inputs
% Set of sites for which to calc design values
% site_set = 'usgs_160';
site_set = 'nehrp_34';
% site_set = 'frtc_testbed';

% USGS NSHM and ASCE version
version_id = 'ASCE 7-22';
nshm = 2018;
% version_id = 'ASCE 7-28';
% nshm = 2023;

% Remote save directory
save_dir = ['C:\Users\dtc2\OneDrive - NIST\PUC\TS5\FREr\' site_set '\' version_id];
plot_dir = [save_dir '\Spectra_Plots'];
mkdir(plot_dir)

%% Optional Analysis Inputs
% Set risk targeted variables
rt_params.FR.beta = 0.6; 
rt_params.FR.PoE = 0.25;
rt_params.FR.rt = 0.1;
rt_params.FR.t = 50; % years
rt_params.FR.maxdir = false; % convert from RotD50 to RotD100

% ASCE 7-28 DE gm parameters
rt_params.LS.beta = 0.6; 
rt_params.LS.PoE = 0.1;
rt_params.LS.rt = 0.015;
rt_params.LS.t = 50; % years
rt_params.LS.maxdir = true; % convert from RotD50 to RotD100

% Set matlab colors
cmap = [0 0.4470 0.7410;
        0.8500 0.3250 0.0980;
        0.9290 0.6940 0.1250;
        0.4940 0.1840 0.5560;
        0.4660 0.6740 0.1880;
        0.3010 0.7450 0.9330;
        0.6350 0.0780 0.1840];
    
%% Imports
% functions
import usgs.*

% Previously saved Matlab data
if exist('DATA_usgs_api.mat','file')
    load('DATA_usgs_api.mat')
else
    error('Missing Data: Call USGS API')
end

% Filter out data to just sites fo interest
out_filt = out(~isnan([out.([site_set '_id'])]));
sites = table;

%% Load hazard data for RTGM
for s = 1:length(out_filt)
    s
    
    %Set Site values
    site = out_filt(s);
    sites.id(s,1) = site.([site_set '_id']);
    sites.name{s,1} = site.name;
    sites.lat(s,1) = site.lat;
    sites.lon(s,1) = site.lon;
    sites.region{s,1} = site.region;
    
    for v = 1:length(site.vs30)
        % Load pre-called values from the saved matlab data
        hazard_curve = site.hazard_curve.(['nshm_' num2str(nshm)]).(['vs30_' num2str(site.vs30(v))]);
        design_values_7_22 = site.design_values.(['site_class_' site.site_class{v}]); % These are strictly ASCE 7-22 design values
        
        %% Pull MPS
        % Develop design spectrum based on FR risk target values
        [rt_spectrum_FR] = fn_develop_RT_spectrum(hazard_curve, rt_params.FR);

        % Develop design spectrum based on LS risk target values
        if strcmp(version_id,'ASCE 7-28')
            % develop MPS from hazard call (fully probabilistic)
            [rt_spectrum_LS] = fn_develop_RT_spectrum(hazard_curve, rt_params.LS);
            rt_spectrum_LS.sa_DE = (2/3)*rt_spectrum_LS.sa;
        elseif strcmp(version_id,'ASCE 7-22')
            % use MPS spectrum from design value call diredctly
            rt_spectrum_LS.period = design_values_7_22.MPS.designSpectrum.periods';
            rt_spectrum_LS.sa_DE = design_values_7_22.MPS.designSpectrum.ordinates';
        else
            error('Not set up for this code version')
        end
        
        %% Calculate 2-point design values
        % for FREr
        [sdsfr, sd1fr] = fn_develop_2point_design_values(rt_spectrum_FR.sa, rt_spectrum_FR.period);

        % for DE
        if strcmp(version_id,'ASCE 7-28')
            % Pull values from MPS developed using hazard call
            [sds, sd1] = fn_develop_2point_design_values(rt_spectrum_LS.sa_DE, rt_spectrum_LS.period);
        elseif strcmp(version_id,'ASCE 7-22')
            % Pull values directly from design values call
            sds = design_values_7_22.sds;
            sd1 = design_values_7_22.sd1;
        else
            error('Not set up for this code version')
        end

        % Save code spectra to table
        sites.(['Sds_' site.site_class{v}])(s) = sds;
        sites.(['Sd1_' site.site_class{v}])(s) = sd1;
        sites.(['Sdsfr_' site.site_class{v}])(s) = sdsfr;
        sites.(['Sd1fr_' site.site_class{v}])(s) = sd1fr;

        %% Plot Spectra
        % FR - two point spectrum vector
        Ts = sd1fr/sdsfr;
        T0 = 0.2*(sd1fr/sdsfr);
        t_range = (Ts+0.1):0.1:5;
        FR_code_spectrum_T = [0 T0 Ts t_range];
        FR_code_spectrum = [0.4*sdsfr sdsfr sdsfr sd1fr./t_range];
        
        % LS - two point spectrum vector
        Ts = sd1/sds;
        T0 = 0.2*(sd1/sds);
        t_range = (Ts+0.1):0.1:5;
        LS_code_spectrum_T = [0 T0 Ts t_range];
        LS_code_spectrum = [0.4*sds sds sds sd1./t_range];
        
        % Plot and save
        hold on
        plot(rt_spectrum_LS.period, rt_spectrum_LS.sa_DE, '-k', 'DisplayName', 'DE: MPS')
        plot(LS_code_spectrum_T, LS_code_spectrum, '--k', 'DisplayName', 'DE: 2-point')
        plot(rt_spectrum_FR.period, rt_spectrum_FR.sa,'-','color',cmap(1,:),'DisplayName','FRE_R: MPS')
        plot(FR_code_spectrum_T, FR_code_spectrum,'--','color',cmap(1,:), 'DisplayName', 'FRE_R: 2-point')
        
        xlim([0 3])
        xlabel('Period (s)')
        ylabel('Sa (g)')
        legend('location','northeast')
        box on
        grid on
        title([version_id ': Site Class ' site.site_class{v}])
        saveas(gcf,[plot_dir filesep sites.name{s} ' - Site Class '  site.site_class{v} '.png'])
        close
    end
end

% Write table to remote save location
writetable(sites,[save_dir filesep 'frer_design_vals.csv'])







