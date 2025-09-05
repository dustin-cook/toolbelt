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
% site_set = 'nehrp_34';
site_set = 'frtc_testbed';

% USGS NSHM and ASCE version
version_id = 'ASCE 7-22';
nshm = 2018;
% version_id = 'ASCE 7-28';
% nshm = 2023;

% Remote save directory
save_dir = ['C:\Users\dtc2\OneDrive - NIST\PUC\TS5\FREr\' site_set '\' version_id];

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
sites_RP = table;

%% Load hazard data for RTGM
for s = 1:length(out_filt)
    site = out_filt(s);

    %Set Site values
    sites.id(s,1) = site.([site_set '_id']);
    sites.name{s,1} = site.name;
    sites.lat(s,1) = site.lat;
    sites.lon(s,1) = site.lon;
    sites.region{s,1} = site.region;
    sites.site_class{s,1} = site.site_class;

    sites_RP.id(s,1) = site.([site_set '_id']);
    sites_RP.name{s,1} = site.name;
    sites_RP.lat(s,1) = site.lat;
    sites_RP.lon(s,1) = site.lon;
    sites_RP.region{s,1} = site.region;
    sites_RP.site_class{s,1} = site.site_class;

    % Load pre-called values from the saved matlab data
    hazard_curve = site.hazard_curve.(['nshm_' num2str(nshm)]);
    design_values_7_22 = site.design_values; % These are strictly ASCE 7-22 design values

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
    sites.Sds(s,1) = sds;
    sites.Sd1(s,1) = sd1;
    sites.Sdsfr(s,1) = sdsfr;
    sites.Sd1fr(s,1) = sd1fr;

    % add MPS values for testbed cases
    if strcmp(site_set,'frtc_testbed')
        for t = 1:length(rt_spectrum_FR.period)
            t_str = strrep(sprintf('%4.2f',rt_spectrum_FR.period(t)),'.','p');
            sites.(['frer_mps_' t_str])(s,1) = rt_spectrum_FR.sa(t);
        end
    end

    %% Quantify Return Periods of Spectral Values at specific periods
%     if strcmp(site_set,'usgs_160')
%         if strcmp(version_id,'ASCE 7-28')
%             sites_RP.frer_0p2(s,1) = rt_spectrum_FR.return_period(rt_spectrum_FR.period == 0.2);
%             sites_RP.frer_1p0(s,1) = rt_spectrum_FR.return_period(rt_spectrum_FR.period == 1.0);
% 
%             sites_RP.mcer_0p2(s,1) = rt_spectrum_LS.return_period(rt_spectrum_LS.period == 0.2);
%             sites_RP.mcer_1p0(s,1) = rt_spectrum_LS.return_period(rt_spectrum_LS.period == 1.0);
% 
%             sites_RP.DE_0p2(s,1) = rt_spectrum_LS.return_period_67(rt_spectrum_LS.period == 0.2);
%             sites_RP.DE_1p0(s,1) = rt_spectrum_LS.return_period_67(rt_spectrum_LS.period == 1.0);
%         else
%             % NOTE: currenly not able to pull return periods for ASCE 7-22
%             % MCE and design values. I would need to interpolate from the
%             % hazard curve using the final mixed RT and derteministic numbers
%             % (be sure to check the risk target is 1% and not 1.5%)
%             error('Update tool to calculate RP for ASCE 7-22')
%         end
%     end
    
    %% Plot Spectra for NEHRP Sites
    if strcmp(site_set,'nehrp_34') || strcmp(site_set,'frtc_testbed')
        plot_dir = [save_dir '\Spectra_Plots - loglog'];
        if ~exist(plot_dir,'dir')
            mkdir(plot_dir)
        end

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
        loglog(rt_spectrum_LS.period, rt_spectrum_LS.sa_DE, '-k', 'DisplayName', 'DE: MPS')
        hold on
        LS_code_spectrum_T(1) = 0.01;
        FR_code_spectrum_T(1) = 0.01;
        loglog(LS_code_spectrum_T, LS_code_spectrum, '--k', 'DisplayName', 'DE: 2-point')
        loglog(rt_spectrum_FR.period, rt_spectrum_FR.sa,'-','color',cmap(1,:),'DisplayName','FRE_R: MPS')
        loglog(FR_code_spectrum_T, FR_code_spectrum,'--','color',cmap(1,:), 'DisplayName', 'FRE_R: 2-point')

        xlim([0 5])
        xlabel('Period (s)')
        ylabel('Sa (g)')
        legend('location','southwest')
        box on
        grid on
        title([version_id ': Site Class ' site.site_class])
        saveas(gcf,[plot_dir filesep sites.name{s} ' - Site Class '  site.site_class '.png'])
        close
    end
end

%% Write table to remote save location
writetable(sites,[save_dir filesep 'frer_design_vals.csv'])
if strcmp(site_set,'usgs_160')
    writetable(sites_RP,[save_dir filesep 'frer_return_periods.csv'])
end

%% Plot Stacked Return Period Plots
% T_lab = {'0p2', '1p0'};
% HL_lab = {'frer', 'mcer', 'DE'};
% HL_tit = {'FRE_R', 'MCE_R', 'DE'};
% 
% if strcmp(site_set,'usgs_160')
%     plot_dir = [save_dir '\Return Period Plots'];
%     if ~exist(plot_dir,'dir')
%         mkdir(plot_dir)
%     end
%     
%     site_class = unique(sites_RP.site_class);
%     regions = unique(sites_RP.region);
%     for s = 1:length(site_class)
%         filt_SC = strcmp(sites_RP.site_class,site_class{s});
%         for t = 1:length(T_lab)
%             for gm = 1:length(HL_lab)
%                 data_rp = sites_RP.([HL_lab{gm} '_' T_lab{t}])(filt_SC);
%                 min_rp = floor(min(data_rp)/100)*100;
%                 max_rp = ceil(max(data_rp)/100)*100;
% 
%                 bin_width = (max_rp - min_rp)/20;
%                 edges = min_rp:bin_width:max_rp;
%                 centers = edges(2:end) - bin_width/2;
% 
%                 counts = [];
%                 for reg = 1:length(regions)
%                     filt_reg = strcmp(sites_RP.region(filt_SC),regions{reg});
%                     counts(reg,:) = histcounts(data_rp(filt_reg), edges);    
%                 end
% 
%                 hold on
%                 for reg = 1:length(regions)
%                     bar(centers, sum(counts(reg:end,:),1),'FaceColor',cmap(reg,:), 'displayName', regions{reg})
%                 end
%                 title([HL_tit{gm} ' Return Period - T = ' strrep(T_lab{t},'p','.') 's - Site Class ' site_class{s}])
%                 box on
%                 xlabel('Return Period (years)')
%                 ylabel('Number of Sites')
%                 legend('location','northeast')
%                 saveas(gcf, [plot_dir filesep 'rp_' HL_lab{gm} '_t' T_lab{t} '_siteclass_' site_class{s} '.png'])
%         %         set(gcf, 'Color', 'None'
%         %         set(gca, 'color', 'none');
%                 close
%             end
%         end
%     end
% end

%% Plot scatter plots comparing FREr with DE
if strcmp(site_set,'usgs_160')
    plot_dir = [save_dir '\Scatter Plots'];
    if ~exist(plot_dir,'dir')
        mkdir(plot_dir)
    end
            
    site_class = unique(sites_RP.site_class);
    regions = unique(sites.region);
    for s = 1:length(site_class)
        filt_SC = strcmp(sites.site_class,site_class{s});
        
        % Short Period Spectral Acceleration
        figure
        hold on
        for r = 1:length(regions)
            filt_reg = strcmp(sites.region,regions{r});
            scatter(sites.Sds(filt_reg & filt_SC), sites.Sdsfr(filt_reg & filt_SC), 'filled', 'displayName', regions{r})
        end
        lim = ceil(max([max(sites.Sds(filt_SC)), max(sites.Sds(filt_SC))])*10)/10;
        plot([0 lim],[0 lim],'--k','handlevisibility','off')
        plot([0 lim],[0 lim]*0.8,':k','handlevisibility','off')
        plot([0 lim],[0 lim]*1.2,':k','displayName','+/- 20%')
        xlim([0 lim])
        ylim([0 lim])
        box on
        title(['Site Class ' site_class{s}])
        legend('location','northwest')
        xlabel('S_{DS} (g)')
        ylabel('S_{DSfr} (g)')
        saveas(gcf, [plot_dir filesep 'Sds_v_Sdsfr_' site_class{s} '.png'])
        close
        
        % 1 second Spectral Acceleration
        figure
        hold on
        for r = 1:length(regions)
            filt_reg = strcmp(sites.region,regions{r});
            scatter(sites.Sd1(filt_reg & filt_SC), sites.Sd1fr(filt_reg & filt_SC), 'filled', 'displayName', regions{r})
        end
        lim = ceil(max([max(sites.Sd1(filt_SC)), max(sites.Sd1(filt_SC))])*10)/10;
        plot([0 lim],[0 lim],'--k','handlevisibility','off')
        plot([0 lim],[0 lim]*0.8,':k','handlevisibility','off')
        plot([0 lim],[0 lim]*1.2,':k','displayName','+/- 20%')
        xlim([0 lim])
        ylim([0 lim])
        box on
        title(['Site Class ' site_class{s}])
        legend('location','northwest')
        xlabel('S_{D1} (g)')
        ylabel('S_{D1fr} (g)')
        saveas(gcf, [plot_dir filesep 'Sd1_v_Sd1fr_' site_class{s} '.png'])
        close
    end
end





