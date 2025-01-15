clear all
close all
clc

%% Inputs
read_hazard_data = false;

T = [0.2 1 2 5];
T_lab = {'0p2' '1p0' '2p0' '5p0'};
maxdir_factor = 1./[1.2, 1.25, 1.256, 1.272]; % based on USGS conversion from Nico

vs30 = [760 265];
site_class = {'BC', 'D'};
% site_class = 'BC';

% vs30 = 265;
% site_class = 'D';

% return_periods = [14, 43, 72, 108, 144, 224, 475, 975, 2475, 5000];
% afe = 1 ./ return_periods;

regions = {'so cal' 'nor cal' 'PNW' 'WUS' 'CEUS'};
reg_tits = {'SoCal' 'NorCal' 'Pacific NW' 'Other WUS' 'CEUS'};

% Imports
import usgs.*

% Load model table
sites = readtable('sites.csv');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Filter to just plots of interest to look for outliers
% [~,w] = size(sites);
% filt = zeros(1,w);
% filt(1:4) = 1;
% sites = sites(:,contains(sites.Properties.VariableNames,'rtgm') | filt);
% [~,w] = size(sites);
% sites = sites(:,contains(sites.Properties.VariableNames,'beta_0p7') | filt(1:w));
% [~,w] = size(sites);
% sites = sites(:,contains(sites.Properties.VariableNames,'rt_5_PoE_10') | filt(1:w));
% 
% writetable(sites,'sites-file-trimmed.csv')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set web read options
options = weboptions('Timeout', 30);

%% Load hazard data for RTGM
rp = [375, 475];
beta = [0.5 0.7 0.9];
beta_label = {'0p5' '0p7' '0p9'};
PoE = [0.25 0.25];
poe_lab = {'25' '25'};
rt = [0.1, 0.15];
rt_lab = {'10', '15'};
t = 50;
nshm = [2023];
if read_hazard_data
    for s = 1:height(sites)
        s
        for m = 1:length(nshm)
            for v = 1:length(vs30)
                % Read hazard data for this site
                DATA = webread( sprintf('https://staging-earthquake.usgs.gov/ws/nshmp/conus-%i/dynamic/hazard/%0.2f/%0.2f/%i',...
                                nshm(m), ...
                                sites.lon(s), ...
                                sites.lat(s), ...
                                vs30(v)), ...
                                options);

                % Pull spectra data
                for t1 = 1:length(T)
                    for i = 1:length(DATA.response.hazardCurves)
                        spectral_label = DATA.response.hazardCurves(i).imt.value;          
                        if contains(spectral_label,'SA')
                            spectral_value = str2double(DATA.response.hazardCurves(i).imt.display(1:end-2));
                            if spectral_value == T(t1)
                                gms = DATA.response.hazardCurves(i).data(1).values.xs';
                                afes = DATA.response.hazardCurves(i).data(1).values.ys';
                                filt = afes > 0;
                                gms = gms(filt);
                                afes = afes(filt);
                                rps = 1./afes;
                            end
                        end
                    end
% 
%                     % Pull ground motion for this return period
%                     for r = 1:length(rp)
%                         sites.(['gm_' num2str(nshm(m)) '_' num2str(rp(r)) '_vs30_' num2str(vs30(v)) '_T_' T_lab{t1}])(s) = interp1(log(rps),gms,log(rp(r)));
%                     end

                    % Calculate and save risk targeted ground motions
                    for b = 1:length(beta)
                        for r1 = 1:length(rt)
                            % Calculate risk targeted ground motion
                            try
                            [rtgm] = fn_calc_rtgm(gms, afes, beta(b), PoE(r1), rt(r1), t);
                            catch
                                test = 5;
                            end

                            % Save rtgm to table
                            var_lab = [num2str(nshm(m)) '_beta_' beta_label{b} '_vs30_' num2str(vs30(v)) '_T_' T_lab{t1} '_rt_' rt_lab{r1} '_PoE_' poe_lab{r1}];
                            sites.(['rtgm_' var_lab])(s) = rtgm/maxdir_factor(t1);
                            sites.(['rp_' var_lab])(s) = 1/exp(interp1(gms,log(afes),rtgm));
                        end
                    end
                end
            end
        end
        
        % Write table
        writetable(sites,'sites-2023NSHM.csv')
    end
else
   sites = readtable('sites-2023NSHM.csv');
end


%% Plot Figures
plot_dir = 'C:\Users\dtc2\OneDrive - NIST\PUC\TS5\FRE Explore';

% Set matlab colors
cmap = [0 0.4470 0.7410;
        0.8500 0.3250 0.0980;
        0.9290 0.6940 0.1250;
        0.4940 0.1840 0.5560;
        0.4660 0.6740 0.1880;
        0.3010 0.7450 0.9330;
        0.6350 0.0780 0.1840];

    
for r1 = 1:length(rt)
    rt_var_lab = ['rt_' rt_lab{r1} '_PoE_' poe_lab{r1}];
    for b = 1:length(beta)
        save_dir = [plot_dir filesep rt_var_lab '_beta_' beta_label{b}];
        
        for t1 = 1:length(T)
            for v = 1:length(vs30)
                var_lab_save = ['vs30_' num2str(vs30(v)) '_T_' T_lab{t1}];
                var_lab = (['beta_' beta_label{b} '_vs30_' num2str(vs30(v)) '_T_' T_lab{t1} '_' rt_var_lab]);
                
                % Plot stacked histograms of RP
                min_rp = floor(min(sites.(['rp_2023_' var_lab]))/100)*100;
                max_rp = ceil(max(sites.(['rp_2023_' var_lab]))/100)*100;

                edges = min_rp:20:max_rp;
                centers = edges(2:end) - 10;

                counts = [];
                for reg = 1:length(regions)
                    filt = strcmp(sites.region,regions{reg});
                    counts(reg,:) = histcounts(sites.(['rp_2023_' var_lab])(filt), edges);    
                end

                hold on
                for reg = 1:length(regions)
                    bar(centers, sum(counts(reg:end,:),1),'FaceColor',cmap(reg,:), 'displayName', reg_tits{reg})
                end
                title(['Site Class ' site_class{v} ' - T = ' num2str(T(t1)) ' s'])
                box on
                xlabel('Return Period (years)')
                ylabel('Number of Sites')
                legend('location','northeast')

                mkdir([save_dir filesep 'ReturnPeriod_Hist' filesep '2023NSHM'])
                saveas(gcf, [save_dir filesep 'ReturnPeriod_Hist' filesep '2023NSHM' filesep var_lab_save '.png'])
                close


    %             % Scatter plot comparing NSHM models for specific risk targets
    %             figure
    %             title(['5% Risk in 50 years - Vs30 = ' num2str(vs30(v)) ' m/s'])
    %             hold on
    %             % scatter(sites.gm_2018_475, sites.gm_2023_475, 'filled')
    %             % regions = unique(models.region);
    %             regions = {'so cal' 'nor cal' 'PNW' 'WUS' 'CEUS'};
    %             reg_tits = {'SoCal' 'NorCal' 'Pacific NW' 'Other WUS' 'CEUS'};
    %             var_lab = (['beta_' beta_label{b} '_vs30_' num2str(vs30(v)) '_T_' T_lab{t1} '_rt_5_PoE_10']);
    %             for reg = 1:length(regions)
    %                 filt = strcmp(sites.region,regions{reg});
    %                 scatter(sites.(['rtgm_2018_' var_lab])(filt), sites.(['rtgm_2023_' var_lab])(filt),  'filled', 'displayName', reg_tits{reg})
    %             end
    %             lim = ceil(max([max(sites.(['rtgm_2018_' var_lab])), max(sites.(['rtgm_2023_' var_lab]))])*10)/10;
    %             plot([0 lim],[0 lim],'--k','handlevisibility','off')
    %             plot([0 lim],[0 lim]*0.8,':k','handlevisibility','off')
    %             plot([0 lim],[0 lim]*1.2,':k','displayName','+/- 20%')
    %             xlim([0 lim])
    %             ylim([0 lim])
    %             box on
    %             legend('location','northwest')
    %             xlabel(['2018 NSHM - Sa (T = ' num2str(T(t1)) ' s) (g)'])
    %             ylabel(['2023 NSHM - Sa (T = ' num2str(T(t1)) ' s) (g)'])
    % 
    %             saveas(gcf,['qq_rtgm_' var_lab '.png'])
    %             close
            end
        end
    end
end
% % Scatter plot comparing NSHM models with fixed return periods
% for r = 1:length(rp)
%     figure
%     title([num2str(rp(r)) ' Year Ground Motion'])
%     hold on
%     % scatter(sites.gm_2018_475, sites.gm_2023_475, 'filled')
%     % regions = unique(models.region);
%     regions = {'so cal' 'nor cal' 'PNW' 'WUS' 'CEUS'};
%     reg_tits = {'SoCal' 'NorCal' 'Pacific NW' 'Other WUS' 'CEUS'};
%     for reg = 1:length(regions)
%         filt = strcmp(sites.region,regions{reg});
%         scatter(sites.(['gm_2018_' num2str(rp(r))])(filt), sites.(['gm_2023_' num2str(rp(r))])(filt),  'filled', 'displayName', reg_tits{reg})
%     end
%     lim = ceil(max([max(sites.(['gm_2018_' num2str(rp(r))])), max(sites.(['gm_2023_' num2str(rp(r))]))])*10)/10;
%     plot([0 lim],[0 lim],'--k','handlevisibility','off')
%     plot([0 lim],[0 lim]*0.85,':k','handlevisibility','off')
%     plot([0 lim],[0 lim]*1.15,':k','displayName','+/- 15%')
%     xlim([0 lim])
%     ylim([0 lim])
%     box on
%     legend('location','northwest')
%     xlabel('2018 NSHM')
%     ylabel('2023 NSHM')
% 
%     saveas(gcf,['qq' num2str(rp(r)) '.png'])
%     close
% end





%%%%%%%%%%%% OLD LOGIC TO SAVE %%%%%%%%%%%%%%%%%%%%%
% 
% % Collect Design Info from USGS
% [design_values] = fn_call_USGS_design_value_API(1, 'asce7-22', models.lat(m), models.lng(m), 'II', site_class);
% 
% % Collect Hazard and Design Info
% [sa_spectra, ~] = fn_call_USGS_hazard_API('E2014', models.lat(m), models.lng(m), vs30, afe);
% 
% % Convert design values to RotD50
% coversion_factor = 1/1.3; % for T = 1 based on JWB code 
% coversion_factor = 1/1.2; % for T = 1 based on USGS conversion from Nico
% 
% % Calculate risk targeted 
% for s = 1:height(models)
%     sa = sa_vals{s,2:end};
%     filt = ~isnan(sa);
%     models.rp_DE(s,1) = 1/interp1(sa(filt),afe(filt),models.sdsrt(s));
%     for a = 1:length(PoE)
%         
%     end
% end



   




