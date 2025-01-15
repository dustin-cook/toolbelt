clear all
close all
clc

%% Inputs
read_hazard_data = false;

T = [0.2 1 2 5];
T_lab = {'0p2' '1p0' '2p0' '5p0'};
maxdir_factor = 1./[1.2, 1.25, 1.256, 1.272]; % based on USGS conversion from Nico

maxdir_factor_s = 1.2;
maxdir_factor_1 = 1.25;


vs30 = [530 760 265];
site_class = {'C' 'BC', 'D'};
% vs30 = [760 265];
% site_class = {'BC', 'D'};

% vs30 = [530];
% site_class = {'C'};

regions = {'so cal' 'nor cal' 'PNW' 'WUS' 'CEUS'};
reg_tits = {'SoCal' 'NorCal' 'Pacific NW' 'Other WUS' 'CEUS'};

% vs30 = 265;
% site_class = 'D';

% return_periods = [14, 43, 72, 108, 144, 224, 475, 975, 2475, 5000];
% afe = 1 ./ return_periods;

% Imports
import usgs.*

% Load model table
% sites = readtable('sites_fre_geomean_2024_06_13.csv');
% sites = sites(151,:);
% sites = readtable('sites_fre_2024_04_15.csv');
% sites = readtable('sites_fre.csv');
% sites = readtable('sites4_EERI.csv');
% sites = readtable('sites_21_cities.csv');
% sites = readtable('sites_34cities.csv');

% sites = readtable('sites_fre_geomean_2024_08_20.csv'); % This is just the 34 cities - both of these I need to rerun for beta = 0.6...
sites = readtable('sites_fre_geomean_2024_08_13.csv'); % This is all of them

% 
% nehrp = readtable('nehrp34_sites.csv');
% filt_rows = ismember(sites.lat,nehrp.lat) & ismember(sites.lon,nehrp.lng);
% sites = sites(filt_rows,:);
% writetable(sites,'sites_fre_geomean_34cities.csv')

% var_names = {'name','lon','lat','region',...
%             'Sa_design_vs30_760_T_0p2',...   
%             'rtgm_2018_beta_0p7_vs30_760_T_0p2_rt_10_PoE_25'};
% filt_cols = ismember(sites.Properties.VariableNames,var_names);
% sites_filt_table = sites(filt_rows,filt_cols);
% sites_filt_table.Properties.VariableNames{5} = 'Sa_DE';
% sites_filt_table.Properties.VariableNames{6} = 'Sa_FRER';
% sites_filt_table.Sa_FRER = min(sites_filt_table.Sa_FRER,sites_filt_table.Sa_DE);
% sites_filt_table.R_RFR_1p0 = sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_1p5 = 1.5 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_2p0 = 2 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_3p0 = 3 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% writetable(sites_filt_table,'sites34_vs30_760_T_0p2.csv')
% 
% var_names = {'name','lon','lat','region',...  
%             'Sa_design_vs30_760_T_1p0',...
%             'rtgm_2018_beta_0p7_vs30_760_T_1p0_rt_10_PoE_25'};
% filt_cols = ismember(sites.Properties.VariableNames,var_names);
% sites_filt_table = sites(filt_rows,filt_cols);
% sites_filt_table.Properties.VariableNames{5} = 'Sa_DE';
% sites_filt_table.Properties.VariableNames{6} = 'Sa_FRER';
% sites_filt_table.Sa_FRER = min(sites_filt_table.Sa_FRER,sites_filt_table.Sa_DE);
% sites_filt_table.R_RFR_1p0 = sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_1p5 = 1.5 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_2p0 = 2 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_3p0 = 3 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% writetable(sites_filt_table,'sites34_vs30_760_T_1p0.csv')
% 
% var_names = {'name','lon','lat','region',...
%             'Sa_design_vs30_265_T_0p2',...   
%             'rtgm_2018_beta_0p7_vs30_265_T_0p2_rt_10_PoE_25'};
% filt_cols = ismember(sites.Properties.VariableNames,var_names);
% sites_filt_table = sites(filt_rows,filt_cols);
% sites_filt_table.Properties.VariableNames{5} = 'Sa_DE';
% sites_filt_table.Properties.VariableNames{6} = 'Sa_FRER';
% sites_filt_table.Sa_FRER = min(sites_filt_table.Sa_FRER,sites_filt_table.Sa_DE);
% sites_filt_table.R_RFR_1p0 = sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_1p5 = 1.5 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_2p0 = 2 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_3p0 = 3 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% writetable(sites_filt_table,'sites34_vs30_265_T_0p2.csv')
% 
% var_names = {'name','lon','lat','region',...  
%             'Sa_design_vs30_265_T_1p0',...
%             'rtgm_2018_beta_0p7_vs30_265_T_1p0_rt_10_PoE_25'};
% filt_cols = ismember(sites.Properties.VariableNames,var_names);
% sites_filt_table = sites(filt_rows,filt_cols);
% sites_filt_table.Properties.VariableNames{5} = 'Sa_DE';
% sites_filt_table.Properties.VariableNames{6} = 'Sa_FRER';
% sites_filt_table.Sa_FRER = min(sites_filt_table.Sa_FRER,sites_filt_table.Sa_DE);
% sites_filt_table.R_RFR_1p0 = sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_1p5 = 1.5 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_2p0 = 2 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% sites_filt_table.R_RFR_3p0 = 3 * sites_filt_table.Sa_FRER ./ sites_filt_table.Sa_DE;
% writetable(sites_filt_table,'sites34_vs30_265_T_1p0.csv')

% tmp = 5;
% mps_test = sites(:,1:4);

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
options = weboptions('Timeout', 60);

%% Load hazard data for RTGM
% beta = [0.5 0.6 0.7 1.0];
beta = [0.6];
beta_label = {'0p6'};
% beta_label = {'0p5' '0p6' '0p7' '1p0'};
% PoE = [0.1 0.25 0.1 0.25 0.1 0.25 0.1 0.25];
% poe_lab = {'10' '25' '10' '25' '10' '25' '10' '25'};
PoE = [0.25];
poe_lab = {'25'};
% rt = [0.1 0.1 0.15 0.15 0.2 0.2 0.25 0.25];
% rt_lab = {'10' '10' '15' '15' '20' '20' '25' '25'};
% rt = [0.1, 0.15, 0.2, 0.25, 0.3, 0.5];
% rt_lab = {'10', '15', '20', '25', '30', '50'};
rt = [0.1];
rt_lab = {'10'};
t = 50;
nshm = [2018];
if read_hazard_data
    for s = 6:height(sites)
        s
        for m = 1:length(nshm)
            for v = 1:length(vs30)
                
                % Pull design values for this site
                [design_values, MPS] = fn_call_USGS_design_value_API_v2(1, 'asce7-22', sites.lat(s), sites.lon(s), 'II', site_class{v});
                    
                % Read hazard data for this site
%                 'https://staging-earthquake.usgs.gov/ws/nshmp/conus-2018/dynamic/hazard/-97.35/37.7/760'
                DATA = webread( sprintf('https://earthquake.usgs.gov/ws/nshmp/conus-%i/dynamic/hazard/%0.2f/%0.2f/%i',...
                                nshm(m), ...
                                sites.lon(s), ...
                                sites.lat(s), ...
                                vs30(v)), ...
                                options);

                % Save hazard data to compare MPS at T = 1s
                mps_test.(['s1_' site_class{v}])(s) = design_values.s1;
                mps_test.(['sd1_' site_class{v}])(s) = design_values.sd1;
                mps_test.(['mps_' site_class{v}])(s) = design_values.mps;
                mps_test.(['mps_rt_' site_class{v}])(s) = design_values.mps_rt;
                mps_test.(['mps_84_' site_class{v}])(s) = design_values.mps_84;
                
%                 MPS.designSpectrum.ordinates(MPS.designSpectrum.periods==1)
                
                
                % Save meta data
                vs30_lab = ['vs30_' num2str(vs30(v))];
                sites.(['SDC_' vs30_lab]){s} = design_values.SDC;
                sites.(['s1_' vs30_lab])(s) = design_values.s1;
                sites.(['sds_' vs30_lab])(s) = design_values.sds;
                sites.(['sd1_' vs30_lab])(s) = design_values.sd1;
                    
                T = MPS.designSpectrum.periods;
                
                % Pull spectra data
                frer_spectra = [];
                for t1 = 1:length(T)
                    found_vals = 0;
                    gms = [];
                    afes = [];
                    rps = []; 
                    
%                     % Pull hazard curve
%                     for i = 1:length(DATA.response.hazardCurves)
%                         spectral_label = DATA.response.hazardCurves(i).imt.value;         
%                         if contains(spectral_label,'SA')
%                             spectral_value = str2double(DATA.response.hazardCurves(i).imt.display(1:end-2));
%                             if spectral_value == T(t1)
%                                 found_vals = 1;
%                                 
%                                 gms = DATA.response.hazardCurves(i).data(1).values.xs';
%                                 afes = DATA.response.hazardCurves(i).data(1).values.ys';
%                                 filt = afes > 0;
%                                 gms = gms(filt);
%                                 afes = afes(filt);
%                                 rps = 1./afes;
%                             end
%                         end
%                     end
                    
                    % Pull hazard curve for this period
                    for i = 1:length(DATA.response.hazardCurves)
                        spectral_label = DATA.response.hazardCurves(i).imt.value;   
                        spectral_value = str2double(DATA.response.hazardCurves(i).imt.display(1:end-2));
          
                        if (strcmp(spectral_label,'PGA') && T(t1) == 0) || (contains(spectral_label,'SA') && spectral_value == T(t1))
                            found_vals = 1;
                            gms = DATA.response.hazardCurves(i).data(1).values.xs';
                            afes = DATA.response.hazardCurves(i).data(1).values.ys';
                            filt = afes > 0;
                            gms = gms(filt);
                            afes = afes(filt);
                            rps = 1./afes;

                        end
                    end
                    
                    if found_vals == 0
                        herestheproblem = 1;
                    end
                    
                    % Pull design values - I would like to use the MPS
                    % here, but something strange is going on with their
                    % site class correction compared to the two point 
                    % design values
%                     de_lab = ['vs30_' num2str(vs30(v)) '_T_' T_lab{t1}];
                    
%                     sites.(['DE_' de_lab])(s) = MPS.designSpectrum.ordinates(MPS.designSpectrum.periods == T(t1));
%                     sites.(['Sa_rt_' de_lab])(s) = MPS.riskTargetedSpectrum.ordinates(MPS.designSpectrum.periods == T(t1));
                    
%                     Ts_design = design_values.sd1 / design_values.sds;
%                     if T(t1) < Ts_design
%                         sites.(['Sa_design_' de_lab])(s) = design_values.sds;
%                     else
%                         sites.(['Sa_design_' de_lab])(s) = design_values.sd1/T(t1);
%                     end
%                     sites.(['RP_design_' de_lab])(s) = 1/exp(interp1(gms,log(afes),maxdir_factor(t1) * sites.(['Sa_design_' de_lab])(s)));
%                     
%                     Ts_rt = design_values.sd1rt / design_values.sdsrt;
%                     if T(t1) < Ts_rt
%                         sites.(['Sa_rt_' de_lab])(s) = design_values.sdsrt;
%                     else
%                         sites.(['Sa_rt_' de_lab])(s) = design_values.sd1rt/T(t1);
%                     end
%                     
%                     test_rp = 1/interp1(gms,afes,sites.(['Sa_rt_' de_lab])(s));
%                     if test_rp > 10000
%                         pause = 1;
%                     end
                    
%                     sites.(['RP_rt_' de_lab])(s) = 1/exp(interp1(gms,log(afes),maxdir_factor(t1) * sites.(['Sa_rt_' de_lab])(s)));
                    
                    % Calculate and save risk targeted ground motions
                    for b = 1:length(beta)
                        for r1 = 1:length(rt)
                            % Calculate risk targeted ground motion
                            try
                            [rtgm] = fn_calc_rtgm(gms, afes, beta(b), PoE(1), rt(r1), t);
                            catch
                                test = 5;
                            end
                            
                            frer_spectra(t1,b) = rtgm;
%                             figure
%                             hold on
%                             yyaxis left
%                             plot(gms,afes,'linewidth',2)
%                             xlim([0 2.5])
%                             xlabel('Spectral Acceleration (g)')
%                             ylim([1e-4, 1e0])
%                             ax1 = gca;
%                             ax1.YScale = 'log';
%                             ylabel('1/Return Period (1/years)')
%                             yyaxis right
%                             plot(gms,logncdf(gms,log(rtgm),beta(b)),'linewidth',2);
%                             ylim([0, 0.8])
%                             ylabel('P[Performance > Target]')
%                             box on
%                             set(gcf,'position',[0,0,500,350])
%                             close

                            % Save rtgm to table
%                             var_lab = [num2str(nshm(m)) '_beta_' beta_label{b} '_vs30_' num2str(vs30(v)) '_T_' T_lab{t1} '_rt_' rt_lab{r1} '_PoE_' poe_lab{1}];
%                             sites.(['rtgm_' var_lab])(s) = rtgm / maxdir_factor(t1);
%                             sites.(['FREr_' var_lab])(s) = rtgm;
%                             sites.(['rp_' var_lab])(s) = 1/exp(interp1(gms,log(afes),rtgm));
                        end
                    end
                end
                
                for b = 1:length(beta)
                    % Find two point FREr design spectrum 
                    T_range_s = MPS.designSpectrum.periods >= 0.2 & MPS.designSpectrum.periods <= 5.0;
                    sdsfr = 0.9*max(frer_spectra(T_range_s, b));
                    T_range_1 = MPS.designSpectrum.periods >= 1.0 & MPS.designSpectrum.periods <= 2.0;
                    sd1fr = max(frer_spectra(MPS.designSpectrum.periods == 1.0, b),...
                                0.9*max(MPS.designSpectrum.periods(T_range_1).*frer_spectra(T_range_1, b)));
                    Ts = sd1fr/sdsfr;
                    T0 = 0.2*(sd1fr/sdsfr);
                    t_range = (Ts+0.1):0.1:5;
                    code_spectrum_T = [0 T0 Ts t_range];
                    code_spectrum = [0.4*sdsfr sdsfr sdsfr sd1fr./t_range];

                    sites.(['Sdsfr_' vs30_lab '_geomean'])(s) = sdsfr;
                    sites.(['Sd1fr_' vs30_lab '_geomean'])(s) = sd1fr;
                    sites.(['Sdsfr_' vs30_lab '_maxdir'])(s) = sdsfr*maxdir_factor_s;
                    sites.(['Sd1fr_' vs30_lab '_maxdir'])(s) = sd1fr*maxdir_factor_1;
                    
                end
            end
        end
        
        % Write table
        writetable(sites,'sites_fre_geomean_2024_08_20.csv')
%         writetable(sites,'sites_fre_Dmax_2024_07_25.csv')
%         writetable(mps_test,'mps_test.csv')
        
%         histogram(mps_test.(['mps_' site_class{1}])*2/3 ./ mps_test.(['sd1_' site_class{1}]))
    end
else
%     sites = readtable('sites_fre_geomean_2024_06_13.csv');
end

%% Scatter plots to compare various beta values
% plot_dir = 'C:\Users\dtc2\OneDrive - NIST\PUC\TS5\Beta Study\FREr-scatter';
% mkdir(plot_dir)
% 
% 
% % Beta of 0.5 to 0.7
% for v = 1:length(vs30)
%     % Sdsfr
%     hold on
%     scatter(sites.(['Sdsfr_vs30_' num2str(vs30(v)) '_beta_0p6']), sites.(['Sdsfr_vs30_' num2str(vs30(v)) '_beta_0p5']),'filled','MarkerFaceAlpha', 0.5,'displayname','\beta = 0.5')
%     scatter(sites.(['Sdsfr_vs30_' num2str(vs30(v)) '_beta_0p6']), sites.(['Sdsfr_vs30_' num2str(vs30(v)) '_beta_0p7']),'filled','MarkerFaceAlpha', 0.5,'displayname','\beta = 0.7')
%     plot([0 2],[0 2],'-k','handlevisibility','off')
%     plot([0 2],[0 2.2],'--k','displayname','+/- 10%')
%     plot([0 2],[0 1.8],'--k','handlevisibility','off')
%     box on
%     grid on
%     max_val = ceil(max(sites.(['Sdsfr_vs30_' num2str(vs30(v)) '_beta_0p6']))*10)/10;
%     xlim([0 max_val])
%     ylim([0 max_val])
%     xlabel('FRE_R | \beta = 0.6')
%     ylabel('FRE_R | modified \beta')
%     legend('location','northwest')
%     title(['S_{dsfr} [Site Class ' site_class{v} ']'])
%     saveas(gcf,[plot_dir filesep 'Sdsfr SiteClass_' site_class{v} '.png'])
%     close
%     
%     % Sdsfr
%     hold on
%     scatter(sites.(['Sd1fr_vs30_' num2str(vs30(v)) '_beta_0p6']), sites.(['Sd1fr_vs30_' num2str(vs30(v)) '_beta_0p5']),'filled','MarkerFaceAlpha', 0.5,'displayname','\beta = 0.5')
%     scatter(sites.(['Sd1fr_vs30_' num2str(vs30(v)) '_beta_0p6']), sites.(['Sd1fr_vs30_' num2str(vs30(v)) '_beta_0p7']),'filled','MarkerFaceAlpha', 0.5,'displayname','\beta = 0.7')
%     plot([0 2],[0 2],'-k','handlevisibility','off')
%     plot([0 2],[0 2.2],'--k','displayname','+/- 10%')
%     plot([0 2],[0 1.8],'--k','handlevisibility','off')
%     box on
%     grid on
%     max_val = ceil(max(sites.(['Sd1fr_vs30_' num2str(vs30(v)) '_beta_0p6']))*10)/10;
%     xlim([0 max_val])
%     ylim([0 max_val])
%     xlabel('FRE_R | \beta = 0.6')
%     ylabel('FRE_R | modified \beta')
%     legend('location','northwest')
%     title(['S_{d1fr} [Site Class ' site_class{v} ']'])
%     saveas(gcf,[plot_dir filesep 'Sd1fr SiteClass_' site_class{v} '.png'])
%     close
% end
% 
% % Beta of 0.6 to 1.0
% for v = 1:length(vs30)
%     % Sdsfr
%     hold on
%     scatter(sites.(['Sdsfr_vs30_' num2str(vs30(v)) '_beta_0p6']), sites.(['Sdsfr_vs30_' num2str(vs30(v)) '_beta_1p0']),'filled','MarkerFaceAlpha', 0.5,'handlevisibility','off')
%     plot([0 2],[0 2],'-k','handlevisibility','off')
%     plot([0 2],[0 2.2],'--k','displayname','+/- 10%')
%     plot([0 2],[0 1.8],'--k','handlevisibility','off')
%     box on
%     grid on
%     max_val = ceil(max(sites.(['Sdsfr_vs30_' num2str(vs30(v)) '_beta_0p6']))*10)/10;
%     xlim([0 max_val])
%     ylim([0 max_val])
%     xlabel('FRE_R | \beta = 0.6')
%     ylabel('FRE_R | \beta = 1.0')
%     legend('location','northwest')
%     title(['S_{dsfr} [Site Class ' site_class{v} ']'])
%     saveas(gcf,[plot_dir filesep 'b1 - Sdsfr SiteClass_' site_class{v} '.png'])
%     close
%     
%     % Sdsfr
%     hold on
%     scatter(sites.(['Sd1fr_vs30_' num2str(vs30(v)) '_beta_0p6']), sites.(['Sd1fr_vs30_' num2str(vs30(v)) '_beta_1p0']),'filled','MarkerFaceAlpha', 0.5,'handlevisibility','off')
%     plot([0 2],[0 2],'-k','handlevisibility','off')
%     plot([0 2],[0 2.2],'--k','displayname','+/- 10%')
%     plot([0 2],[0 1.8],'--k','handlevisibility','off')
%     box on
%     grid on
%     max_val = ceil(max(sites.(['Sd1fr_vs30_' num2str(vs30(v)) '_beta_0p6']))*10)/10;
%     xlim([0 max_val])
%     ylim([0 max_val])
%     xlabel('FRE_R | \beta = 0.6')
%     ylabel('FRE_R | \beta = 1.0')
%     legend('location','northwest')
%     title(['S_{d1fr} [Site Class ' site_class{v} ']'])
%     saveas(gcf,[plot_dir filesep 'b1 - Sd1fr SiteClass_' site_class{v} '.png'])
%     close
% end


%% Post process to find Dmin and Dmax sites
% sites_Dmax = table;
% sites_Dmin = table;
% for v = 1:length(vs30)
%     for r1 = 1:length(rt)
%         tmp_tab = table;
%         vs30_lab = ['vs30_' num2str(vs30(v))];
% 
%         % new Vs30 Specific table
%         tmp_tab.name = sites.name;
%         tmp_tab.lon = sites.lon;
%         tmp_tab.lat = sites.lat;
%         tmp_tab.region = sites.region;
%         tmp_tab.vs30 = ones(height(tmp_tab),1)*vs30(v);
%         tmp_tab.SDC = sites.(['SDC_' vs30_lab]);
%         tmp_tab.risk_target = ones(height(tmp_tab),1)*rt(r1);
% 
%         % Dmax ratio
%         tmp_tab.ratio_Dmax = sites.(['s1_' vs30_lab])/0.75;
% 
%         % Dmin ratio
%         dmin_sds = sites.(['sds_' vs30_lab])/0.5;
%         dmin_sd1 = sites.(['sd1_' vs30_lab])/0.2;
%         tmp_tab.ratio_Dmin = max(dmin_sds,dmin_sd1);
% 
%         % FREr Ratio
%         tmp_tab.frer_ratio_sds = sites.(['Sdsfr_' vs30_lab '_RT_' rt_lab{r1}]) ./ sites.(['sds_' vs30_lab]);
%         tmp_tab.frer_ratio_sd1 = sites.(['Sd1fr_' vs30_lab '_RT_' rt_lab{r1}]) ./ sites.(['sd1_' vs30_lab]);
% 
%         % populate rest of values
%         tmp_tab.Sdsfr = sites.(['Sdsfr_' vs30_lab '_RT_' rt_lab{r1}]);
%         tmp_tab.Sd1fr = sites.(['Sd1fr_' vs30_lab '_RT_' rt_lab{r1}]);
%         tmp_tab.s1 = sites.(['s1_' vs30_lab]);
%         tmp_tab.sds = sites.(['sds_' vs30_lab]);
%         tmp_tab.sd1 = sites.(['sd1_' vs30_lab]);
% 
%         % Filter to Dmax sites
%         dmax_filt = tmp_tab.ratio_Dmax >= 0.9 & tmp_tab.ratio_Dmax <= 1.1;
%         sites_Dmax = [sites_Dmax; tmp_tab(dmax_filt,:)];
% 
%         % Filter to Dmin sites
%         dmin_filt = tmp_tab.ratio_Dmin >= 0.9 & tmp_tab.ratio_Dmin <= 1.1;
%         sites_Dmin = [sites_Dmin; tmp_tab(dmin_filt,:)];
%     end
% end
% 
% % Save tables
% writetable(sites_Dmax,'fre_Dmax_2024_07_25.csv')
% writetable(sites_Dmin,'fre_Dmin_2024_07_25.csv')

%% Plot Figures
% Set plot dir
% plot_dir = 'C:\Users\dtc2\OneDrive - NIST\Conferences\2024 EERI AM';
plot_dir = 'C:\Users\dtc2\OneDrive - NIST\PUC\TS5\FRE Explore\2024_05_24_geomean';
summary_table = table;


% Set matlab colors
cmap = [0 0.4470 0.7410;
        0.8500 0.3250 0.0980;
        0.9290 0.6940 0.1250;
        0.4940 0.1840 0.5560;
        0.4660 0.6740 0.1880;
        0.3010 0.7450 0.9330;
        0.6350 0.0780 0.1840];

% New Scatter plots
for v = 1:length(vs30)
    % Plot sds compare
    figure
    title(['Site Class ' site_class{v}])
    hold on
    for reg = 1:length(regions)
        filt = strcmp(sites.region,regions{reg});
        scatter(sites.(['sds_vs30_' num2str(vs30(v))])(filt), sites.(['Sdsfr_vs30_' num2str(vs30(v)) '_beta_0p6'])(filt),  'filled', 'displayName', reg_tits{reg})
    end
    lim = ceil(max([max(sites.(['sds_vs30_' num2str(vs30(v))])), max(sites.(['Sdsfr_vs30_' num2str(vs30(v)) '_beta_0p6']))])*10)/10;
    plot([0 lim],[0 lim],'--k','handlevisibility','off')
%             plot([0 lim],[0 lim]*0.8,':k','handlevisibility','off')
%             plot([0 lim],[0 lim]*1.2,':k','displayName','+/- 20%')
    xlim([0 lim])
    ylim([0 lim])
    box on
    legend('location','northwest')
    xlabel('S_{DS} (g)')
    ylabel('S_{DSfr} (g)')
    saveas(gcf, [plot_dir filesep 'Sds_vs30_' num2str(vs30(v)) '.png'])
%     set(gcf, 'Color', 'None')
%     set(gca, 'color', 'none');
    close
    
    % Plot sd1 compare
    figure
    title(['Site Class ' site_class{v}])
    hold on
    for reg = 1:length(regions)
        filt = strcmp(sites.region,regions{reg});
        scatter(sites.(['sd1_vs30_' num2str(vs30(v))])(filt), sites.(['Sd1fr_vs30_' num2str(vs30(v)) '_beta_0p6'])(filt),  'filled', 'displayName', reg_tits{reg})
    end
    lim = ceil(max([max(sites.(['sd1_vs30_' num2str(vs30(v))])), max(sites.(['Sd1fr_vs30_' num2str(vs30(v)) '_beta_0p6']))])*10)/10;
    plot([0 lim],[0 lim],'--k','handlevisibility','off')
%             plot([0 lim],[0 lim]*0.8,':k','handlevisibility','off')
%             plot([0 lim],[0 lim]*1.2,':k','displayName','+/- 20%')
    xlim([0 lim])
    ylim([0 lim])
    box on
    legend('location','northwest')
    xlabel('S_{D1} (g)')
    ylabel('S_{D1fr} (g)')
    saveas(gcf, [plot_dir filesep 'Sd1_vs30_' num2str(vs30(v)) '.png'])
%     set(gcf, 'Color', 'None')
%     set(gca, 'color', 'none');
    close
end

    
%% OLD PLOTTERS
% Scatter plot comparing NSHM models for specific risk targets
c = 0;
for r1 = 1:length(rt)
    rt_var_lab = ['rt_' rt_lab{r1} '_PoE_' poe_lab{r1}];
    for b = 1:length(beta)
        save_dir = [plot_dir filesep rt_var_lab '_beta_' beta_label{b}];

        fre_table = sites(:,1:4);
        rp_table = sites(:,1:4);
        cr_de_table = sites(:,1:4);
        cr_rt_table = sites(:,1:4);

        for v = 1:length(vs30)
            for t1 = 1:length(T)
                tab_lab = (['vs30_' num2str(vs30(v)) '_T_' T_lab{t1}]);
                var_lab_save = ['vs30_' num2str(vs30(v)) '_T_' T_lab{t1}];
                var_lab = (['beta_' beta_label{b} '_vs30_' num2str(vs30(v)) '_T_' T_lab{t1} '_' rt_var_lab]);

                % Add data to table
                fre_table.(tab_lab) = sites.(['rtgm_2018_' var_lab]);
                rp_table.(tab_lab) = sites.(['rp_2018_' var_lab]);
                cr_de_table.(tab_lab) = sites.(['rtgm_2018_' var_lab]) ./ sites.(['Sa_design_' tab_lab]) ;
                cr_rt_table.(tab_lab) = sites.(['rtgm_2018_' var_lab]) ./ sites.(['Sa_rt_' tab_lab]);

                % Cap values at DE
%                 data_tmp = min(sites.(['rp_2018_' var_lab]),sites.(['RP_design_' tab_lab]),'includenan');
                
                % DO NOT CAP RETURN PERIODS...
                data_tmp = sites.(['rp_2018_' var_lab]);
                
                % Plot stacked histograms of RP
%                 min_rp = floor(min(sites.(['rp_2018_' var_lab]))/100)*100;
%                 max_rp = ceil(max(sites.(['rp_2018_' var_lab]))/100)*100;
                min_rp = floor(min(data_tmp)/100)*100;
                max_rp = ceil(max(data_tmp)/100)*100;

                edges = min_rp:20:max_rp;
                centers = edges(2:end) - 10;

                counts = [];
                for reg = 1:length(regions)
                    filt = strcmp(sites.region,regions{reg});
                    counts(reg,:) = histcounts(data_tmp(filt), edges);    
                end

                hold on
                for reg = 1:length(regions)
                    bar(centers, sum(counts(reg:end,:),1),'FaceColor',cmap(reg,:), 'displayName', reg_tits{reg})
                end
                title(['Site Class ' site_class{v} ', T = ' num2str(T(t1)) ' s'])
%                 title(['Risk Target = ' rt_lab{r1} '% in 50'])
                box on
                xlabel('Return Period (years)')
                ylabel('Number of Sites')
                legend('location','northeast')

%                 mkdir([save_dir filesep 'ReturnPeriod_Hist'])
%                 saveas(gcf, [save_dir filesep 'ReturnPeriod_Hist' filesep var_lab_save '.png'])
%                 mkdir(save_dir)
%                 saveas(gcf, [save_dir filesep var_lab_save '.png'])
                set(gcf, 'Color', 'None')
                set(gca, 'color', 'none');
                close
                
                % capture summary data for return periods
                if T(t1) <= 1
                    c = c + 1;
                    summary_table.risk_target(c,1) = str2double(rt_lab{r1});
                    summary_table.reliability_target(c,1) = str2double(poe_lab{r1});
                    summary_table.period(c,1) = T(t1);
                    summary_table.site_class{c,1} = site_class{v};
                    summary_table.RP_med(c,1) = nanmedian(data_tmp);
                    summary_table.RP_min(c,1) = min(data_tmp);
                    summary_table.RP_max(c,1) = max(data_tmp);
                end
                

                % Plot stacked histograms of Cr - Design
                Cr = sites.(['rtgm_2018_' var_lab]) ./ sites.(['Sa_design_' tab_lab]);
                min_rp = floor(min(Cr)*10)/10;
                max_rp = ceil(max(Cr)*10)/10;

                edges = min_rp:0.1:max_rp;
                centers = edges(2:end) - 0.05;

                counts = [];
                for reg = 1:length(regions)
                    filt = strcmp(sites.region,regions{reg});
                    counts(reg,:) = histcounts(Cr(filt), edges);    
                end

                hold on
                for reg = 1:length(regions)
                    bar(centers, sum(counts(reg:end,:),1),'FaceColor',cmap(reg,:), 'displayName', reg_tits{reg})
                end
                title(['Site Class ' site_class{v} ' - T = ' num2str(T(t1)) ' s'])
                box on
                xlabel('C_R')
                ylabel('Number of Sites')
                legend('location','northeast')

%                 mkdir([save_dir filesep 'Cr_Hist_design'])
%                 saveas(gcf, [save_dir filesep 'Cr_Hist_design' filesep var_lab_save '.png'])
                close

                % Plot stacked histograms of Cr - Prob
                Cr = sites.(['rtgm_2018_' var_lab]) ./ sites.(['Sa_rt_' tab_lab]);
                min_rp = floor(min(Cr)*10)/10;
                max_rp = ceil(max(Cr)*10)/10;

                edges = min_rp:0.1:max_rp;
                centers = edges(2:end) - 0.05;

                counts = [];
                for reg = 1:length(regions)
                    filt = strcmp(sites.region,regions{reg});
                    counts(reg,:) = histcounts(Cr(filt), edges);    
                end

                hold on
                for reg = 1:length(regions)
                    bar(centers, sum(counts(reg:end,:),1),'FaceColor',cmap(reg,:), 'displayName', reg_tits{reg})
                end
                title(['Site Class ' site_class{v} ' - T = ' num2str(T(t1)) ' s'])
                box on
                xlabel('C_R')
                ylabel('Number of Sites')
                legend('location','northeast')

%                 mkdir([save_dir filesep 'Cr_Hist_prob'])
%                 saveas(gcf, [save_dir filesep 'Cr_Hist_prob' filesep var_lab_save '.png'])
                close


                % Plot ground motion comparison - Probabilistic GMs
                figure
                title(['Site Class ' site_class{v} ' - T = ' num2str(T(t1)) ' s'])
                hold on
                for reg = 1:length(regions)
                    filt = strcmp(sites.region,regions{reg});
                    scatter(sites.(['Sa_rt_' tab_lab])(filt), sites.(['rtgm_2018_' var_lab])(filt),  'filled', 'displayName', reg_tits{reg})
                end
                lim = ceil(max([max(sites.(['Sa_rt_' tab_lab])), max(sites.(['rtgm_2018_' var_lab]))])*10)/10;
                plot([0 lim],[0 lim],'--k','handlevisibility','off')
    %             plot([0 lim],[0 lim]*0.8,':k','handlevisibility','off')
    %             plot([0 lim],[0 lim]*1.2,':k','displayName','+/- 20%')
                xlim([0 lim])
                ylim([0 lim])
                box on
                legend('location','northwest')
                xlabel(['DE_{Probabilistic} (T = ' num2str(T(t1)) ' s) (g)'])
                ylabel(['FRE_R (T = ' num2str(T(t1)) ' s) (g)'])

%                 mkdir([save_dir filesep 'qq_DErt'])
%                 saveas(gcf, [save_dir filesep 'qq_DErt' filesep var_lab_save '.png'])
                close

                % Plot return period comparison - Probabilistic GMs
                figure
                title(['Site Class ' site_class{v} ' - T = ' num2str(T(t1)) ' s'])
                hold on
                for reg = 1:length(regions)
                    filt = strcmp(sites.region,regions{reg});
                    scatter(sites.(['RP_rt_' tab_lab])(filt), sites.(['rp_2018_' var_lab])(filt),  'filled', 'displayName', reg_tits{reg})
                end
                lim = ceil(max([max(sites.(['RP_rt_' tab_lab])), max(sites.(['rp_2018_' var_lab]))])*10)/10;
                plot([0 lim],[0 lim],'--k','handlevisibility','off')
    %             plot([0 lim],[0 lim]*0.8,':k','handlevisibility','off')
    %             plot([0 lim],[0 lim]*1.2,':k','displayName','+/- 20%')
                xlim([0 lim])
                ylim([0 lim])
                box on
                legend('location','northwest')
                xlabel('Return Period DE_{probabilistic} (Years)')
                ylabel('Return Period FRE_R (Years)')

%                 mkdir([save_dir filesep 'qq_RPrt'])
%                 saveas(gcf, [save_dir filesep 'qq_RPrt' filesep var_lab_save '.png'])
                close

                
                % Cap values at DE
                data_tmp = min(sites.(['rtgm_2018_' var_lab]),sites.(['Sa_design_' tab_lab]),'includenan');
                
                
                % Plot ground motion comparison - DE
                figure
                title(['Site Class ' site_class{v} ' - T = ' num2str(T(t1)) ' s'])
                hold on
                for reg = 1:length(regions)
                    filt = strcmp(sites.region,regions{reg});
                    scatter(sites.(['Sa_design_' tab_lab])(filt), data_tmp(filt),  'filled', 'displayName', reg_tits{reg})
                end
                lim = ceil(max([max(sites.(['Sa_design_' tab_lab])), max(sites.(['rtgm_2018_' var_lab]))])*10)/10;
                plot([0 lim],[0 lim],'--k','handlevisibility','off')
    %             plot([0 lim],[0 lim]*0.8,':k','handlevisibility','off')
    %             plot([0 lim],[0 lim]*1.2,':k','displayName','+/- 20%')
                xlim([0 lim])
                ylim([0 lim])
                box on
                legend('location','northwest')
                xlabel(['DE (T = ' num2str(T(t1)) ' s) (g)'])
                ylabel(['FRE_R (T = ' num2str(T(t1)) ' s) (g)'])

%                 mkdir([save_dir filesep 'qq_DE'])
%                 saveas(gcf, [save_dir filesep 'qq_DE' filesep var_lab_save '.png'])
                set(gcf, 'Color', 'None')
                set(gca, 'color', 'none');
                close

                % Plot return period comparison - DE
                figure
                title(['Site Class ' site_class{v} ' - T = ' num2str(T(t1)) ' s'])
                hold on
                for reg = 1:length(regions)
                    filt = strcmp(sites.region,regions{reg});
                    scatter(sites.(['RP_design_' tab_lab])(filt), sites.(['rp_2018_' var_lab])(filt),  'filled', 'displayName', reg_tits{reg})
                end
                lim = ceil(max([max(sites.(['RP_design_' tab_lab])), max(sites.(['rp_2018_' var_lab]))])*10)/10;
                plot([0 lim],[0 lim],'--k','handlevisibility','off')
    %             plot([0 lim],[0 lim]*0.8,':k','handlevisibility','off')
    %             plot([0 lim],[0 lim]*1.2,':k','displayName','+/- 20%')
                xlim([0 lim])
                ylim([0 lim])
                box on
                legend('location','northwest')
                xlabel('Return Period DE (Years)')
                ylabel('Return Period FRE_R (Years)')

%                 mkdir([save_dir filesep 'qq_RP'])
%                 saveas(gcf, [save_dir filesep 'qq_RP' filesep var_lab_save '.png'])
                close



            end
        end

        % Save tabulated values
%         writetable(fre_table,[save_dir filesep 'fre_table_' tab_lab '.csv'])
%         writetable(rp_table,[save_dir filesep 'rp_table_' tab_lab '.csv'])
%         writetable(cr_de_table,[save_dir filesep 'cr_de_table_' tab_lab '.csv'])
%         writetable(cr_rt_table,[save_dir filesep 'cr_rt_table_' tab_lab '.csv'])
    end
end
% 
% for v = 1:length(vs30)
%     cr_de_table_v2 = sites(:,1:4);
%     
%     for t1 = 1:length(T)
%     
%         tab_lab = (['vs30_' num2str(vs30(v)) '_T_' T_lab{t1}]);
%         var_lab = (['beta_0p7' '_vs30_' num2str(vs30(v)) '_T_' T_lab{t1} '_rt_10_PoE_25']);
%         
%         cr_de_table_v2.(['T_' T_lab{t1}]) = sites.(['rtgm_2018_' var_lab]) ./ sites.(['Sa_design_' tab_lab]) ;
%     end
%     
%     % Save tabulated values
%     writetable(cr_de_table_v2,[plot_dir filesep 'cr_table_160sites' '_vs30_' num2str(vs30(v)) '.csv'])
% end



