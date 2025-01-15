clear all
close all
clc

% Inputs

vs30 = 760;
site_class = 'BC';

return_periods = [14, 43, 72, 108, 144, 224, 475, 975, 2475, 5000];
afe = 1 ./ return_periods;

% Imports
import usgs.*

% Load model table
models = readtable('nehrp34_RCSW_models.csv');

% for m = 1:height(models)
%     options = weboptions('Timeout', 30);
%     
%     DATA = webread( sprintf('https://earthquake.usgs.gov/nshmp-haz-ws/deagg/E2014/COUS/%0.2f/%0.2f/SA1P0/760/%0.0f',...
%                     models.lng(m),...
%                     models.lat(m), ...
%                     models.rp_10_2(m)), ...
%                   options);
% 
%     models.mean_Mw(m) = DATA.response.data(1).summary(4).data(1).value;
%     models.mode_Mw(m) = DATA.response.data(1).summary(5).data(1).value;
% end

% % Go through each model and pull design values 
% for m = 1:height(models)
%     m
%     % Collect Design Info from USGS
%     [design_values] = fn_call_USGS_design_value_API(1, 'asce7-22', models.lat(m), models.lng(m), 'II', site_class);
%     
%     % Save design values to model table
%     models.sd1(m) = design_values.sd1;
%     models.sd184(m) = design_values.sd184;
%     models.sd1rt(m) = design_values.sd1rt;
%     models.sds(m) = design_values.sds;
%     models.sds84(m) = design_values.sds84;
%     models.sdsrt(m) = design_values.sdsrt;
% end
% 
% % Save model table with design values
% writetable(models,'nehrp34_RCSW_models.csv');
    
% Load hazard table
% sa_vals = readtable('nehrp34_hazard_1sec.csv');

sa_vals = readtable('nehrp34_hazard_0p2sec.csv');

% % Go through each model and pull hazard values 
% for m = 1:height(models)
%     m
%     % Collect Hazard and Design Info
%     [sa_spectra, ~] = fn_call_USGS_hazard_API('E2014', models.lat(m), models.lng(m), vs30, afe);
% 
%     % Get 1 second hazard curve
% %     idx = sa_periods == 1;
%     sa_vals{m,2:end} = sa_spectra;
% end
% 
% % Save model table with design values
% writetable(sa_vals,'nehrp34_hazard_0p2sec.csv');

% Convert design values to RotD50
% coversion_factor = 1/1.3; % for T = 1 based on JWB code 
coversion_factor = 1/1.2; % for T = 1 based on USGS conversion from Nico

% Calculate risk targeted 
% betas = [0.5 0.55 0.6 0.7 0.9];
betas = 0.7;
PoE = [0.10 0.10 0.10 0.25 0.25];
rt = [0.01 0.02 0.05 0.05 0.1];
for m = 1:height(models)
    sa = sa_vals{m,2:end};
    filt = ~isnan(sa);
    models.rp_DE(m,1) = 1/interp1(sa(filt),afe(filt),models.sdsrt(m));
    for a = 1:length(PoE)
        [rtgm] = fn_calc_rtgm(sa(filt), afe(filt), betas, PoE(a), rt(a), 50);
        rtgm = 2/3*rtgm;
        models.(['rtgm_' num2str(100*PoE(a)) '_' num2str(100*rt(a))])(m,1) = rtgm;
        models.(['cr_' num2str(100*PoE(a)) '_' num2str(100*rt(a))])(m,1) = rtgm / (models.sdsrt(m)*coversion_factor);
        models.(['rp_' num2str(100*PoE(a)) '_' num2str(100*rt(a))])(m,1) = 1/interp1(sa(filt),afe(filt),rtgm);
    end
end

% Save model table with design values
% writetable(models,'nehrp34_RCSW_models_0p2.csv');

% Load model table
% models = readtable('nehrp34_RCSW_models_0p2.csv');
% models = readtable('nehrp34_RCSW_models.csv');
   
% 10 and 2
figure
title('10% Reliability - 1% Risk')
% regions = unique(models.region);
regions = {'so cal' 'nor cal' 'WUS' 'PNW' 'CEUS'};
reg_tits = {'SoCal' 'NorCal' 'Other WUS' 'Pacific NW' 'CEUS'};
hold on
for r = 1:length(regions)
    filt = strcmp(models.region,regions{r});
    scatter(models.rp_DE(filt), models.rp_10_1(filt), 'filled', 'displayName',reg_tits{r})
%     scatter(models.rp_DE(filt) .* coversion_factor, models.rp_25_5(filt), 'filled', 'displayName',reg_tits{r})
end
lim = 1200;
plot([0 lim],[0 lim],'--k','handlevisibility','off')
plot([0 lim],[0 lim]*0.85,':k','handlevisibility','off')
plot([0 lim],[0 lim]*1.15,':k','displayName','+/- 15%')
xlim([0 lim])
ylim([0 lim])
box on
legend('location','northwest')
% xlabel('S_{DS} (g)')
% ylabel('FRE_R (g)')
% legend('location','southeast')
xlabel('RP_{2/3 MCE_R} (years)')
ylabel('RP_{Functional Recovery} (years)')

% % 25 and 5
% figure
% title('25% Reliability - 5% Risk')
% % regions = unique(models.region);
% hold on
% for r = 1:length(regions)
%     filt = strcmp(models.region,regions{r});
%     scatter(models.sdsrt(filt) .* coversion_factor, models.rtgm_25_5(filt), 'filled', 'displayName',regions{r})
% end
% lim = 2;
% plot([0 lim],[0 lim],'--k','handlevisibility','off')
% plot([0 lim],[0 lim]*0.85,':k','handlevisibility','off')
% plot([0 lim],[0 lim]*1.15,':k','displayName','+/- 15%')
% % xlim([0 0.7])
% % ylim([0 0.7])
% box on
% legend('location','northwest')
% xlabel('S_{DS} (g)')
% ylabel('FRE_R (g)')

mean(models.cr_10_2)
mean(models.cr_10_5)
mean(models.cr_25_5)
mean(models.cr_25_10)

mean(models.rp_10_2)
mean(models.rp_10_5)
mean(models.rp_25_5)
mean(models.rp_25_10)

tmp1 = mean(cr,2)
tmp2 = mean(rp,2)

tmp = cr'

figure
h = histogram(models.cr_10_2);
h.EdgeColor = 'none';
ax = gca;
ax.YAxis.Visible = 'off';
box off
xlabel('C_r')
xlim([0 1.5])
set(gcf,'position',[0,0,300,250])

figure
h = histogram(models.cr_25_5);
h.EdgeColor = 'none';
ax = gca;
ax.YAxis.Visible = 'off';
box off
xlabel('C_r')
xlim([0 1.5])
set(gcf,'position',[0,0,300,250])



Sds = design_values.sds;
Sd1 = design_values.sd1;

% Sds = 1;
% Sd1 = 0.6;
Tds = Sd1 / Sds;
design_spectra = Sd1 ./ sa_periods;
design_spectra(sa_periods <= Tds) = Sds;


% Plot Values
plot(sa_periods,sa_spectra)
hold on
plot(sa_periods,design_spectra)

ratio = sa_spectra./design_spectra

