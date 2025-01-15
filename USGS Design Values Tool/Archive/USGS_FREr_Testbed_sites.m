clear all
close all
clc

%% Inputs
% read_hazard_data = true;

% maxdir_factor = [1.200	1.200	1.200	1.200	1.200	1.200	1.200	1.200	1.200	1.203	1.206	1.213	1.219	1.234	1.250	1.253	1.256	1.261	1.267	1.272	1.286	1.300]; % based on USGS conversion from Nico

% vs30 = [530 760 265];
% site_class = {'C' 'BC', 'D'};

% vs30 = [530];
% site_class = {'C'};


% Imports
import usgs.*

% Load model table
sites = readtable('testbed_sites.csv');

% Set web read options
options = weboptions('Timeout', 30);

% Set plot dir
% plot_dir = 'C:\Users\dtc2\OneDrive - NIST\PUC\TS5\FRE Explore\Spectra_Plots';
% mkdir(plot_dir)

% Set matlab colors
cmap = [0 0.4470 0.7410;
        0.8500 0.3250 0.0980;
        0.9290 0.6940 0.1250;
        0.4940 0.1840 0.5560;
        0.4660 0.6740 0.1880;
        0.3010 0.7450 0.9330;
        0.6350 0.0780 0.1840];

%% Load hazard data for RTGM
beta = 0.6;
PoE = 0.25;
rt = 0.1;
t = 50;
nshm = 2018;
for s = 1:height(sites)
    s

    % Pull design values for this site
    [design_values, MPS] = fn_call_USGS_design_value_API_v2(1, 'asce7-22', sites.lat(s), sites.lon(s), 'II', sites.site_class{s});

    % Read hazard data for this site
    DATA = webread( sprintf('https://staging-earthquake.usgs.gov/ws/nshmp/conus-%i/dynamic/hazard/%0.2f/%0.2f/%i',...
                    nshm, ...
                    sites.lon(s), ...
                    sites.lat(s), ...
                    sites.vs30(s)), ...
                    options);

    T = MPS.designSpectrum.periods;

    % Pull spectra data
    for t1 = 1:length(T)
        found_vals = 0;
        gms = [];
        afes = [];
        rps = []; 

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

        % Calculate and save risk targeted ground motions
        try
            [rtgm] = fn_calc_rtgm(gms, afes, beta, PoE, rt, t);
            MPS.designSpectrum.fre_r(t1,1) = rtgm;
        catch
            error('TEST ERROR');
        end
    end
       
    % Develop the DE 2-point spectrum
    sds = design_values.sds;
    sd1 = design_values.sd1;
    Ts = sd1/sds;
    T0 = 0.2*(sd1/sds);
    t_range = (Ts+0.1):0.1:5;
    code_spectrum_T = [0 T0 Ts t_range];
    code_spectrum = [0.4*sds sds sds sd1./t_range];
    
    % Develop The FREr 2-point spectrum
    frer_spectra = MPS.designSpectrum.fre_r;
    T_range_s = MPS.designSpectrum.periods >= 0.2 & MPS.designSpectrum.periods <= 5.0;
    sdsfr = 0.9*max(frer_spectra(T_range_s));
    T_range_1 = MPS.designSpectrum.periods >= 1.0 & MPS.designSpectrum.periods <= 2.0;
    sd1fr = max(frer_spectra(MPS.designSpectrum.periods == 1.0),...
                0.9*max(MPS.designSpectrum.periods(T_range_1).*frer_spectra(T_range_1)));
    Ts = sd1fr/sdsfr;
    T0 = 0.2*(sd1fr/sdsfr);
    t_range = (Ts+0.1):0.1:5;
    code_spectrum_T_frer = [0 T0 Ts t_range];
    code_spectrum_frer = [0.4*sdsfr sdsfr sdsfr sd1fr./t_range];

    % Save or plot spectra
    figure
    hold on
    plot(MPS.designSpectrum.periods,MPS.designSpectrum.fre_r,'-k','DisplayName','FRE_R - MPS')
    plot(code_spectrum_T_frer,code_spectrum_frer, '--k', 'DisplayName','FRE_R - 2-point Spectrum')
    plot(MPS.designSpectrum.periods,MPS.designSpectrum.ordinates,'color',cmap(1,:),'DisplayName','DE - MPS')
    plot(code_spectrum_T,code_spectrum, '-','color',cmap(1,:),'DisplayName','DE - 2-point Spectrum')
%             plot(MPS{2}.designSpectrum.periods,MPS{2}.designSpectrum.ordinates,'color',cmap(2,:),'DisplayName','DE - Site Class D')
%             plot(MPS{2}.designSpectrum.periods,MPS{2}.designSpectrum.fre_r,'--','color',cmap(2,:),'DisplayName','FRE_R - Site Class D')
    xlim([0 5])
    xlabel('Period (s)')
    ylabel('Sa (g)')
    legend('location','northeast')
    box on
    grid on
%     title('FRE_R Geomean')
%     saveas(gcf,[plot_dir filesep sites.name{s} '-geomean.png'])
    close

    % Save rtgm to table
    sites.sds(s) = sds;
    sites.sd1(s) = sd1;
    sites.sdsfr(s) = sdsfr;
    sites.sd1fr(s) = sd1fr;
    for T1 = 1:length(MPS.designSpectrum.periods)
        t_num = strrep(sprintf('%0.2f',MPS.designSpectrum.periods(T1)),'.','p');
        sites.(['frer_mps_' t_num])(s) = MPS.designSpectrum.fre_r(T1);
    end

    % Write table
    writetable(sites,'testbed_sites.csv')
end








