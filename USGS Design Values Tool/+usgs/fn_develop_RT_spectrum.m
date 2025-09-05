function [rt_spectrum] = fn_develop_RT_spectrum(hazard_curve, rt_params)
% Description: Calculate a risk targeted spectrum from USGS hazard curve.

% Created by: Dustin Cook
% Date Created: 1/15/2025

% Inputs:
%   hazard_curve - hazard curves from USGS api call
%   rt_params - data structure of risk targetting inputs / control points.

% Ouputs:
%   rt_spectrum = Risk targeted spectrum, defined by Spectral accleration
%   values at a given set of period points.

% functions
import usgs.fn_calc_rtgm

% Set fixed Period vector
T_vec =  [0,   0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25,  0.3,   0.4,   0.5,   0.75,  1,     1.5,   2,     3,     4,     5,     7.5,   10];
maxdir = [1.2, 1.2, 1.2,   1.2,  1.2,  1.2,   1.2, 1.2,  1.2, 1.203, 1.206, 1.213, 1.219, 1.234, 1.250, 1.253, 1.256, 1.261, 1.267, 1.272, 1.286, 1.3]; % based on USGS conversion from Nico

% Pull spectrum from hazard curve
for t1 = 1:length(T_vec)
    found_vals = 0;
    gms = [];
    afes = [];
    rtgm = [];

    % Pull hazard curve for this period
    for i = 1:length(hazard_curve)
        spectral_label = hazard_curve(i).imt.value;   
        spectral_value = str2double(hazard_curve(i).imt.display(1:end-2));

        if (strcmp(spectral_label,'PGA') && T_vec(t1) == 0) || (contains(spectral_label,'SA') && spectral_value == T_vec(t1))
            found_vals = 1;
            gms = hazard_curve(i).data(1).values.xs';
            afes = hazard_curve(i).data(1).values.ys';
            filt = afes > 0;
            gms = gms(filt);
            afes = afes(filt);
        end
    end

    if found_vals == 0
        error('couldnt find hazard values for given period')
    end

    % Calculate and save risk targeted ground motions
    try
        [rtgm] = fn_calc_rtgm(gms, afes, rt_params.beta, rt_params.PoE, rt_params.rt, rt_params.t);
    catch
        error('failed to converge on risk target')
    end

    if rt_params.maxdir
        rt_spectrum.sa(t1) = maxdir(t1)*rtgm; % Convert to max direction
    else
        rt_spectrum.sa(t1) = rtgm;
    end
    
    rt_spectrum.return_period(t1) = 1/exp(interp1(gms,log(afes),rtgm));
    rt_spectrum.return_period_67(t1) = 1/exp(interp1(gms,log(afes),rtgm*2/3));
end

rt_spectrum.period = T_vec;