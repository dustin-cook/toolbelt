function [sds, sd1] = fn_develop_2point_design_values(sa, period)
% Description: Calculate two-point spectral design values according to ASCE
% 7-22 21.4

% Created by: Dustin Cook
% Date Created: 1/15/2025

% Inputs:
%   sa - vector of sa values of design spectrum
%   period - vector of period values of design spectrum

% Ouputs:
%   sds = ASCE 7 short period spectral design value
%   sd1 = ASCE 7 1-second period spectral design value

% Sds
T_range_s = period >= 0.2 & period <= 5.0;
sds = 0.9*max(sa(T_range_s));

% Sd1
T_range_1 = period >= 1.0 & period <= 2.0;
sd1 = max(sa(period == 1.0),...
            0.9*max(period(T_range_1).*sa(T_range_1)));
                
end