function [rtgm] = fn_calc_rtgm(sa, afe, beta, PoE, risk_target, t)
% Description: Integrate a fragility shape with the hazard curve to develop
% risk targeted ground motions defined as a specific anchor point

% Created by: Dustin Cook
% Date Created: 1/15/2025

% Inputs:
%   sa - sa values of hazard curve (vector; period dependant)
%   afe - annaul frequence of exceedance corresponding to each Sa value
%   beta - shape of the fragility curve, as a lognormal standard deviation
%   PoE - probability of exceedance given targe shaking
%   (e.g., 0.10 for 10% chance given rtgm)
%   risk_target - target risk of exceedance across specified timespan
%   (e.g., 0.01 for 1% risk)
%   t - timespan of risk in years (e.g. design life of building)

% Ouputs:
%   rtgm = Sa value of the risk targeted ground motion (period dependant)

% Notes:
%   I dont think this function is currenlty being used...

%% Intial Parameters
% Define Search Grid
x = sa(1):0.01:sa(end);

%% Define fragility
z = norminv(PoE,0,1);
theta = exp(log(x) - z*beta);
frag_pdf = lognpdf(sa,log(theta'),beta);

%% Integrate hazard curve and fragility to calculate RTGM
haz_x_frag_pdf = frag_pdf .* afe;

% end_point_contribution = phaz_x_frag_pdf(:,end) ./ sum(haz_x_frag_pdf(:,1:(end-1)),2);
% check = end_point_contribution > 0.01; % make sure it doesnt represent more than 1% of the relative risk
% if any(check)
%     error('NOOOO LAMA NOOOOOOOOOO!')
% end

cum_int = cumtrapz(sa,haz_x_frag_pdf');
cum_afe = cum_int(end,:);

% Calc the lifetime rist
risk = 1 - (1 - cum_afe).^t;

% Find RTGM
[~,filt,~] = unique(risk,'stable'); % Find unique risk points for interpolation
rtgm = exp(interp1(risk(filt),log(x(filt)),risk_target,'linear','extrap'));

end