function [ coversion_factor ] = fn_max_to_geomean( T )
% Description: find the conversion factor between max direction Sa and
% geomean Sa values (period dependant)

% Created by: Dustin Cook
% Date Created: 1/15/2025

% Inputs:
%   T - period of interest in seconds

% Ouputs:
%   coversion_factor = ratio of max direction Sa to geomean Sa 

% Notes:
%   I dont think this function is currenlty being used...

%% Begin Function
% Calc the coversion factor currently used in SP3 (developed by jack, curt
% or katie)
if T<= 0.3
    coversion_factor = 1.1;
elseif T<= 0.5
    coversion_factor = 1.1 + 0.5*(T-0.3);
elseif T <= 1.0
    coversion_factor = 1.2 + 0.2*(T-0.5);
elseif T <= 2.0
    coversion_factor = 1.3;
elseif T <= 4.0
    coversion_factor = 1.3 + 0.05*(T-2);
else
    coversion_factor = 1.4;
end


end

