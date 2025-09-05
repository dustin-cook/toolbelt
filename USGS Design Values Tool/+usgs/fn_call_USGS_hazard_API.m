function [DATA, status] = fn_call_USGS_hazard_API(edition, lat, lon, vs30)
% Description: Function to Call the USGG Hazard Curve API.

% Created by: Dustin Cook
% Date Created: 1/15/2025

% Inputs:
%   edition - NSHM release year
%   lat - numeric lattitide of the site
%   lng - numeric longitude of the site
%   vs30 - average shear wave velocity at 30m for site (in m/s)

% Ouputs:
%   DATA = returns hazard curve data. Empty if call fails.
%   status = [sting] status of webcall

% Set parameters
options = weboptions('Timeout', 30);         

try
    %% Call USGS Hazards API
    DATA = webread( sprintf('https://earthquake.usgs.gov/ws/nshmp/conus-%i/dynamic/hazard/%0.2f/%0.2f/%i',...
                                    edition, ...
                                    lon, ...
                                    lat, ...
                                    vs30), ...
                                    options);
    status = DATA.status;
catch ME
    DATA = [];
    status = ME.message;
end

end

