%% Calls USGS API ro pulls hazard curves and design values for specified sites
% Saves as matlab data (probably should save as JSON). Function will bypass
% site that fail the API call and story them as FALSE in the site.log
% field. Run the script again to retry all failures.
clear all
close all
clc

%% Inputs
% Load model table
sites = readtable('sites.csv'); 
if exist('DATA_usgs_api.mat','file')
    load('DATA_usgs_api.mat')
end

% Set pause time after each call (rate limit calls b/c usgs will start
% throwing "too many requests" errors.
pause_time = 30; % seconds

%% Imports
import usgs.*

%% Go through each site and call USGS API
for s = 1:height(sites)
    if any(ismember(sites.Properties.VariableNames,'log')) && sites.log(s) == 0
        s

        call_success = 1;
        
        % build output data structure
        out(s).name = sites.name{s};
        out(s).lat = sites.lat(s);
        out(s).lon = sites.lon(s);
        out(s).region = sites.region{s};
        out(s).sdc = sites.sdc{s};
        out(s).usgs_160_id = sites.usgs_160_id(s);
        out(s).nehrp_34_id = sites.nehrp_34_id(s);
        out(s).frtc_testbed_id = sites.frtc_testbed_id(s);

        % Set default site classes to pull
        if ~isempty(sites.site_class{s})
            site_class = sites.site_class(s);
        else
            site_class = {'BC', 'C', 'D'};
        end
        
        % Set default vs30s to pull
        if ~isnan(sites.vs30(s))
            vs30 = sites.vs30(s); %must be the same length as site clas
        else
            vs30 = [760 530 265];
        end
        
        % Save site class and vs30 to out data
        out(s).site_class = site_class;
        out(s).vs30 = vs30;
        
        for v = 1:length(vs30)
            if call_success
                % Pull design values for this site
                [design_values, MPS, status_design] = fn_call_USGS_design_value_API(1, 'asce7-22', sites.lat(s), sites.lon(s), 'II', site_class{v});
                pause(pause_time)
                
                % Read 2018 hazard data for this site
                [DATA_2018, status_2018] = fn_call_USGS_hazard_API(2018, sites.lat(s), sites.lon(s), vs30(v));
                pause(pause_time)
                
                % Read 2018 hazard data for this site
                [DATA_2023, status_2023] = fn_call_USGS_hazard_API(2023, sites.lat(s), sites.lon(s), vs30(v));
                pause(pause_time)

                if strcmp(status_design,'success') && strcmp(status_2018,'success') && strcmp(status_2023,'success') % If all succeed
                    % Save data into matlab struture
                    out(s).design_values.(['site_class_' site_class{v}]) = design_values;
                    out(s).design_values.(['site_class_' site_class{v}]).MPS = MPS;
                    out(s).hazard_curve.nshm_2018.(['vs30_' num2str(vs30(v))]) = DATA_2018.response.hazardCurves;
                    out(s).hazard_curve.nshm_2023.(['vs30_' num2str(vs30(v))]) = DATA_2023.response.hazardCurves;

                    save('DATA_usgs_api.mat','out')
                else
                    % Calls failed, come back later to finish this off
                    call_success = 0;
                    
                    fprintf('ASCe 7-22 Desgin Call: %s\n', status_design)
                    fprintf('2018 Hazard Call: %s\n', status_2018)
                    fprintf('2023 Hazard Call: %s\n', status_2023)
                    
                    if ~strcmp(status_design,'success')
                        pause(4*pause_time) % pause extra long when erros occur
                    end
                    
                    if contains([status_hazard2018,status_hazard2023],'Too Many Requests')
                        pause(4*pause_time) % pause extra long when erros occur
                    end
                end
            end
        end
        
        % log success rate
        sites.log(s) = call_success;
        writetable(sites,'sites.csv')
        save('DATA_usgs_api.mat','out')
    end
end


