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

%%%%%% NEED TO RUN THIS FROM AN ID SCRIPT INSTEAD OF OVERWRITTING EACH TIME %%%%%%%%%

%% Go through each site and call USGS API
id = 0;
for s = 1:height(sites)

    % Set default site classes to pull
    if ~isempty(sites.site_class{s})
        site_class = sites.site_class(s);
        vs30 = sites.vs30(s); %must be the same length as site clas
    else
        site_class = {'BC', 'C', 'D', 'E'};
        vs30 = [760 530 265 150];

%             site_class = {'A'  'B'  'BC' 'C' 'CD' 'D' 'DE' 'E'};
%             vs30 =       [1500 1080 760  530 365  265 180  150];
    end

    % Call USGS for each applicable site class and save
    for v = 1:length(vs30)
        id = id + 1;

        if any(ismember(sites.Properties.VariableNames,'log')) && sites.log(s) == 0
        call_success = 1;

        % build output data structure
        out(id).name = sites.name{s};
        out(id).lat = sites.lat(s);
        out(id).lon = sites.lon(s);
        out(id).region = sites.region{s};
        out(id).sdc = sites.sdc{s};
        out(id).usgs_160_id = sites.usgs_160_id(s);
        out(id).nehrp_34_id = sites.nehrp_34_id(s);
        out(id).frtc_testbed_id = sites.frtc_testbed_id(s);
        out(id).site_class = site_class{v};
        out(id).vs30 = vs30(v);

        % Pull design values for this site
        if ~isfield(out,'design_values') || isempty(out(id).design_values)
            [design_values, MPS, status_design] = fn_call_USGS_design_value_API(1, 'asce7-22', sites.lat(s), sites.lon(s), 'II', site_class{v});
            fprintf('ASCe 7-22 Desgin Call: %s\n', status_design)
            if strcmp(status_design,'success')
                out(id).design_values = design_values;
                out(id).design_values.MPS = MPS;
                pause(pause_time)
            else
                pause(4*pause_time) % pause extra long when erros occur (come back and finish later)
            end
        end

        % Read 2018 hazard data for this site
        if ~isfield(out,'hazard_curve') || ~isfield(out(id).hazard_curve,'nshm_2018') || isempty(out(id).hazard_curve.nshm_2018)
            [DATA_2018, status_2018] = fn_call_USGS_hazard_API(2018, sites.lat(s), sites.lon(s), vs30(v));
            fprintf('2018 Hazard Call: %s\n', status_2018)
            if strcmp(status_2018,'success')
                out(id).hazard_curve.nshm_2018 = DATA_2018.response.hazardCurves;
                pause(pause_time)
            elseif contains(status_2018,'Too Many Requests')
                pause(4*pause_time) % pause extra long when this erros occur
            else
                pause(pause_time)
            end
        end

        % Read 2023 hazard data for this site
        if ~isfield(out,'hazard_curve') || ~isfield(out(id).hazard_curve,'nshm_2023') || isempty(out(id).hazard_curve.nshm_2023)
            [DATA_2023, status_2023] = fn_call_USGS_hazard_API(2023, sites.lat(s), sites.lon(s), vs30(v));
            fprintf('2023 Hazard Call: %s\n', status_2023)
            if strcmp(status_2023,'success')
                out(id).hazard_curve.nshm_2023 = DATA_2023.response.hazardCurves;
                pause(pause_time)
            elseif contains(status_2023,'Too Many Requests')
                pause(4*pause_time) % pause extra long when this erros occur
            else
                pause(pause_time)
            end
        end

        end
    end
%     % Save data for this site 
%     tic
%     save('DATA_usgs_api.mat','out')
%     toc
end

% Save data for this site 
% tic
% save('DATA_usgs_api.mat','out')
% toc

