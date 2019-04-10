%% Script to Pull Design Hazard info from USGS web services
clear
close all
clc
rehash
rng('shuffle')

%% Load inputs data
site = readtable([pwd filesep 'inputs.csv'],'ReadVariableNames',true);

%% Pull Down info for each site in table
tic
for i = 1:height( site )
    [ design_values ] = fn_call_USGS_design_value_API( i, site.reference_doc{i}, site.lat(i), site.lng(i), site.risk_cat{i}, site.site_class{i} );
    design_table(i,:) = struct2table(design_values);
end
toc

%% Merge Site and Design Value Table
site = [site, design_table];

%% Save data back to table
writetable(site,[pwd filesep 'outputs.csv']);

