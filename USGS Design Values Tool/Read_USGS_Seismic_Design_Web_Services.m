%% Script to Pull Design Hazard info from USGS web services
clear
close all
clc
rehash
rng('shuffle')

%% Load inputs data
site = readtable([pwd filesep 'Design Values.csv'],'ReadVariableNames',true);

%% Pull Down info for each site in table
for i = 1:length( site.id )
    disp( [ int2str(i) ') ' site.city{i} ] )
    
    DATA = webread( [ 'https://earthquake.usgs.gov/ws/designmaps/' ...
                      site.reference_doc{i} '.json?' ], ...
                      'latitude', site.lat(i), ...
                      'longitude', site.lng(i), ...
                      'riskCategory', site.risk_cat{i}, ...
                      'siteClass', site.site_class{i}, ...
                      'title', site.city{i} );
            
    site.ss(i) = DATA.response.data.ss;
    site.s1(i) = DATA.response.data.s1;
    site.fa(i) = DATA.response.data.fa;
    site.sms(i) = DATA.response.data.sms;
    site.sds(i) = DATA.response.data.sds;
    
    % Check if Fv is empty and why
    if isempty(DATA.response.data.fv)
        site.fv{i} = DATA.response.data.fv_note;
        site.sm1{i} = 'NA';
        site.sd1{i} = 'NA';
        site.sdc{i} = 'NA';
    else
        site.fv{i} = DATA.response.data.fv;
        site.sm1{i} = DATA.response.data.sm1;
        site.sd1{i} = DATA.response.data.sd1;
        site.sdc{i} = DATA.response.data.sdc;
    end
    
    % Check to see if risk targeted values exist
    if isfield(DATA.response.data,'ssrt') && isfield(DATA.response.data,'s1rt')
        site.ssrt{i} = DATA.response.data.ssrt;
        site.s1rt{i} = DATA.response.data.s1rt;
    else
        site.ssrt{i} = 'NA';
        site.s1rt{i} = 'NA';
    end
    
    % Check to see if Uniform Hazard values exist
    if isfield(DATA.response.data,'ssuh') && isfield(DATA.response.data,'s1uh')
        site.ssuh{i} = DATA.response.data.ssuh;
        site.s1uh{i} = DATA.response.data.s1uh;
    else
        site.ssuh{i} = 'NA';
        site.s1uh{i} = 'NA';
    end

end

%% Save data back to table
writetable(site,[pwd filesep 'Design Values.csv']);

