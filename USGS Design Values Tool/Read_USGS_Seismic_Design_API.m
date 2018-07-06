%% Script to Pull Design Hazard info from USGS web services
clear
close all
clc
rehash
rng('shuffle')

%% Load inputs data
site_inputs = readtable([pwd filesep 'inputs.csv'],'ReadVariableNames',true);
site_coefs = readtable([pwd filesep 'asce_site_coefficients.csv'],'ReadVariableNames',true);

% Set up outputs table
site = table(site_inputs.id,site_inputs.lat,site_inputs.lng,'VariableNames',{'id','lat','lng'});

%% Pull Down info for each site in table
options = weboptions;
options.Timeout = 30;
for i = 1:length( site_inputs.id )
    DATA = webread( [ 'https://earthquake.usgs.gov/ws/designmaps/' ...
                      site_inputs.reference_doc{i} '.json?' ], ...
                      'latitude', site_inputs.lat(i), ...
                      'longitude', site_inputs.lng(i), ...
                      'riskCategory', site_inputs.risk_cat{i}, ...
                      'siteClass', site_inputs.site_class{i}, ...
                      'title', site_inputs.id(1), options );
            
    
    site.ss(i) = DATA.response.data.ss;
    site.s1(i) = DATA.response.data.s1;

    % Check if Fa is empty and why
    if isempty(DATA.response.data.fa)
        % If USGS API fails to calculate site coeffs then calc manually for ASCE references
        if contains(site_inputs.reference_doc{i},'asce') && any(strcmp(site_coefs.version,site_inputs.reference_doc{i}))            
            filter = strcmp(site_coefs.site_coef,'f_a') & strcmp(site_coefs.version,site_inputs.reference_doc{i});
            fa_array = site_coefs.(site_inputs.site_class{i})(filter);
            ss_array = site_coefs.(site_inputs.site_class{i})(filter);
            site.fa(i) = interp1(ss_array,fa_array,site.ss(i));
            site.sms(i) = site.fa(i)*site.ss(i);
            site.sds(i) = (2/3)*site.sms(i);
        else
            site.fa{i} = DATA.response.data.fa_note;
            site.sms{i} = 'NA';
            site.sds{i} = 'NA';
        end
    else
        site.fa{i} = DATA.response.data.fa;
        site.sms{i} = DATA.response.data.sms;
        site.sds{i} = DATA.response.data.sds;
    end
    
    % Check if Fv is empty and why
    if isempty(DATA.response.data.fv)
        % If USGS API fails to calculate site coeffs then calc manually for ASCE reference
        if contains(site_inputs.reference_doc{i},'asce') && any(strcmp(site_coefs.version,site_inputs.reference_doc{i}))
            filter = strcmp(site_coefs.site_coef,'f_v') & strcmp(site_coefs.version,site_inputs.reference_doc{i});
            fv_array = site_coefs.(site_inputs.site_class{i})(filter);
            s1_array = site_coefs.(site_inputs.site_class{i})(filter);
            site.fv(i) = interp1(s1_array,fv_array,site.s1(i));
            site.sm1(i) = site.fv(i)*site.s1(i);
            site.sd1(i) = (2/3)*site.sms(i);
        else
            site.fv{i} = DATA.response.data.fv_note;
            site.sm1{i} = 'NA';
            site.sd1{i} = 'NA';
        end
    else
        site.fv{i} = DATA.response.data.fv;
        site.sm1{i} = DATA.response.data.sm1;
        site.sd1{i} = DATA.response.data.sd1;
      
    end
    
    % Check if SDC is empty
    if isempty(DATA.response.data.sdc)
        site.sdc{i} = 'NA';
    else
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
writetable(site,[pwd filesep 'outputs.csv']);

