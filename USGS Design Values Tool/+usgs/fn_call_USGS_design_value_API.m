function [design_values, MPS, status] = fn_call_USGS_design_value_API(id, reference_doc, lat, lng, risk_cat, site_class)
% Description: Function to Call the USGG Design Values API and return,
% code short and 1 second spectra values for the base line code, risk
% targeted approach, and uniform hazard approach

% Created by: Dustin Cook
% Date Created: 2/15/2019

% Inputs:
%   id - numeric or integer unique identifier for the site
%   reference_doc - string represention of the code edition to be called
%                   options include [ asce7-05, asce7-10, asce7-16, 
%                                     nehrp-2009, nehrp-2015, nehrp-2020
%                                     asce41-13, asce41-17, and others]
%   lat - numeric lattitide of the site
%   lng - numeric longitude of the site
%   risk_cat - string representation of the risk category (I, II, III or IV)
%   site_class - sting representation of the site class (A, B, C, D, E, or F)

% Ouputs:
%   design_values = data structure of design parameters returned by call
%   (likely specific to ASCE 7-22 web format)
%   MPS = data structure of multi point spectra values returned by call
%   status = string status message to indicate if call was sucessful or not

if ~strcmp(reference_doc,'asce7-22')
    warning('Standard other than ASCE 7-22 was called. Please to make sure design value call is functioning as expected')
end

try
    %% Call USGS API
    options = weboptions('Timeout', 30);
    DATA = webread( [ 'https://earthquake.usgs.gov/ws/designmaps/' ...
                      reference_doc '.json?' ], ...
                      'latitude', lat, ...
                      'longitude', lng, ...
                      'riskCategory', risk_cat, ...
                      'siteClass', site_class, ...
                      'title', num2str(id), ...
                      options);

    %% Post Process Return
    % Conditioned on ASCE 7-22 output data structure (I think its different
    % for ASCE 7-16)
    design_values.s1 = DATA.response.data.s1;
    design_values.sm1 = DATA.response.data.sm1;
    design_values.sd1 = DATA.response.data.sd1;
    design_values.fa = DATA.response.data.sm1 / DATA.response.data.s1;
    design_values.mps = DATA.response.data.multiPeriodMCErSpectrum.ordinates(DATA.response.data.multiPeriodMCErSpectrum.periods == 1);
    design_values.ss = DATA.response.data.ss;    
    design_values.sds = DATA.response.data.sds;
    design_values.fv = DATA.response.data.sms / DATA.response.data.ss;
    
    % Set Multi Period Spectra return
    MPS = DATA.response.data.underlyingData;
    MPS.designSpectrum = DATA.response.data.multiPeriodDesignSpectrum;
    
    % Save SDC
    design_values.SDC = DATA.response.data.sdc;
    
    status = 'success';
    
    
catch ME
    design_values = [];
    MPS = [];
    status = ME.message;
end

end

