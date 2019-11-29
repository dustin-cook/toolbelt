function [ spectra ] = fn_single_spectra( dt, ag, ag_2 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

spectra = table;

% Define Period Range that get more coarse as period increases
spectra.T = [0.01:0.01:1,1.02:.02:2,2.05:.05:3,3.1:.1:5]';
damp_ratio = [0.01, 0.02, 0.025, 0.03, 0.05];

%% Run Single Spectra
for i = 1:length(damp_ratio)
    for j = 1:height(spectra)
        [psuedoAccelerationTH, ~, ~, ~] = fn_sdof_th(spectra.T(j), damp_ratio(i), ag, dt);
        if exist('ag_2','var')
            [psuedoAccelerationTH_2, ~, ~] = fn_sdof_th(spectra.T(j), damp_ratio(i), ag_2, dt);
            min_dim = min([length(psuedoAccelerationTH),length(psuedoAccelerationTH_2)]);
            spectra.(['psa_' strrep(num2str(damp_ratio(i)),'.','')])(j) = max(abs(sqrt(psuedoAccelerationTH(1:min_dim).^2 + psuedoAccelerationTH_2(1:min_dim).^2))); 
        else
            spectra.(['psa_' strrep(num2str(damp_ratio(i)),'.','')])(j) = max(abs(psuedoAccelerationTH));    
        end
    end
end

end

