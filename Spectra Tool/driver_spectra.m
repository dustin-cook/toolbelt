clear 
close all
clc

%% Notes

%% Define Inputs Directory
gm_set_name = 'FEMA far field set';
g_factor = 1; % denominator to convert raw data to g
dt_idx = 3; % location of the dt value
eq_data_idx = 5; % location of the start of the ground motion
plot_spectra = 1;

%% Run Method
main_calc_spectra(gm_set_name, g_factor, dt_idx, eq_data_idx, plot_spectra)