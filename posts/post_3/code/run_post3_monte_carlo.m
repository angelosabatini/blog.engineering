%% run_post3_monte_carlo
% Monte Carlo driver for the balloon + barometer simulation.
%
% This script defines the common simulation parameters, runs the four
% scenarios considered in the post, and saves the resulting structure to
% disk. No plotting or statistical analysis is performed here.

clear
clc
close all

rng(2028)

%% Common parameters

params.T  = 80;           % Physical duration of the tracking task [s]
params.MC = 250;          % Monte Carlo runs

params.h0 = 250;          % Nominal altitude [m]

% Barometric correlated disturbance: GM1 / AR(1)
params.tau_b   = 5;       % Correlation time [s]
params.sigma_b = 0.5;     % Stationary standard deviation [m]

% White measurement noise
params.sigma_v = 0.1;     % Standard deviation [m]

% Moving average
params.T_ma = 5;          % Causal moving-average window [s]

% Altitude random-walk model, no-gust scenarios
params.sigma_h_rw = 0.01; % Random-walk standard deviation per sample [m]
params.Q_h_rw = params.sigma_h_rw^2;

% Vertical gust model, CV scenarios
params.q_h = 1e-4;        % Acceleration noise intensity [m^2/s^3]

%% Four experiments

configs = struct([]);

configs(1).name = "no gusts, dt = 1 s";
configs(1).dt = 1.0;
configs(1).use_gusts = false;

configs(2).name = "no gusts, dt = 0.1 s";
configs(2).dt = 0.1;
configs(2).use_gusts = false;

configs(3).name = "gusts, dt = 1 s";
configs(3).dt = 1.0;
configs(3).use_gusts = true;

configs(4).name = "gusts, dt = 0.1 s";
configs(4).dt = 0.1;
configs(4).use_gusts = true;

%% Run Monte Carlo simulations

results = cell(numel(configs), 1);

for i = 1:numel(configs)
    fprintf('\nRunning case: %s\n', configs(i).name)
    results{i} = run_balloon_case(params, configs(i));
end

post3.params  = params;
post3.configs = configs;
post3.results = results;

save('post3_balloon_case.mat', 'post3')