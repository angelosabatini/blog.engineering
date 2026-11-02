%% analyze_post3_results
% Post-processing script for the balloon + barometer simulation.
%
% This script loads the saved Monte Carlo results, prints summary tables,
% produces exploratory plots, and runs paired statistical comparisons of
% trajectory RMSE across estimators.

clear
clc
close all

load post3_balloon_case.mat

results = post3.results;

%% Summary table

fprintf('\nSummary: final RMSE\n')
fprintf('-------------------------------------------------------------------------------\n')
fprintf('%24s | %8s | %8s | %8s | %8s\n', ...
    'case', 'average', 'mov avg', 'std KF', 'aug KF')
fprintf('-------------------------------------------------------------------------------\n')

for i = 1:numel(results)
    fprintf('%24s | %8.3f | %8.3f | %8.3f | %8.3f\n', ...
        results{i}.name, ...
        results{i}.rmse_avg_final, ...
        results{i}.rmse_mov_final, ...
        results{i}.rmse_kf_final, ...
        results{i}.rmse_aug_final)
end

%% Example plots

idx = 3;
r = results{idx};

figure
plot(r.t, r.K_kf_mean, 'LineWidth', 1.8)
hold on
plot(r.t, r.K_aug_mean, 'LineWidth', 1.8)
grid on
xlabel('Time [s]')
ylabel('Mean Kalman gain norm')
legend('standard KF', 'augmented KF', 'Location', 'best')
title(r.name)
set(gca, 'FontSize', 14, 'LineWidth', 1.5)

figure
plot(r.t, r.x_true_one(1,:), 'k', 'LineWidth', 2)
hold on
plot(r.t, r.z_one, '.', 'MarkerSize', 10)
plot(r.t, r.avg_one, 'LineWidth', 1.8)
plot(r.t, r.mov_one, 'LineWidth', 1.8)
plot(r.t, r.kf_one(1,:), 'LineWidth', 1.8)
plot(r.t, r.aug_one(1,:), 'LineWidth', 1.8)
grid on
xlabel('Time [s]')
ylabel('Altitude [m]')
legend('true altitude', ...
       'barometric measurements', ...
       'recursive average', ...
       'moving average', ...
       'standard KF', ...
       'augmented KF', ...
       'Location', 'best')
%title(r.name)
set(gca, 'FontSize', 14, 'LineWidth', 1.5)

figure
%plot(r.t, r.rmse_avg, 'LineWidth', 1.8)
plot(r.t, r.rmse_mov, 'LineWidth', 1.8)
hold on
plot(r.t, r.rmse_kf, 'LineWidth', 1.8)
plot(r.t, r.rmse_aug, 'LineWidth', 1.8)
grid on
xlabel('Time [s]')
ylabel('Altitude RMSE [m]')
legend('moving average', ...
       'standard KF', ...
       'augmented KF', ...
       'Location', 'best')
%title(r.name)
set(gca, 'FontSize', 14, 'LineWidth', 1.5)

%% Final RMSE comparison across all scenarios

case_names = strings(numel(results), 1);
rmse_final = zeros(numel(results), 4);

for i = 1:numel(results)
    case_names(i) = results{i}.name;

    rmse_final(i,:) = [ ...
        results{i}.rmse_avg_final, ...
        results{i}.rmse_mov_final, ...
        results{i}.rmse_kf_final, ...
        results{i}.rmse_aug_final ];
end

figure
bar(rmse_final)
grid on
ylabel('Final altitude RMSE [m]')
set(gca, ...
    'XTickLabel', case_names, ...
    'XTickLabelRotation', 30)
legend('recursive average', ...
       'moving average', ...
       'standard KF', ...
       'augmented KF', ...
       'Location', 'best')
set(gca, 'FontSize', 12, 'LineWidth', 1.5)

fprintf('\nFinal RMSE comparison\n')
fprintf('---------------------------------------------------------------------------------------------------\n')
fprintf('%24s | %12s | %12s | %12s | %12s\n', ...
    'case', 'average', 'mov avg', 'std KF', 'aug KF')
fprintf('---------------------------------------------------------------------------------------------------\n')

for i = 1:numel(results)
    fprintf('%24s | %12.3f | %12.3f | %12.3f | %12.3f\n', ...
        results{i}.name, ...
        results{i}.rmse_avg_final, ...
        results{i}.rmse_mov_final, ...
        results{i}.rmse_kf_final, ...
        results{i}.rmse_aug_final)
end

%% Consistency diagnostics

fprintf('\nConsistency diagnostics\n')
fprintf('-----------------------------------------------------------------------------------------------------------------------------\n')
fprintf('%24s | %10s | %10s | %10s | %10s\n', ...
    'case', 'KF NIS', 'KF NEES', 'AUG NIS', 'AUG NEES')
fprintf('-----------------------------------------------------------------------------------------------------------------------------\n')

for i = 1:numel(results)
    fprintf('%24s | %10.3f | %10.3f | %10.3f | %10.3f\n', ...
        results{i}.name, ...
        results{i}.mean_NIS_kf, ...
        results{i}.mean_NEES_kf, ...
        results{i}.mean_NIS_aug, ...
        results{i}.mean_NEES_aug)
end

%% Statistical analysis

for i = 1:numel(results)
    compare_estimators_rmse(results{i});
end