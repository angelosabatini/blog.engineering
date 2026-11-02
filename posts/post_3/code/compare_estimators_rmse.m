function compare_estimators_rmse(r)
% COMPARE_ESTIMATORS_RMSE
% Paired statistical comparison of trajectory RMSE across estimators.
%
% The comparison is performed across Monte Carlo runs for one simulation
% scenario. The recursive average is omitted from the formal test because it
% serves as a naive baseline; the comparison focuses on the moving average,
% the standard Kalman filter, and the augmented Kalman filter.

X = [ ...
    r.traj_rmse_mov_mc(:), ...
    r.traj_rmse_kf_mc(:), ...
    r.traj_rmse_aug_mc(:) ];

[p_friedman, ~, ~] = friedman(X, 1, 'off');

[p_mov_kf,  ~] = signrank(X(:,1), X(:,2));
[p_mov_aug, ~] = signrank(X(:,1), X(:,3));
[p_kf_aug,  ~] = signrank(X(:,2), X(:,3));

p_raw  = [p_mov_kf, p_mov_aug, p_kf_aug];
p_bonf = min(3*p_raw, 1);

D = [ ...
    X(:,1) - X(:,2), ...
    X(:,1) - X(:,3), ...
    X(:,2) - X(:,3) ];

fprintf('\nFriedman test, %s\n', r.name)
fprintf('p = %.3g\n', p_friedman)

fprintf('\nPairwise signed-rank tests, Bonferroni-corrected\n')
fprintf('MOV vs KF  : p = %.3g, p_corr = %.3g\n', p_mov_kf,  p_bonf(1))
fprintf('MOV vs AUG : p = %.3g, p_corr = %.3g\n', p_mov_aug, p_bonf(2))
fprintf('KF  vs AUG : p = %.3g, p_corr = %.3g\n', p_kf_aug,  p_bonf(3))

fprintf('\nMedian paired differences in trajectory RMSE [m]\n')
labels = ["MOV - KF", "MOV - AUG", "KF - AUG"];

for j = 1:3
    interval = prctile(D(:,j), [2.5 97.5]);

    fprintf('%s : %.4f, central 95%% empirical interval [%.4f, %.4f]\n', ...
        labels(j), ...
        median(D(:,j)), ...
        interval(1), ...
        interval(2))
end

end