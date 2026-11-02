function results = run_balloon_case(params, config)
% RUN_BALLOON_CASE
% Run one Monte Carlo simulation scenario for the balloon + barometer model.
%
% The function compares four estimators:
% 1. recursive average;
% 2. causal moving average;
% 3. standard Kalman filter tracking the balloon state only;
% 4. augmented Kalman filter tracking the balloon state and barometric bias.
%
% The no-gust scenario uses a scalar random-walk model for altitude, whereas
% the gust scenario uses a constant-velocity model driven by white acceleration
% noise. In both cases, the augmented filter includes the barometric disturbance
% as an additional Gauss-Markov / AR(1) state.

dt = config.dt;
N  = round(params.T / dt) + 1;
MC = params.MC;

t = (0:N-1) * dt;

h0 = params.h0;

% Barometric correlated disturbance: GM1 / AR(1)
tau_b   = params.tau_b;
sigma_b = params.sigma_b;

phi = exp(-dt/tau_b);
Q_b = sigma_b^2 * (1 - phi^2);

% White measurement noise
sigma_v = params.sigma_v;
R    = sigma_v^2;
R_kf = R + sigma_b^2;

use_gusts = config.use_gusts;

% *********************
% True process model
% *********************

if use_gusts

    q_h = params.q_h;

    Q_cv = q_h * [dt^3/3  dt^2/2;
        dt^2/2  dt];

    A_true = [1 dt 0;
        0  1 0;
        0  0 phi];

    H_true = [1 0 1];

    Q_true = [Q_cv(1,1) Q_cv(1,2) 0;
        Q_cv(2,1) Q_cv(2,2) 0;
        0         0         Q_b];

    nx_true = 3;

else

    A_true = [1 0;
        0 phi];

    H_true = [1 1];

    % no-gust true model
    Q_h_rw = params.Q_h_rw*dt;

    Q_true = [Q_h_rw 0;
              0      Q_b];

    nx_true = 2;

end

[L_true, p] = chol(Q_true,'lower');

if p ~= 0
    [V,D]  = eig(Q_true);
    L_true = V*sqrt(max(D,0));
end

K_kf_all  = zeros(N, MC);
K_aug_all = zeros(N, MC);

% ***************************************
% Standard KF model: no barometric state
% ***************************************

if use_gusts

    q_h = params.q_h;

    Q_cv = q_h * [dt^3/3  dt^2/2;
        dt^2/2  dt];

    A_kf = [1 dt;
        0  1];

    H_kf = [1 0];

    Q_kf = Q_cv;

    x0_kf = [h0; 0];
    P0_kf = diag([10^2, 1^2]);

else

    A_kf = 1;
    H_kf = 1;
    Q_kf = params.Q_h_rw*dt;
    x0_kf = h0;
    P0_kf = 10^2;

end

nx_kf = size(A_kf, 1);

% *********************************************
% Augmented KF model: includes barometric state
% *********************************************

A_aug = A_true;
H_aug = H_true;
Q_aug = Q_true;

nx_aug = nx_true;
nz = 1;

% ********
% Storage
% ********

err_avg = zeros(N, MC);
err_kf  = zeros(N, MC);
err_aug = zeros(N, MC);
err_mov = zeros(N, MC);

NIS_kf   = zeros(N, MC);
NEES_kf  = zeros(N, MC);
NIS_aug  = zeros(N, MC);
NEES_aug = zeros(N, MC);

x_true_one = [];
z_one      = [];
avg_one    = [];
mov_one    = [];
kf_one     = [];
aug_one    = [];

% ****************
% Monte Carlo loop
% ****************

for m = 1:MC

    if rem(m, 50) == 0
        fprintf('  MC %d / %d\n', m, MC);
    end

    x_true = zeros(nx_true, N);
    z = zeros(1, N);

    if use_gusts
        x_true(:, 1) = [h0; 0; sigma_b*randn];
    else
        x_true(:, 1) = [h0; sigma_b*randn];
    end

    % Simulate true trajectory
    for k = 2:N
        w = L_true*randn(nx_true,1);
        x_true(:,k) = A_true*x_true(:,k-1) + w;
    end

    % Simulate measurements
    for k = 1:N
        z(k) = H_true*x_true(:,k) + sigma_v*randn;
    end

    % *******************************
    % Estimator 1: recursive average
    % *******************************

    zbar = zeros(1, N);
    zbar(1) = z(1);

    for k = 2:N
        zbar(k) = zbar(k-1) + (z(k) - zbar(k-1))/k;
    end

    err_avg(:,m) = zbar(:) - x_true(1,:).';

    % ***************************
    % Estimator 2: moving average
    % ***************************

    W = max(1, round(params.T_ma / dt));
    zmov = movmean(z, [W-1 0]);   % causal moving average

    err_mov(:,m) = zmov(:) - x_true(1,:).';

    % ************************
    % Estimator 3: standard KF
    % ************************

    x_hat_kf = zeros(nx_kf, N);

    x_pred = x0_kf;
    P_pred = P0_kf;

    for k = 1:N

        innov = z(k) - H_kf*x_pred;
        S = H_kf*P_pred*H_kf' + R_kf;

        K = P_pred*H_kf'/S;
        K_kf_all(k,m) = norm(K);

        x_upd = x_pred + K*innov;

        I = eye(nx_kf);
        P_upd = (I - K*H_kf)*P_pred*(I - K*H_kf)' + K*R_kf*K';

        % Diagnostics for standard KF
        if use_gusts
            e_kf = x_true(1:2,k) - x_upd;
        else
            e_kf = x_true(1,k) - x_upd;
        end

        NIS_kf(k,m)  = innov'*(S\innov);
        NEES_kf(k,m) = e_kf'*(P_upd\e_kf);

        x_hat_kf(:,k) = x_upd;

        x_pred = A_kf*x_upd;
        P_pred = A_kf*P_upd*A_kf' + Q_kf;

    end

    err_kf(:,m) = x_hat_kf(1,:).' - x_true(1,:).';

    % *************************
    % Estimator 4: augmented KF
    % *************************

    x_hat_aug = zeros(nx_aug, N);

    if use_gusts
        x_pred = [z(1); 0; 0];
        P_pred = diag([10^2, 1^2, sigma_b^2]);
    else
        x_pred = [z(1); 0];
        P_pred = diag([10^2, sigma_b^2]);
    end

    for k = 1:N

        innov = z(k) - H_aug*x_pred;
        S = H_aug*P_pred*H_aug' + R;

        K = P_pred*H_aug'/S;
        K_aug_all(k,m) = norm(K);

        x_upd = x_pred + K*innov;

        I = eye(nx_aug);
        P_upd = (I - K*H_aug)*P_pred*(I - K*H_aug)' + K*R*K';

        x_hat_aug(:,k) = x_upd;

        e = x_true(:,k) - x_upd;

        NIS_aug(k,m)  = innov'*(S\innov);
        NEES_aug(k,m) = e'*(P_upd\e);

        x_pred = A_aug*x_upd;
        P_pred = A_aug*P_upd*A_aug' + Q_aug;

    end

    err_aug(:,m) = x_hat_aug(1,:).' - x_true(1,:).';

    % Store one representative run
    if m == 1
        x_true_one = x_true;
        z_one      = z;
        avg_one    = zbar;
        mov_one    = zmov;
        kf_one     = x_hat_kf;
        aug_one    = x_hat_aug;
    end
end

% ******************
% Summary statistics
% ******************

rmse_avg = sqrt(mean(err_avg.^2, 2));
rmse_mov = sqrt(mean(err_mov.^2, 2));
rmse_kf  = sqrt(mean(err_kf.^2, 2));
rmse_aug = sqrt(mean(err_aug.^2, 2));

final_abs_err_avg_mc = abs(err_avg(end,:));
final_abs_err_mov_mc = abs(err_mov(end,:));
final_abs_err_kf_mc  = abs(err_kf(end,:));
final_abs_err_aug_mc = abs(err_aug(end,:));

traj_rmse_avg_mc = sqrt(mean(err_avg.^2, 1));
traj_rmse_mov_mc = sqrt(mean(err_mov.^2, 1));
traj_rmse_kf_mc  = sqrt(mean(err_kf.^2, 1));
traj_rmse_aug_mc = sqrt(mean(err_aug.^2, 1));

mean_NIS_kf  = mean(NIS_kf, 2);
mean_NEES_kf = mean(NEES_kf, 2);

mean_NIS_aug  = mean(NIS_aug, 2);
mean_NEES_aug = mean(NEES_aug, 2);

T_trans   = 20; % seconds
transient = find(t >= T_trans);

% Expected values for consistency diagnostics
expected_NIS = nz;

expected_NEES_kf  = nx_kf;
expected_NEES_aug = nx_aug;

% ***************
% Output package
% ***************

results.name = config.name;
results.dt = dt;
results.N = N;
results.MC = MC;
results.t = t;

results.use_gusts = use_gusts;
results.phi = phi;
results.Q_b = Q_b;
results.R = R;

results.K_kf_mean  = mean(K_kf_all, 2);
results.K_aug_mean = mean(K_aug_all, 2);

results.rmse_avg = rmse_avg;
results.rmse_mov = rmse_mov;
results.rmse_kf  = rmse_kf;
results.rmse_aug = rmse_aug;

results.rmse_avg_final = rmse_avg(end);
results.rmse_mov_final = rmse_mov(end);
results.rmse_kf_final  = rmse_kf(end);
results.rmse_aug_final = rmse_aug(end);

results.final_abs_err_avg_mc = final_abs_err_avg_mc;
results.final_abs_err_mov_mc = final_abs_err_mov_mc;
results.final_abs_err_kf_mc  = final_abs_err_kf_mc;
results.final_abs_err_aug_mc = final_abs_err_aug_mc;

results.traj_rmse_avg_mc = traj_rmse_avg_mc;
results.traj_rmse_mov_mc = traj_rmse_mov_mc;
results.traj_rmse_kf_mc  = traj_rmse_kf_mc;
results.traj_rmse_aug_mc = traj_rmse_aug_mc;

results.mean_NIS_kf  = mean(mean_NIS_kf(transient));
results.mean_NEES_kf = mean(mean_NEES_kf(transient));

results.expected_NIS_kf  = expected_NIS;
results.expected_NEES_kf = expected_NEES_kf;

results.mean_NIS_aug  = mean(mean_NIS_aug(transient));
results.mean_NEES_aug = mean(mean_NEES_aug(transient));

results.expected_NIS_aug  = expected_NIS;
results.expected_NEES_aug = expected_NEES_aug;

results.x_true_one = x_true_one;
results.z_one      = z_one;
results.avg_one    = avg_one;
results.mov_one    = mov_one;
results.kf_one     = kf_one;
results.aug_one    = aug_one;

end