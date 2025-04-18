%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         January 22, 2024
% CLASS:        ASEN 6080: StatOD
% INSTRUCTOR:   Prof. Jay W. Mcmahon
% ASSIGNMENT:   HW8 - UQ Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close;

%% Create Reference Trajectory

% Constants
coeffs = [0, 1.082626925638815 * 1e-3]; % Zonal harmonic coefficients (i.e., [J1, J2])
R = 6378.1363; % Earth radius 
mu = 398600.4415;  % Gravitational parameter 
l_max = 2; % Zonal Harmonics Terms 
A = 3 * 1e-6; % Cross-Sectional Area
m = 970; % Mass
CD = 2.0; % Drag Coefficient
param_meas = 9;

% Initial state (reference orbit)
state0 = [757.7; 5222.607; 4851.5; 2.21321; 4.67834; -5.3713]; % (km, km/s)
state0 = [state0; reshape(eye(6), [], 1)]; % Adding STM

% Initialize ODEs
[acc_SPH, A_SPH] = precomputeZonalSPH(R, l_max);
[acc_Drag, A_Drag] = precomputeDrag(A, m, R);

% Simulation parameters
T_orbit = 24 * 60 * 60; % Orbital period (s)
t_span = 0:T_orbit;      
options = odeset('RelTol', 2.22045 * 1e-14, 'AbsTol', 2.22045 * 1e-16);

% Integrate reference trajectory
[t_ref, state_ref] = ode45(@(t, state) ...
     sc_ODE(t, state, coeffs, mu, l_max, CD, acc_SPH, A_SPH, acc_Drag, ...
     A_Drag, param_meas), t_span, state0, options);

% Plot Orbit
gca1 = figure(1);
hold on;  
plot3(state_ref(:, 1), state_ref(:, 2), state_ref(:, 3), ...
    'r.', 'MarkerSize', 10);   
plot_earth();  
axis equal;
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('z (km)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Reference Trajectory at $t_{obs}$', 'Earth', ...
    'Location', 'northwest', 'Interpreter', 'latex','FontSize', 14);
grid on;
view(3);
exportgraphics(gca1, 'OrtbitUQ.pdf', 'ContentType','vector');


%% Initialization UQ Methods

% Initial uncertainty (covariance)
sigma_pos = 1;       % km
sigma_vel = 1e-3;    % km/s
P0 = diag([sigma_pos^2 * ones(1,3), sigma_vel^2 * ones(1,3)]);

% Extract nominal state and STM
x0 = state0(1:6);
STM0 = eye(6);

% Time tags for comparison (in seconds)
t_eval_hours = [6, 12, 18, 24];
t_eval = t_eval_hours * 3600;


%% Linear Covariance Propagation

% Preallocate
lin_states = zeros(6, length(t_eval));
lin_covs = zeros(6, 6, length(t_eval));

fprintf('\n===== Propagating Reference Trajectory for LinCov =====\n');

% Propagate nominal trajectory and STM directly at t_eval times
[~, state_nom_all] = ode113(@(t, state) ...
    sc_ODE(t, state, coeffs, mu, l_max, CD, acc_SPH, A_SPH, ...
           acc_Drag, A_Drag, param_meas), ...
    t_eval, [x0; reshape(STM0, 36, 1)], options);

% Now we know state_nom_all is already sampled at t_eval

fprintf('Extracting nominal states and propagating covariance...\n');
for i = 1:length(t_eval)
    x_nom = state_nom_all(i, 1:6)';
    STM_t = reshape(state_nom_all(i, 7:end), 6, 6);
    lin_states(:, i) = x_nom;
    lin_covs(:, :, i) = STM_t * P0 * STM_t';
    fprintf('  → t = %2d hr...\n', t_eval_hours(i));
end

fprintf('===== LinCov propagation complete =====\n');

%% Monte Carlo Propagation 

% Number of Monte Carlo samples
N_MC = 1000;

fprintf('\n===== Monte Carlo Propagation (%d samples, %d epochs) =====\n', N_MC, length(t_eval));

% Preallocate
MC_states = zeros(6, N_MC, length(t_eval));

% Sample Monte Carlo ICs
rng(42);  % for reproducibility
X_samples = mvnrnd(x0, P0, N_MC)';

% Loop over each sample
for j = 1:N_MC
    if mod(j, 10) == 0 || j == 1
        fprintf('[%3d%%] Propagating sample %4d / %d\n', round(100*j/N_MC), j, N_MC);
    end

    x_mc = X_samples(:, j);

    % Direct propagation at desired epochs
    [~, state_mc_all] = ode113(@(t, state) ...
        sc_ODE(t, state, coeffs, mu, l_max, CD, acc_SPH, A_SPH, ...
               acc_Drag, A_Drag, param_meas), ...
        t_eval, [x_mc; reshape(STM0, 36, 1)], options);

    for i = 1:length(t_eval)
        MC_states(:, j, i) = state_mc_all(i, 1:6)';
    end
end

% Allocate mean/cov storage
MC_mean = zeros(6, length(t_eval));
MC_cov  = zeros(6, 6, length(t_eval));

fprintf('\n===== Computing MC Mean and Covariance at Each Evaluation Time =====\n');
for i = 1:length(t_eval)
    fprintf('  → t = %2d hr...\n', t_eval_hours(i));
    X = squeeze(MC_states(:, :, i));  % size: 6 x N_MC
    MC_mean(:, i) = mean(X, 2);
    dx = X - mean(X, 2);
    MC_cov(:, :, i) = (dx * dx.') / (N_MC - 1);
end

fprintf('===== Monte Carlo Complete =====\n');

%% Unscented Transform Propagation 

fprintf('\n===== Unscented Transform Propagation =====\n');

% UT parameters
n = length(x0);
alpha = 1; beta = 2; kappa = 3 - n;
lambda = alpha^2 * (n + kappa) - n;
gamma = sqrt(n + lambda);

% Generate sigma points
S = sqrtm(P0);
X_ut = [x0, x0 + gamma * S, x0 - gamma * S];  % shape: 6 x 13

% Weights
w_m = [lambda / (n + lambda), repmat(1 / (2 * (n + lambda)), 1, 2*n)];
w_c = w_m;
w_c(1) = w_c(1) + (1 - alpha^2 + beta);

% Preallocate
ut_states = zeros(6, length(t_eval));
ut_covs = zeros(6, 6, length(t_eval));
ut_trajectories = zeros(6, size(X_ut, 2), length(t_eval));

% Propagate each sigma point
for j = 1:size(X_ut, 2)
    fprintf('[%2d/%d] Propagating sigma point...\n', j, size(X_ut, 2));
    [~, state_ut] = ode113(@(t, state) ...
        sc_ODE(t, state, coeffs, mu, l_max, CD, acc_SPH, A_SPH, ...
               acc_Drag, A_Drag, param_meas), ...
        t_eval, [X_ut(:, j); reshape(STM0, 36, 1)], options);

    for i = 1:length(t_eval)
        ut_trajectories(:, j, i) = state_ut(i, 1:6)';
    end
end

% Compute UT mean and covariance at each epoch
fprintf('\n===== Computing UT Mean and Covariance =====\n');
for i = 1:length(t_eval)
    fprintf('  → t = %2d hr\n', t_eval_hours(i));
    X = squeeze(ut_trajectories(:, :, i));  % shape: 6 x 13
    mu_ut = X * w_m';                       % mean
    dx = X - mu_ut;
    P_ut = zeros(n);
    for j = 1:size(X_ut, 2)
        P_ut = P_ut + w_c(j) * (dx(:, j) * dx(:, j)');
    end
    ut_states(:, i) = mu_ut;
    ut_covs(:, :, i) = P_ut;
end

fprintf('===== Unscented Transform Complete =====\n');

%% Gaussian Mixture Model Propagation (K = 5)

% K = 100;
% fprintf('\n===== GMM Propagation with %d Components =====\n', K);
% 
% % Sample K Gaussian components from initial distribution
% mu_gmm_ic = mvnrnd(x0, P0, K);   % size: K x 6
% w_gmm = ones(K, 1) / K;       % uniform weights
% 
% % Preallocate
% gmm_states = zeros(6, K, length(t_eval));
% gmm_mean   = zeros(6, length(t_eval));
% gmm_cov    = zeros(6, 6, length(t_eval));
% 
% % Propagate each Gaussian component
% for k = 1:K
%     fprintf('  → Propagating GMM component %d / %d\n', k, K);
%     [~, state_gmm] = ode113(@(t, state) ...
%         sc_ODE(t, state, coeffs, mu, l_max, CD, acc_SPH, A_SPH, ...
%                acc_Drag, A_Drag, param_meas), ...
%         t_eval, [mu_gmm_ic(k, :)'; reshape(STM0, 36, 1)], options);
% 
%     for i = 1:length(t_eval)
%         gmm_states(:, k, i) = state_gmm(i, 1:6)';
%     end
% end
% 
% % Compute weighted mean and covariance at each evaluation time
% fprintf('\n===== Computing GMM Mean and Covariance =====\n');
% for i = 1:length(t_eval)
%     fprintf('  → t = %2d hr\n', t_eval_hours(i));
%     Xk = squeeze(gmm_states(:, :, i));   % shape: 6 x K
%     mu_i = Xk * w_gmm;                   % weighted mean
%     dx = Xk - mu_i;
%     P_i = zeros(6, 6);
%     for k = 1:K
%         P_i = P_i + w_gmm(k) * (dx(:, k) * dx(:, k)');
%     end
% 
%     gmm_mean(:, i) = mu_i;
%     gmm_cov(:, :, i) = P_i;
% end
% 
% fprintf('===== GMM Propagation Complete =====\n');

%% Gaussian Mixture Model Propagation (Fixed version)

K = 100;
fprintf('\n===== Corrected GMM Propagation with %d Components =====\n', K);

% Step 1: Sample GMM means from initial distribution
mu_gmm_ic = mvnrnd(x0, P0, K)';  % 6 x K

% Step 2: Compute weights based on initial PDF
pdf_vals = mvnpdf(mu_gmm_ic', x0', P0);   % size K x 1
w_gmm = pdf_vals / sum(pdf_vals);         % Normalize

% Step 3: Assign same covariance to each component (could scale P0 if desired)
P_gmm_ic = repmat(P0, 1, 1, K);           % 6 x 6 x K
 
% % Draw N samples from the initial Gaussian
% N_samples = 1000;
% X_samples = mvnrnd(x0, P0, N_samples);  % N x 6
% 
% % Fit GMM and choose best number of components based on BIC
% max_K = 10;
% bic_vals = zeros(1, max_K);
% gmm_models = cell(1, max_K);
% 
% fprintf('\n===== Fitting GMMs and Selecting Optimal Number of Components =====\n');
% for k = 1:max_K
%     fprintf('  → Trying GMM with %d components...\n', k);
%     gmm_models{k} = fitgmdist(X_samples, k, ...
%         'CovarianceType','full', 'RegularizationValue',1e-6, ...
%         'Options', statset('MaxIter',500, 'Display','off'));
%     bic_vals(k) = gmm_models{k}.BIC;
% end
% 
% % Select best model (lowest BIC)
% [~, best_K] = min(bic_vals);
% gmm = gmm_models{best_K};
% fprintf('Best number of components based on BIC: %d\n', best_K);
% K = best_K;
% 
% % Extract GMM parameters
% w_gmm     = gmm.ComponentProportion;              % 1 x K
% mu_gmm_ic = gmm.mu';                              % 6 x K
% P_gmm_ic  = reshape(gmm.Sigma, 6, 6, best_K);     % 6 x 6 x K

% Preallocate storage
gmm_states = zeros(6, K, length(t_eval));
gmm_covs   = zeros(6, 6, K, length(t_eval));
gmm_mean   = zeros(6, length(t_eval));
gmm_cov    = zeros(6, 6, length(t_eval));

% Step 4: Propagate each component linearly
for k = 1:K
    fprintf('  → Propagating GMM component %d / %d\n', k, K);

    % Flatten STM for integration
    IC = [mu_gmm_ic(:, k); reshape(eye(6), 36, 1)];

    [~, traj_k] = ode113(@(t, state) ...
        sc_ODE(t, state, coeffs, mu, l_max, CD, acc_SPH, A_SPH, ...
               acc_Drag, A_Drag, param_meas), ...
        t_eval, IC, options);

    for i = 1:length(t_eval)
        xk = traj_k(i, 1:6)';
        Phi = reshape(traj_k(i, 7:end), 6, 6);
        Pk = Phi * P_gmm_ic(:,:,k) * Phi';

        gmm_states(:, k, i) = xk;
        gmm_covs(:, :, k, i) = Pk;
    end
end

% Step 5: Fuse GMM statistics at each time step
fprintf('\n===== Computing GMM Mixture Statistics =====\n');
for i = 1:length(t_eval)
    fprintf('  → t = %2d hr\n', t_eval_hours(i));

    mu_i = zeros(6, 1);
    P_i = zeros(6, 6);
    
    for k = 1:K
        mu_k = gmm_states(:, k, i);
        P_k  = gmm_covs(:, :, k, i);
        mu_i = mu_i + w_gmm(k) * mu_k;
    end

    for k = 1:K
        mu_k = gmm_states(:, k, i);
        P_k  = gmm_covs(:, :, k, i);
        dmu  = mu_k - mu_i;
        P_i = P_i + w_gmm(k) * (P_k + dmu * dmu');
    end

    gmm_mean(:, i) = mu_i;
    gmm_cov(:, :, i) = P_i;
end

fprintf('===== GMM Propagation Complete =====\n');


%% Corner Plots for Monte Carlo States at Multiple Time Steps

% Initialize
n = 6;
t_indices = 1:4;
labels = {'$x$ [km]', '$y$ [km]', '$z$ [km]', ...
          '$v_x$ [km/s]', '$v_y$ [km/s]', '$v_z$ [km/s]'};

for t_idx = t_indices
    % Extract Monte Carlo samples at current time step
    X_mc = MC_states(:, :, t_idx);  % 6 x 1000
    mu_mc = mean(X_mc, 2);          % 6x1
    C_mc = cov(X_mc');              % 6x6

    % Extract other methods’ estimates
    mu_lin = lin_states(:, t_idx);
    C_lin  = lin_covs(:, :, t_idx);

    mu_ut  = ut_states(:, t_idx);
    C_ut   = ut_covs(:, :, t_idx);

    mu_gmm = gmm_mean(:, t_idx);
    C_gmm  = gmm_cov(:, :, t_idx);

    % Create figure
    figure('Color','w','Position',[100 100 1200 1000]);
    t = tiledlayout(n,n,'Padding','compact','TileSpacing','compact');
    set(gcf, 'DefaultAxesFontSize', 10)
    title(t, sprintf('Corner Plot - Time Index %d', t_idx), 'FontSize', 14, 'Interpreter', 'latex');
    legend_shown = false;
    for i = 1:n
        for j = 1:i
            nexttile((i-1)*n + j); hold on;

            if i == j
                % Histogram of MC samples
                histogram(X_mc(i,:), 50, 'Normalization','pdf', ...
                    'FaceColor',[0.3 0.6 1], 'DisplayName','MC');

                % Plot Gaussian PDFs
                x_vals = linspace(min(X_mc(i,:)), max(X_mc(i,:)), 200);
                plot(x_vals, normpdf(x_vals, mu_lin(i), sqrt(C_lin(i,i))), 'r-', 'LineWidth', 1.2, 'DisplayName','LinCov');
                plot(x_vals, normpdf(x_vals, mu_ut(i), sqrt(C_ut(i,i))), 'g-', 'LineWidth', 1.2, 'DisplayName','UT');
                plot(x_vals, normpdf(x_vals, mu_gmm(i), sqrt(C_gmm(i,i))), 'm-', 'LineWidth', 1.2, 'DisplayName','GMM');

                title(labels{i}, 'Interpreter', 'latex');
            else
                scatter(X_mc(j,:), X_mc(i,:), 5, 'b', 'filled', ...
                        'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);

                % Ellipses (only add labels once to avoid legend spam)
                if ~legend_shown
                    plot_cov_ellipse(C_mc([j i],[j i]), mu_mc([j i]), 2, 'k', 1.5, 'MC');
                    plot_cov_ellipse(C_lin([j i],[j i]), mu_lin([j i]), 2, 'r', 1.2, 'LinCov');
                    plot_cov_ellipse(C_ut([j i],[j i]), mu_ut([j i]), 2, 'g', 1.2, 'UT');
                    plot_cov_ellipse(C_gmm([j i],[j i]), mu_gmm([j i]), 2, 'm', 1.2, 'GMM');
                    legend_shown = true;
                else
                    plot_cov_ellipse(C_mc([j i],[j i]), mu_mc([j i]), 2, 'k', 1.5);
                    plot_cov_ellipse(C_lin([j i],[j i]), mu_lin([j i]), 2, 'r', 1.2);
                    plot_cov_ellipse(C_ut([j i],[j i]), mu_ut([j i]), 2, 'g', 1.2);
                    plot_cov_ellipse(C_gmm([j i],[j i]), mu_gmm([j i]), 2, 'm', 1.2);
                end

                if i == n, xlabel(labels{j}, 'Interpreter', 'latex'); end
                if j == 1, ylabel(labels{i}, 'Interpreter', 'latex'); end
            end
            axis tight; grid on;
        end
    end

    % Single legend
    lgd = legend({'MC','LinCov','UT','GMM'}, ...
        'Location','southoutside','Orientation','horizontal', ...
        'Interpreter','latex','FontSize',12);
    lgd.Layout.Tile = 'south';

    % Save figure as PDF
    filename = sprintf('corner_plot_t%02d.pdf', t_idx);
    exportgraphics(gcf, filename, 'ContentType','vector');
end

% Function: Plot Covariance Ellipse
function plot_cov_ellipse(C, mu, nsig, color, lw, label)
    theta = linspace(0, 2*pi, 100);
    unit_circle = [cos(theta); sin(theta)];
    [V, D] = eig(C);
    ellipse = V * sqrt(D) * unit_circle * nsig;
    if nargin < 6
        plot(mu(1) + ellipse(1,:), mu(2) + ellipse(2,:), ...
            'Color', color, 'LineWidth', lw);
    else
        plot(mu(1) + ellipse(1,:), mu(2) + ellipse(2,:), ...
            'Color', color, 'LineWidth', lw, 'DisplayName', label);
    end
end

%% Print Stats

fprintf('\n===== Mean & Covariance Comparison vs Monte Carlo =====\n');

for t_idx = t_indices
    fprintf('\n--- Time Index %d ---\n', t_idx);

    % Monte Carlo
    X_mc = MC_states(:, :, t_idx);  % 6 x N
    mu_mc = mean(X_mc, 2);          % 6x1
    C_mc = cov(X_mc');              % 6x6

    % Other method estimates
    method_names = {'LinCov', 'UT', 'GMM'};
    mus = [lin_states(:, t_idx), ut_states(:, t_idx), gmm_mean(:, t_idx)];
    covs = cat(3, lin_covs(:, :, t_idx), ut_covs(:, :, t_idx), gmm_cov(:, :, t_idx));

    for m = 1:length(method_names)
        mu_i = mus(:, m);
        C_i = covs(:, :, m);

        % Compute normed errors
        mean_err = norm(mu_i - mu_mc);
        cov_err  = norm(C_i - C_mc, 'fro');

        % Optionally compute % error
        mean_err_pct = 100 * mean_err / norm(mu_mc);
        cov_err_pct  = 100 * cov_err  / norm(C_mc, 'fro');

        % KL divergence D_KL(P_mc || P_i)
        delta_mu = mu_mc - mu_i;
        inv_Ci = inv(C_i);
        tr_term = trace(inv_Ci * C_mc);
        mahal_term = delta_mu' * inv_Ci * delta_mu;
        logdet_term = log(det(C_i) / det(C_mc));
        D_KL = 0.5 * (tr_term + mahal_term - n + logdet_term);

        % Print with high precision
        fprintf('\n[%s]\n', method_names{m});
        fprintf('  Mean Error       : %.14e\n', mean_err);
        fprintf('  Covariance Error : %.14e\n', cov_err);
        fprintf('  Mean Error %%      : %.10f %%\n', mean_err_pct);
        fprintf('  Cov Error %%       : %.10f %%\n', cov_err_pct);
        fprintf('  KL Divergence     : %.10f\n', D_KL);
    end
end

%% Print mean and std for all methods at each evaluation time

fprintf('\n===== Method-wise μ and σ Summary (at each t_eval) =====\n');

method_names = {'MC', 'LinCov', 'UT', 'GMM'};

for i = 1:length(t_eval)
    fprintf('\n--- Time %2d hr ---\n', t_eval_hours(i));
    
    % Monte Carlo stats
    X_mc = squeeze(MC_states(:, :, i));     % 6 x N_MC
    mu_MC = mean(X_mc, 2);                  % 6 x 1
    sigma_MC = std(X_mc, 0, 2);             % 6 x 1

    % LinCov (mean is nominal, std is from diag of covariance)
    mu_Lin = lin_states(:, i);
    sigma_Lin = sqrt(diag(lin_covs(:, :, i)));

    % UT
    mu_UT = ut_states(:, i);
    sigma_UT = sqrt(diag(ut_covs(:, :, i)));

    % GMM
    mu_GMM = gmm_mean(:, i);
    sigma_GMM = sqrt(diag(gmm_cov(:, :, i)));

    % Combine into a table: 6 x 4 for each of mean and std
    mu_all = [mu_MC, mu_Lin, mu_UT, mu_GMM];
    sigma_all = [sigma_MC, sigma_Lin, sigma_UT, sigma_GMM];

    for j = 1:6
        fprintf('\n[%s]\n', labels{j});
        fprintf('  %-7s  μ = %+ .6e   σ = %.6e\n', method_names{1}, mu_all(j,1), sigma_all(j,1));
        fprintf('  %-7s  μ = %+ .6e   σ = %.6e\n', method_names{2}, mu_all(j,2), sigma_all(j,2));
        fprintf('  %-7s  μ = %+ .6e   σ = %.6e\n', method_names{3}, mu_all(j,3), sigma_all(j,3));
        fprintf('  %-7s  μ = %+ .6e   σ = %.6e\n', method_names{4}, mu_all(j,4), sigma_all(j,4));
    end
end


%% Containment within 2σ Ellipsoids

methods = {'LinCov', 'UT', 'GMM'};
mu_all   = {lin_states, ut_states, gmm_mean};
cov_all  = {lin_covs,   ut_covs,   gmm_cov};

fprintf('\n===== Percentage of Monte Carlo Samples within 2σ Ellipsoids =====\n');

for i = 1:length(t_eval)
    fprintf('\n--- Time %2d hr ---\n', t_eval_hours(i));
    X = squeeze(MC_states(:, :, i));  % 6 x N_MC

    for m = 1:length(methods)
        mu_i = mu_all{m}(:, i);
        P_i  = cov_all{m}(:, :, i);

        % Mahalanobis distances for all MC samples
        dx = X - mu_i;
        mahal_sq = sum((P_i \ dx) .* dx, 1);  % 1 x N_MC

        % Count how many are within 2σ
        inside = mahal_sq <= chi2inv(0.954, 6);  % 95.4% if Gaussian
        pct_inside = 100 * sum(inside) / N_MC;

        fprintf('[%s] %% MC inside 2σ ellipsoid: %.4f %%\n', methods{m}, pct_inside);
    end
end


%% Helper Functions

function R = DCM_eciToEcef(t)
    % Earth's rotation rate 
    % omega_earth = 2 * pi / (24 * 3600);
    omega_earth = 7.2921158553 * 1e-5; % OSS: Given my assignment.
    
    % Rotation angle
    theta = omega_earth * t;
    
    % Rotation matrix for ECI to ECEF
    R = [cos(theta), sin(theta), 0;
        - sin(theta), cos(theta), 0;
         0,          0,          1];
end

function [acceleration_func, A_func] = precomputeZonalSPH(R_val, l_max)
   % Symbolic definitions
    syms x y z vx vy vz mu_sym real
    coeffs = sym('J', [l_max, 1], 'real');       

    % Compute spherical coordinates
    r = sqrt(x^2 + y^2 + z^2);
    phi = asin(z / r);                              % Latitude

    % Gravitational potential
    U = mu_sym / r;
    if l_max > 0
        % Add zonal harmonic contributions
        for l = 1:l_max
            U = U - mu_sym / r * coeffs(l) * (R_val / r)^l * legendreP(l, sin(phi));
        end
    end

    % Compute symbolic acceleration in Cartesian coordinates
    acc_state = gradient(U, [x, y, z]);

    % Define augmented acceleration including parameter derivatives
    zero_acc_mu = sym(0);                           % Placeholder for mu derivative
    zero_acc_Jl = sym(zeros(l_max, 1));             % Placeholder for coefficients
    augmented_acceleration = [acc_state; zero_acc_mu; zero_acc_Jl];

    % Convert to a numeric function
    acceleration_func = matlabFunction(augmented_acceleration, 'Vars', {x, y, z, mu_sym, coeffs});

    % Include velocity for Jacobian computation
    augmented_acceleration_full = [vx; vy; vz; augmented_acceleration];

    % Compute full Jacobian of augmented acceleration
    A_sym = jacobian(augmented_acceleration_full, [x, y, z, vx, vy, vz, mu_sym, coeffs']);

    % Convert to a numeric function
    A_func = matlabFunction(A_sym, 'Vars', {x, y, z, vx, vy, vz, mu_sym, coeffs});
end

function [acceleration_func, A_func] = precomputeDrag(A, m, REarth)
    % Constants
    omega_E = 7.2921158553 * 1e-5; % Earth's angular velocity (rad/s)
    rho0 = 3.614 * 1e-4; % Reference Density(kg/km^3)
    r0 = (700.00 + REarth); % Reference Radius (km)
    H = 88.6670; % Scale Height(km)

    % Symbolic definitions
    syms x y z vx vy vz CD real

    % Define position and velocity vectors
    r_vec = [x; y; z];
    v_vec = [vx; vy; vz];

    % Earth's rotation vector (assumed along z-axis)
    omega_vec = [0; 0; omega_E];

    % Compute velocity relative to atmosphere
    v_rel = v_vec - cross(omega_vec, r_vec);
    v_rel_mag = norm(v_rel);

    % Compute atmospheric density at current altitude
    r = norm(r_vec);
    rho = rho0 * exp(-(r - r0) / H);

    % Drag acceleration using relative velocity
    a_drag = -0.5 * (CD * A / m) * rho * v_rel_mag * v_rel;
    zero_acc_CD = sym(0); % Placeholder for CD derivative
    augmented_acceleration = [a_drag; zero_acc_CD];

    % Convert to a numeric function
    acceleration_func = matlabFunction(augmented_acceleration, ...
        'Vars', {x, y, z, vx, vy, vz, CD});

    % Compute the augmented acceleration and Jacobian of the system
    augmented_acceleration_full = [vx; vy; vz; augmented_acceleration];
    A_sym = jacobian(augmented_acceleration_full, [x, y, z, vx, vy, vz, CD]);

    % Convert to a numeric function
    A_func = matlabFunction(A_sym, 'Vars', {x, y, z, vx, vy, vz, CD});
end

function [acceleration, A] = sc_Dynamics(position_eci, velocity_eci, coeffs_vals, mu_val, ...
    ~, CD, t, acc_SPH, A_SPH, acc_Drag, A_Drag, ~)

    % Convert ECI to ECEF coordinates
    R_eci_to_ecef = DCM_eciToEcef(t);               % Rotation matrix ECI to ECEF
    position_ecef = R_eci_to_ecef * position_eci;   % Transform position to ECEF
 
    % Evaluate acceleration and parameter derivatives at the given inputs
    acceleration_SPH = acc_SPH(position_ecef(1), position_ecef(2), ...
                                     position_ecef(3), mu_val, coeffs_vals');
    acceleration_SPH(1:3) = R_eci_to_ecef' * acceleration_SPH(1:3);
    acceleration_Drag = acc_Drag(position_eci(1), position_eci(2), ...
                                     position_eci(3), velocity_eci(1), velocity_eci(2), ...
                                     velocity_eci(3), CD);

    % Sum acceleration contributions
    acceleration = acceleration_SPH(1:3) + acceleration_Drag(1:3);
    
   if nargout > 1
    % Evaluate Jacobian of zonal harmonics (SPH) at the given position and parameters
    A_SPH_ecef = A_SPH(position_ecef(1), position_ecef(2), position_ecef(3), ...
                    0, 0, 0, mu_val, coeffs_vals');

    % Evaluate Jacobian of drag in ECI
    A_Drag_eci = A_Drag(position_eci(1), position_eci(2), ...
                                 position_eci(3), velocity_eci(1), velocity_eci(2), ...
                                 velocity_eci(3), CD);

    % Transform zonal harmonics Jacobian back to ECI
    Ar_eci = R_eci_to_ecef' * A_SPH_ecef(4:6, 1:3) * R_eci_to_ecef; 

    % Extract drag contributions in ECI
    Ar_drag = A_Drag_eci(4:6, 1:3);  % Partial derivatives wrt position
    Av_drag = A_Drag_eci(4:6, 4:6);  % Partial derivatives wrt velocity

    % Sum the contributions from both models
    Ar_total = Ar_eci + Ar_drag;
    Av_total = Av_drag; % Only drag affects velocity derivatives

    % Construct the full state transition matrix A 
    A = [zeros(3, 3), eye(3);      % Velocity derivatives
         Ar_total, Av_total];       % Acceleration derivatives
   end
end

function dstate = sc_ODE(t, state, coeffs, mu, l_max, CD, acc_SPH, A_SPH, ...
    acc_Drag, A_Drag, params_meas)
    % Extract from state
    position_eci = state(1:3);
    velocity_eci = state(4:6);

    % Compute acceleration and A matrix
    [acceleration_eci, A_eci] = sc_Dynamics(position_eci, velocity_eci, coeffs, mu, ...
         l_max, CD, t, acc_SPH, A_SPH, acc_Drag, A_Drag, params_meas);

    % Initialize derivative of state vector
    n = (-1 + sqrt(1 + 4 * size(state, 1))) / 2;
    dstate = zeros(n + n^2, 1);
    dstate(1:3) =  velocity_eci;          
    dstate(4:n) = acceleration_eci;       

    % Propagate STM if requested
    STM = reshape(state(n + 1:end), n, n); % Reshape flattened STM
    Phi_dot = A_eci * STM;                 % STM dynamics
    dstate(n + 1:end) = Phi_dot(:);        % Flatten and append 
end

function plot_earth()
    % Earth's radius in kilometers
    earth_radius = 6378;  % km

    % Create a sphere to represent the Earth
    [x, y, z] = sphere(50);          % Create a sphere with a 50x50 grid
    
    % Scale the sphere to Earth's radius
    x = x * earth_radius;
    y = y * earth_radius;
    z = z * earth_radius;

    % Load the topographic data provided by MATLAB (you can use your own Earth texture if needed)
    load topo;                       % Load MATLAB's topographic data

    % Plot the spherical surface
    s = surface(x, y, z);
    
    % Apply texture mapping and set topographic data as the color data
    s.FaceColor = 'texturemap';      % Use texture mapping
    s.CData = topo;                  % Set color data to topographic data
    s.EdgeColor = 'none';            % Remove edges
    s.FaceLighting = 'gouraud';      % Set Gouraud lighting for smooth curved surface lighting
    s.SpecularStrength = 0.4;        % Adjust strength of the specular reflection

    % Add a light source
    light('Position', [-1, 0, 1]);   % Position the light
end
