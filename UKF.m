function results = ukf(trajectory_ref, sorted_measurements, P0, settings_UKF, settings_PN)
    % Extract state and parameter dimensions dynamically
    n = size(trajectory_ref(1).state, 2);  % Number of state variables 
    n_par = size(trajectory_ref(1).parameters, 2);  % Number of parameters in f(x)
    n_full = n + n_par;  % Total number of variables in f(x)

    % UKF scaling parameters from settings
    alpha = settings_UKF.alpha;
    beta = settings_UKF.beta;
    kappa = 3 - n;  
    lambda = alpha^2 * (n + kappa) - n;
    gamma = sqrt(n + lambda);

    % Precompute weights for the UKF
    num_sigma = 2 * n + 1;
    c = 0.5 / (n + lambda);
    W_m = repmat(c, 1, num_sigma);
    W_c = repmat(c, 1, num_sigma);
    W_m(1) = lambda / (n + lambda);
    W_c(1) = W_m(1) + (1 - alpha^2 + beta);

    % Initialization
    T = length(sorted_measurements);
    max_m = max(arrayfun(@(s) length(s.residual), sorted_measurements));
    
    state_corrected_hist = zeros(n, T);  
    state_deviation_hist = zeros(n, T);  
    P_hist = zeros(n, n, T);
    postfit_residuals = NaN(max_m, T);  

    % Extract trajectory times for reference
    trajectory_times = [trajectory_ref.time];
    measurement_times = [sorted_measurements.time];
    [~, traj_indices] = ismember(measurement_times, trajectory_times);

    % Initial state estimate
    x0 = trajectory_ref(traj_indices(1)).state(:);
    par0 = trajectory_ref(traj_indices(1)).parameters(:);
    x_aug = [x0; par0];  
    P = P0;  

    % ODE options
    ode_options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    fprintf('UKF Progress: 0%%');

    % UKF Loop
    for t = 1:T
        t_now = measurement_times(t);
        t_prev = (t > 1) * measurement_times(t - 1) + (t == 1) * trajectory_times(1);
        
        traj_idx = traj_indices(t);
        propagate_func = trajectory_ref(traj_idx).function_handle;

        % Generate Sigma Points
        X_sigma = generate_sigma_points(x_aug, P, gamma);

        % Prediction Step - Batch Process Sigma Points
        fprintf('\nPropagating Sigma Points for Time Step %d/%d: ', t, T);
        
        % Stack sigma points into a matrix to propagate in batch
        t_span = [t_prev, t_now];
        [~, X_propagated] = ode113(@(t, X) batch_propagate_func(t, X, ...
            propagate_func), t_span, X_sigma(:), ode_options);
        
        % Reshape results back into sigma point matrix
        X_pred = reshape(X_propagated(end, :), n_full, num_sigma);
        
        % Compute weighted mean and covariance of transformed sigma points
        [x_pred, P_pred] = unscented_transform(X_pred, P, W_m, W_c);

        % Compute process noise (if available)
        if nargin > 4
            dt = 0;
            if t > 1
                dt = measurement_times(t) - measurement_times(t-1);
            end
            Q_disc = calculate_Q_discrete(dt, x_pred(1:6), settings_PN);
        else
            Q_disc = zeros(n);
        end
        P_pred = P_pred + Q_disc;

        % Recompute Sigma Points for Measurement Update
        X_sigma = generate_sigma_points(x_pred, P_pred, gamma);

        % Measurement Update
        measurement_function = sorted_measurements(t).function_handle;
        R = diag(sorted_measurements(t).covariance);

        % Transform sigma points through measurement function
        Z_sigma = zeros(length(R), num_sigma);
        for i = 1:num_sigma
            Z_sigma(:, i) = measurement_function(X_sigma(:, i));
        end
        
        % Compute predicted measurement mean and covariance
        [z_pred, P_zz] = unscented_transform(Z_sigma, R, W_m, W_c);

        % Compute cross covariance
        P_xz = zeros(n, length(R));
        for i = 1:num_sigma
            P_xz = P_xz + W_c(i) * (X_sigma(:, i) - x_pred) * (Z_sigma(:, i) - z_pred)';
        end
        
        % Kalman Gain
        K = P_xz / P_zz;
        
        % Compute residual
        observed_measurement = sorted_measurements(t).observed;
        postfit_res = observed_measurement - z_pred;
        
        % Update step
        x_upd = x_pred + K * postfit_res;
        P_upd = P_pred - K * P_zz * K';

        % Store results
        state_corrected_hist(:, t) = x_upd;
        state_deviation_hist(:, t) = x_upd - x_pred;
        P_hist(:, :, t) = P_upd;
        postfit_residuals(1:length(observed_measurement), t) = postfit_res;

        % Update for next iteration
        x_aug = x_upd;
        P = P_upd;

        fprintf('\b\b\b\b%3d%%', round((t / T) * 100));
    end

    fprintf('\nUKF Completed!\n');

    % Store all outputs
    results = struct(...
        'state_corrected_hist', state_corrected_hist, ...
        'state_deviation_hist', state_deviation_hist, ... 
        'P_hist', P_hist, ...
        'postfit_residuals', postfit_residuals ...
    );
end

% Batch Propagation Function
function dX = batch_propagate_func(t, X, propagate_func)
    num_sigma = size(X, 1);
    X_reshaped = reshape(X, [], num_sigma);
    dX = zeros(size(X_reshaped));
    
    % Compute propagation for each sigma point in parallel
    for i = 1:num_sigma
        dX(:, i) = propagate_func(t, X_reshaped(:, i));
    end
    
    dX = dX(:);
end

% Sigma Point Generation Function
function X_sigma = generate_sigma_points(x, P, gamma)
    sqrtP = chol(P, 'upper');  % Upper triangular matrix
    X_sigma = [x, x + gamma * sqrtP, x - gamma * sqrtP];
end

% Unscented Transform Function
function [x_mean, P_cov] = unscented_transform(X_sigma, noise_cov, W_m, W_c)
    % Initialize
    num_sigma = size(X_sigma, 2);

    % Compute mean
    x_mean = sum(W_m .* X_sigma, 2);

    % Compute covariance
    P_cov = noise_cov;
    for i = 1:num_sigma
        P_cov = P_cov + W_c(i) * (X_sigma(:, i) - x_mean) * (X_sigma(:, i) - x_mean)';
    end
end
