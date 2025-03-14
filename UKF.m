function results = ukf(trajectory_ref, sorted_measurements, P0, settings_UKF, settings_PN)
    % Extract state and parameter dimensions dynamically
    n = size(trajectory_ref(1).state, 2);  % Number of state variables 
    n_par = size(trajectory_ref(1).parameters, 2);  % Number of parameters in f(x)
    n_full = n + n_par;  % Total number of variables in f(x)

    % UKF scaling parameters from settings
    alpha = settings_UKF.alpha;
    beta = settings_UKF.beta;
    kappa = 3 - n_full;  
    lambda = alpha^2 * (n_full + kappa) - n_full;
    gamma = sqrt(n_full + lambda);

    % Precompute weights for the UKF
    num_sigma = 2 * n_full + 1;
    c = 0.5 / (n_full + lambda);
    W_m = repmat(c, 1, num_sigma);
    W_c = repmat(c, 1, num_sigma);
    W_m(1) = lambda / (n_full + lambda);
    W_c(1) = W_m(1) + (1 - alpha^2 + beta);

    % Initialization
    T = length(sorted_measurements);
    max_m = max(arrayfun(@(s) length(s.residual), sorted_measurements));
    
    state_corrected_hist = zeros(n_full, T);  
    state_deviation_hist = zeros(n_full, T);  
    P_hist = zeros(n_full, n_full, T);
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
        % Extract measurement time
        t_now = measurement_times(t);
        if t > 1
            t_prev = measurement_times(t - 1);
        else
            t_prev = trajectory_times(1);
        end
        
        traj_idx = traj_indices(t);
        propagate_func = trajectory_ref(traj_idx).function_handle;

        % Generate Sigma Points
        X_sigma = generate_sigma_points(x_aug, P, gamma);

        % Prediction Step - Batch Process Sigma Points
        t_span = [t_prev, t_now];
        if t_prev ~= t_now
            % Propagate if not t0
            [~, X_propagated] = ode113(@(t, X) batch_propagate_func(t, X, ...
                propagate_func, n_full), t_span, X_sigma(:), ode_options);
            
            % Reshape results back into sigma point matrix
            X_pred = reshape(X_propagated(end, :), n_full, num_sigma);
        else
            X_pred = X_sigma;
        end

        % Compute process noise (if available)
        if nargin > 4
            dt = 0;
            if t > 1
                dt = measurement_times(t) - measurement_times(t-1);
            end
            Q_disc = zeros(n_full);
            Q_disc(1:6, 1:6) = calculate_Q_discrete(dt, x_aug(1:6), ...
                settings_PN);
        else
            Q_disc = zeros(n_full);
        end

        % Compute weighted mean and covariance of transformed sigma points
        [x_pred, P_pred] = unscented_transform(X_pred, Q_disc, W_m, W_c);

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
        P_xz = zeros(n_full, length(R));
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
function dX = batch_propagate_func(t, X, propagate_func, n_full)
    % Compute the number of sigma points
    num_sigma = size(X, 1) / n_full;

    % Reshape into matrix form where each column is a sigma point of size n_full
    X_reshaped = reshape(X, n_full, num_sigma);
    dX_reshaped = zeros(size(X_reshaped));

    % Compute propagation for each sigma point in batch
    for i = 1:num_sigma
        dX_reshaped(:, i) = propagate_func(t, X_reshaped(:, i), 0);
    end

    % Reshape back into a vector for the ODE solver
    dX = dX_reshaped(:);
end

% Sigma Point Generation Function
function X_sigma = generate_sigma_points(x, P, gamma)
    sqrtP = sqrtm(P);
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

