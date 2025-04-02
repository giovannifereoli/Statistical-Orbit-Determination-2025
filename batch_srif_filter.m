function results = batch_srif_filter(trajectory_ref, sorted_measurements, P0)
    % Initialization
    T = length(sorted_measurements); % Number of measurements
    n = size(P0, 1);  % State dimension
    max_m = max(arrayfun(@(s) length(s.residual), sorted_measurements));

    % Cholesky decomposition to form square root information matrix and info vector
    sqrt_info_matrix = chol(inv(P0), 'upper');
    info_vector = sqrt_info_matrix * zeros(n, 1);

    % Extract trajectory times for reference
    trajectory_times = [trajectory_ref.time];
    measurement_times = [sorted_measurements.time];
    [~, traj_indices] = ismember(measurement_times, trajectory_times);

    % Progress tracking
    fprintf('SRIF-Batch: Processing measurements...\n');

    % Process each measurement
    for t = 1:T
        traj_idx = traj_indices(t);
        STM_t = squeeze(trajectory_ref(traj_idx).STM);

        residual = sorted_measurements(t).residual;
        H_tilde = sorted_measurements(t).partials.wrt_X;
        R = sorted_measurements(t).covariance;

        % Whitening transformation
        sqrt_R_inv = chol(inv(R), 'upper');
        whitened_H = sqrt_R_inv * H_tilde;
        whitened_r = sqrt_R_inv * residual;

        % Form augmented system and perform QR factorization
        augmented = [sqrt_info_matrix, info_vector; whitened_H * STM_t, whitened_r];
        [~, R_aug] = qr(augmented);

        sqrt_info_matrix = R_aug(1:n, 1:n);
        info_vector = R_aug(1:n, end);

        % Print progress
        fprintf('\b\b\b\b%3d%%', round((t / T) * 100));
    end

    % Solve for state correction using final info matrix/vector
    dx0 = sqrt_info_matrix \ info_vector;
    covariance_0 = inv(sqrt_info_matrix' * sqrt_info_matrix);

    % Propagate deviation and covariance
    state_deviation_hist = zeros(n, T);
    state_corrected_hist = zeros(n, T);
    P_hist = zeros(n, n, T);
    postfit_residuals = NaN(max_m, T);

    for t = 1:T
        traj_idx = traj_indices(t);
        STM_t = squeeze(trajectory_ref(traj_idx).STM);
        x_ref = [trajectory_ref(traj_idx).state(1:6)'; trajectory_ref(traj_idx).parameters(:)];

        dx_t = STM_t * dx0;
        P_t = STM_t * covariance_0 * STM_t';

        H_tilde = sorted_measurements(t).partials.wrt_X;
        residual = sorted_measurements(t).residual;
        postfit_res = residual - H_tilde * dx_t;

        state_deviation_hist(:, t) = dx_t;
        state_corrected_hist(:, t) = x_ref + dx_t;
        P_hist(:, :, t) = P_t;
        postfit_residuals(1:length(residual), t) = postfit_res;

        fprintf('\b\b\b\b%3d%%', round((t / T) * 100));
    end

    fprintf('\nSRIF-Batch Filter Completed!\n');

    % Store results
    results = struct(...
        'state_corrected_hist', state_corrected_hist, ...
        'state_deviation_hist', state_deviation_hist, ...
        'P_hist', P_hist, ...
        'postfit_residuals', postfit_residuals ...
    );
end
