function results = SRIF(trajectory_ref, sorted_measurements, P0, settings_PN)
    % Initialization
    T = length(sorted_measurements); % Number of measurements
    n = size(P0, 1);  % State dimension
    max_m = max(arrayfun(@(s) length(s.residual), sorted_measurements));

    state_deviation_hist = zeros(n, T);
    state_corrected_hist = zeros(n, T);
    P_hist = zeros(n, n, T);
    postfit_residuals = NaN(max_m, T);  % Preallocate with NaN

    % Initial state deviation
    dx = zeros(n, 1);
    % Compute initial square root information matrix (upper triangular R)
    R_sqrt = chol(inv(P0), 'upper');
    info_vector = R_sqrt * dx;

    % Identity matrix for stability
    I = eye(n);
    STM_tm = I;  % Initial transition matrix (identity)

    % Extract trajectory times for reference
    trajectory_times = [trajectory_ref.time];
    measurement_times = [sorted_measurements.time];
    [~, traj_indices] = ismember(measurement_times, trajectory_times);

    % Progress tracking
    fprintf('SRIF Progress: 0%%');

    % SRIF loop over measurements
    for t = 1:T
        % Extract corresponding STM and state from trajectory
        traj_idx = traj_indices(t);
        STM_t = squeeze(trajectory_ref(traj_idx).STM);
        x_ref = trajectory_ref(traj_idx).state(1:n)';

        % Compute STM transition matrix
        STM_dt = STM_t / STM_tm;

        % Extract measurement data
        prefit_residual = sorted_measurements(t).residual;
        H_tilde = sorted_measurements(t).partials.wrt_X;
        R = diag(sorted_measurements(t).covariance);
        
        % Compute Process Noise
        if nargin > 3
            dt = 0;
            if t > 1
                dt = measurement_times(t) - measurement_times(t-1);
            end
            Q_disc = calculate_Q_discrete(dt, x_ref, settings_PN);
        else 
            Q_disc = zeros(n);
        end
        
        % Prediction Step: Propagate information matrix 
        % OSS: info_vector prediction is useless!
        R_sqrt = R_sqrt / STM_dt;

        % Apply process noise correction (QR decomposition)
        if any(Q_disc(:) ~= 0)
            [Q, R_sqrt] = qr([chol(Q_disc, 'upper'); R_sqrt]);
            info_vector = Q' * [zeros(n, 1); info_vector];
            info_vector = info_vector(n+1:end);
        end

        % Measurement Update: Use QR decomposition to incorporate new information
        % HP: SRIF assumes whitened observations.
        meas_SRI = chol(inv(R), 'upper');
        H_tilde_whitened = meas_SRI * H_tilde;
        residuals_whitened = meas_SRI * prefit_residual;
        
        % Construct and solve augmented system
        % OSS: QR decomposition acts as Householder transormation.
        update_matrix = [R_sqrt, info_vector; H_tilde_whitened, residuals_whitened];
        [~, R_new] = qr(update_matrix);
        R_sqrt = R_new(1:n, 1:n);
        info_vector = R_new(1:n, n+1);
        
        % Compute updated state deviation
        dx_upd = R_sqrt \ info_vector;

        % Store results
        state_deviation_hist(:, t) = dx_upd;
        state_corrected_hist(:, t) = x_ref + dx_upd;
        P_hist(:, :, t) = inv(R_sqrt' * R_sqrt);
        postfit_residuals(1:length(prefit_residual), t) = prefit_residual - H_tilde * dx_upd;

        % Update for next iteration
        STM_tm = STM_t;
        
        % Print progress
        fprintf('\b\b\b\b%3d%%', round((t / T) * 100));
    end

    fprintf('\nSRIF Completed!\n');

    % Store all outputs in a struct
    results = struct( ...
        'state_corrected_hist', state_corrected_hist, ...
        'state_deviation_hist', state_deviation_hist, ...
        'P_hist', P_hist, ...
        'postfit_residuals', postfit_residuals ...
    );
end