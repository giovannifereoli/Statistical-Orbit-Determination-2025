function results = ekf(trajectory_ref, sorted_measurements, P0, settings_PN)
    % Extract state and parameter dimensions dynamically
    n = size(trajectory_ref(1).state, 2);                   % Number of state variables 
    n_par = size(trajectory_ref(1).parameters, 2);          % Number of parameters in f(x)
    n_full = n + n_par;                                     % Number of total variables in f(x)

    % TODO: initialize with LKF!

    % Initialization
    T = length(sorted_measurements); % Number of measurements
    max_m = max(arrayfun(@(s) length(s.residual), sorted_measurements));

    state_corrected_hist = zeros(n_full, T);  % Corrected state estimates
    state_deviation_hist = zeros(n_full, T);  % Cumulative state deviations
    P_hist = zeros(n_full, n_full, T);
    postfit_residuals = NaN(max_m, T);  % Preallocate with NaN

    % Extract trajectory times for reference
    trajectory_times = [trajectory_ref.time];

    % Find indices of trajectory points that match measurement times
    measurement_times = [sorted_measurements.time];
    [~, traj_indices] = ismember(measurement_times, trajectory_times);

    % Initial state estimate (from first trajectory reference)
    x0 = trajectory_ref(traj_indices(1)).state(:);
    par0 = trajectory_ref(traj_indices(1)).parameters(:);
    STM0 = eye(n_full);  % Initial STM
    x_aug = [x0; par0; STM0(:)];  % Augmented state with STM
    P = P0;  % Initial covariance

    % ODE options for propagation
    ode_options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % Control LKF/EKF switching logic
    LKF_steps = 100;
    gap_threshold = 3 * 3600; % seconds (e.g., 3 hours)

    % Progress tracking
    fprintf('EKF Progress: 0%%');

    % EKF loop over measurements
    for t = 1:T
        % Extract measurement time
        t_now = measurement_times(t);
        if t > 1
            t_prev = measurement_times(t - 1);
        else
            t_prev = trajectory_times(1);
        end

        % Extract stored function handle for state propagation
        traj_idx = traj_indices(t);
        propagate_func = trajectory_ref(traj_idx).function_handle;

        % Prediction Step
        t_span = [t_prev, t_now];  % Propagate to next measurement time
        if t_prev ~= t_now
            % Propagate Dynamics
            [~, x_propagated] = ode113(@(t, x) propagate_func(t, x), t_span,...
                x_aug, ode_options);
    
            % Extract STM (skip parameters)
            STM_t = reshape(x_propagated(end, n_full + 1:end)', n_full, n_full); 

            % Compute Process Noise
            if nargin > 3
                dt = 0;  % Default value for the first iteration
                if t > 1
                    dt = measurement_times(t) - measurement_times(t-1); % Time step
                end
                Q_disc = zeros(n_full);
                Q_disc(1:6, 1:6) = calculate_Q_discrete(dt, x_propagated(end, 1:6), settings_PN);  % Store Q_disc
            else 
                Q_disc = zeros(n_full);
            end
    
            % Do actual Prediction Step
            x_pred = x_propagated(end, 1:n_full)';  
            P_pred = STM_t * P * STM_t' + Q_disc;
        else
            % No propagation needed, use current state
            x_pred = x_aug(1:n_full);  
            P_pred = P;
        end

        % Extract Measurement Model Function Handle and Covariance
        measurement_function = sorted_measurements(t).measurement_function;
        R = sorted_measurements(t).covariance;

        % Compute Computed Measurement and Partial
        [h_pred, H_tilde] = measurement_function(x_pred);   

        % Compute postfit-fit residual
        observed_measurement = sorted_measurements(t).observed;
        postfit_res = observed_measurement - h_pred;

        % Kalman Gain computation
        S = H_tilde * P_pred * H_tilde' + R;
        K = P_pred * H_tilde' / S;

        % Update Step
        dx_upd = K * postfit_res;
        x_upd = x_pred + dx_upd;
        P_upd = (eye(n_full) - K * H_tilde) * P_pred * (eye(n_full) - K * H_tilde)' + K * R * K';

        % Store results
        state_corrected_hist(:, t) = x_upd;
        state_deviation_hist(:, t) = dx_upd; 
        P_hist(:, :, t) = P_upd;
        postfit_residuals(1:length(observed_measurement), t) = postfit_res;

        % Update for next iteration
        use_lkf = (t <= LKF_steps) || (t > 1 && abs(t_now - t_prev) > gap_threshold);
        if use_lkf
            x_aug = [x_pred; STM0(:)]; % Use fixed reference for LKF behavior
        else
            x_aug = [x_upd; STM0(:)]; % Use updated reference
        end
        P = P_upd;

        % Print progress
        fprintf('\b\b\b\b%3d%%', round((t / T) * 100));
    end

    fprintf('\nEKF Completed!\n');

    % Store all outputs in a struct (dictionary)
    results = struct(...
        'state_corrected_hist', state_corrected_hist, ...
        'state_deviation_hist', state_deviation_hist, ... 
        'P_hist', P_hist, ...
        'postfit_residuals', postfit_residuals ...
    );
end

function Q_disc = calculate_Q_discrete(dt, state, settings)
    % Extract settings from the structure
    Q_cont = settings.Q_cont;
    threshold = settings.threshold;
    frame_type = settings.frame_type;
    method = settings.method;

    % Reset dt if it exceeds the threshold
    dt = min(dt, threshold);

    % If the frame_type is not 'ECI', apply the transformation
    if strcmp(frame_type, 'RIC')
        R_eci_to_ric = transform_to_ric_frame(state);
        Q_cont = R_eci_to_ric' * Q_cont * R_eci_to_ric; 
    end

    if strcmp(method, 'SNC')
        % State Noise Compensation (Q_snc)
        Gamma = [dt^3 * eye(3) / 2; dt^2 * eye(3)];
        Q_snc = Gamma * Q_cont * Gamma';
        Q_disc = Q_snc;  % Return Q_snc directly
        
    elseif strcmp(method, 'DMC')
        % Dynamical Model Compensation (Q_dmc)
        % Extract time-constant matrix
        B = settings.B;

        % Initialize matrices for Q_dmc
        Qrr = eye(3);
        Qrv = eye(3);
        Qrw = eye(3);
        Qvv = eye(3);
        Qvw = eye(3);
        Qww = eye(3);

        % Calculate Q_dmc blocks
        for i = 1:3
            % Calculate intermediate terms
            dt2 = dt^2;
            dt3 = dt^3;
            exp_beta_dt = exp(-B(i,i) * dt);
            exp_2beta_dt = exp(-2 * B(i,i) * dt);

            % Qrr block
            Qrr(i,i) = Q_cont(i,i) * ( ...
                (1 / (3 * B(i,i)^2)) * dt3 ...
                - (1 / (B(i,i)^3)) * dt2 ...
                + (1 / (B(i,i)^4)) * dt ...
                - (2 / (B(i,i)^4)) * dt * exp_beta_dt ...
                + (1 / (2 * B(i,i)^5)) * (1 - exp_2beta_dt) ...
            );

            % Qrv block
            Qrv(i,i) = Q_cont(i,i) * ( ...
                (1 / (2 * B(i,i)^2)) * dt2 ...
                - (1 / (B(i,i)^3)) * dt ...
                + (1 / (B(i,i)^3)) * exp_beta_dt * dt ...
                + (1 / (B(i,i)^4)) * (1 - exp_beta_dt) ...
                - (1 / (2 * B(i,i)^4)) * (1 - exp_2beta_dt) ...
            );

            % Qrw block
            Qrw(i,i) = Q_cont(i,i) * ( ...
                (1 / (2 * B(i,i)^3)) * (1 - exp_2beta_dt) ...
                - (1 / (B(i,i)^2)) * exp_beta_dt * dt ...
            );

            % Qvv block
            Qvv(i,i) = Q_cont(i,i) * ( ...
                (1 / (B(i,i)^2)) * dt ...
                - (2 / (B(i,i)^3)) * (1 - exp_beta_dt) ...
                + (1 / (2 * B(i,i)^3)) * (1 - exp_2beta_dt) ...
            );

            % Qvw block
            Qvw(i,i) = Q_cont(i,i) * ( ...
                (1 / (2 * B(i,i)^2)) * (1 + exp_2beta_dt) ...
                - (1 / (B(i,i)^2)) * exp_beta_dt ...
            );

            % Qww block
            Qww(i,i) = Q_cont(i,i) * ( ...
                (1 / (2 * B(i,i))) * (1 - exp_2beta_dt) ...
            );
        end

        % Combine the blocks into the final Q_dmc matrix
        Q_dmc = [
            Qrr, Qrv, Qrw;
            Qrv, Qvv, Qvw;
            Qrw, Qvw, Qww
        ];

        Q_disc = Q_dmc;
        
    else
        error('Invalid method. Choose "SNC" or "DMC".');
    end
end

function R_eci_to_ric = transform_to_ric_frame(state)
    % Extract state
    position = state(1:3);
    velocity = state(4:6);

    % Normalize the position and velocity vectors
    r = position / norm(position);   % Radial direction (normalized position)
    
    % Compute the cross-track  (C) using the cross-product of r and v
    c = cross(r, velocity);  % Cross-track direction
    c = c / norm(c);  % Normalize the cross-track direction

    % Compute the in-track direction (I) by crossing velocity with radial direction
    i = cross(c, r);             % In-track direction (normalized)
    
    % Define the RIC frame axes
    R = r;  % Radial (R)
    I = i;  % In-track (I)
    C = c;  % Cross-track (C)
    
    % Build the rotation matrix from ECI to RIC
    R_eci_to_ric = [R, I, C];  % 3x3 rotation matrix
end