function results = pesf(trajectory_ref, sorted_measurements, P0, settings_PN)
    % Pseudo-epoch State Filter (PESF) implementation
    % Inputs:
    %   trajectory_ref: Reference trajectory data (struct array with time, state, STM, Phi_xp, Phi_xc, parameters, Pc0)
    %   sorted_measurements: Measurement data (struct array with time, residual, partials.wrt_X, partials.wrt_C, covariance)
    %   P0: Initial state covariance matrix
    %   settings_PN: Process noise settings (optional, struct with Q_cont, tau_p, dt_batch)
    % Outputs:
    %   results: Struct containing corrected state, parameters, stochastics, covariances, and residuals

    % Initialization
    T = length(sorted_measurements); % Number of measurements
    n = size(P0, 1); % State dimension
    Np = 3; % Number of stochastic parameters
    c0 = trajectory_ref(1).parameters; % Initial parameter estimate from trajectory_ref
    Nc = length(c0); % Number of parameters
    max_m = max(arrayfun(@(s) length(s.residual), sorted_measurements));

    % Initialize histories
    state_deviation_hist = zeros(n, T);
    state_corrected_hist = zeros(n, T);
    chat_hist = zeros(Nc, T);
    P_xx_hist = zeros(n, n, T);
    P_xc_hist = zeros(n, Nc, T);
    P_cc_hist = zeros(Nc, Nc, T);
    P_pp_hist = zeros(Np, Np, T);
    P_px_hist = zeros(Np, n, T);
    P_pc_hist = zeros(Np, Nc, T);
    postfit_residuals = NaN(max_m, T);

    % Initial state and parameter deviation
    dx = zeros(n, 1);
    dc = zeros(Nc, 1);
    x0 = trajectory_ref(1).state(1:n)';

    % Compute initial square root information matrices (upper triangular)
    Pc0 = trajectory_ref(1).Pc0; % Initial parameter covariance from trajectory_ref
    R_sqrt = chol(inv(P0), 'upper');
    Rchat0 = chol(inv(Pc0), 'upper');
    Rphat0 = chol(inv(1e-20 * eye(Np)), 'upper');
    Rpxhat0 = zeros(Np, n);
    Rpchat0 = zeros(Np, Nc);
    Rxphat0 = zeros(n, Np);
    Rxchat0 = zeros(n, Nc);
    info_vector = R_sqrt * dx;
    zchat0 = Rchat0 * dc;
    zphat0 = zeros(Np, 1);

    % Process noise settings
    if nargin > 3 && ~isempty(settings_PN)
        Rw_cov = settings_PN.Q_cont;
        Rw = chol(inv(Rw_cov), 'upper');
        tau_p = settings_PN.tau_p;
        dt_batch = settings_PN.dt_batch;
    else
        Rw_cov = (1e-10^2) * eye(Np);
        Rw = chol(inv(Rw_cov), 'upper');
        tau_p = 100 * 3600;
        dt_batch = 500 * 3600;
    end
    M = (-1 / tau_p) * eye(Np);

    % Extract trajectory and measurement times
    trajectory_times = [trajectory_ref.time];
    measurement_times = [sorted_measurements.time];
    [~, traj_indices] = ismember(measurement_times, trajectory_times);

    % Compute batch times
    t0 = 0;
    tend = measurement_times(end);
    ts_batches = t0:dt_batch:tend;
    if ts_batches(end) == tend
        ts_batches(end) = [];
    end
    N_batches = length(ts_batches);
    p_all = zeros(Np, N_batches);

    % Assign measurements to batches
    ind_obs_batches = cell(N_batches, 1);
    for ii = 1:N_batches-1
        if ii == 1
            ind_obs_batches{ii} = find(measurement_times >= ts_batches(ii) & measurement_times <= ts_batches(ii+1));
        else
            ind_obs_batches{ii} = find(measurement_times > ts_batches(ii) & measurement_times <= ts_batches(ii+1));
        end
    end
    if N_batches == 1
        ind_obs_batches{N_batches} = find(measurement_times >= ts_batches(end));
    else
        ind_obs_batches{N_batches} = find(measurement_times > ts_batches(end));
    end

    % Progress tracking
    fprintf('PESF Progress: 0%%');

    % Initialize filter variables
    Rphat = Rphat0;
    Rpxhat = Rpxhat0;
    Rpchat = Rpchat0;
    Rxphat = Rxphat0;
    R_sqrt = R_sqrt;
    Rxchat = Rxchat0;
    Rchat = Rchat0;
    info_vector = info_vector;
    zchat = zchat0;
    zphat = zphat0;
    obs_num = 1;

    % Check for measurement at initial time
    special_first_meas = measurement_times(1) == t0;

    % PESF loop over batches
    for ii = 1:N_batches
        zphat = zphat0;

        % Handle special first measurement
        if special_first_meas
            % Measurement Update
            traj_idx = traj_indices(ind_obs_batches{ii}(1));
            x_ref = trajectory_ref(traj_idx).state(1:n)';
            prefit_residual = sorted_measurements(ind_obs_batches{ii}(1)).residual;
            H_tilde = sorted_measurements(ind_obs_batches{ii}(1)).partials.wrt_X;
            Hc_tilde = sorted_measurements(ind_obs_batches{ii}(1)).partials.wrt_C;
            R_meas = diag(sorted_measurements(ind_obs_batches{ii}(1)).covariance);
            meas_SRI = chol(inv(R_meas), 'upper');
            H_tilde_whitened = meas_SRI * H_tilde;
            Hc_tilde_whitened = meas_SRI * Hc_tilde;
            residuals_whitened = meas_SRI * prefit_residual;
            update_matrix = [Rphat Rpxhat Rpchat zphat; ...
                             Rxphat R_sqrt Rxchat info_vector; ...
                             zeros(length(prefit_residual), Np) H_tilde_whitened Hc_tilde_whitened residuals_whitened];
            [~, R_new] = qr(update_matrix);
            Rphat = R_new(1:Np, 1:Np);
            Rpxhat = R_new(1:Np, Np+1:Np+n);
            Rpchat = R_new(1:Np, Np+n+1:Np+n+Nc);
            zphat = R_new(1:Np, end);
            Rxphat = R_new(Np+1:Np+n, 1:Np);
            R_sqrt = R_new(Np+1:Np+n, Np+1:Np+n);
            Rxchat = R_new(Np+1:Np+n, Np+n+1:Np+n+Nc);
            info_vector = R_new(Np+1:Np+n, end);
            Hchat = R_new(Np+n+1:end, Np+n+1:Np+n+Nc);
            pf_resid_hat = R_new(Np+n+1:end, end);
            Sc = [Rchat zchat; Hchat pf_resid_hat];
            [~, Sout] = qr(Sc);
            Rchat = Sout(1:Nc, 1:Nc);
            zchat = Sout(1:Nc, end);
            postfit_residuals(1:length(prefit_residual), obs_num) = Sout(Nc+1:end, end);
            dx_upd = R_sqrt \ (info_vector - Rxchat * (Rchat \ zchat));
            state_deviation_hist(:, obs_num) = dx_upd;
            state_corrected_hist(:, obs_num) = x_ref + dx_upd;
            chat_hist(:, obs_num) = Rchat \ zchat;
            P_xx_hist(:, :, obs_num) = inv(R_sqrt' * R_sqrt);
            obs_num = obs_num + 1;
        end

        % Process measurements in the batch
        for jj = 1:length(ind_obs_batches{ii})
            if special_first_meas && jj == 1
                special_first_meas = 0;
            else
                % Time Update
                traj_idx = traj_indices(ind_obs_batches{ii}(jj));
                x_ref = trajectory_ref(traj_idx).state(1:n)';
                Phi_xx = squeeze(trajectory_ref(traj_idx).STM);
                Phi_xp = squeeze(trajectory_ref(traj_idx).Phi_xp);
                Vp = Phi_xx \ Phi_xp;
                Shat = [-Rw * M Rw zeros(Np, n) zeros(Np, Nc) zeros(Np, 1); ...
                        (Rphat - Rpxhat * Vp) zeros(Np, Np) Rpxhat Rpchat zphat; ...
                        (Rxphat - R_sqrt * Vp) zeros(n, Np) R_sqrt Rxchat info_vector];
                [~, Sout] = qr(Shat);
                Rphat = Sout(Np+1:Np+Np, Np+1:Np+Np);
                Rpxhat = Sout(Np+1:Np+Np, Np+Np+1:Np+Np+n);
                Rpchat = Sout(Np+1:Np+Np, Np+Np+n+1:Np+Np+n+Nc);
                zphat = Sout(Np+1:Np+Np, end);
                Rxphat = Sout(Np+Np+1:end, Np+1:Np+Np);
                R_sqrt = Sout(Np+Np+1:end, Np+Np+1:Np+Np+n);
                Rxchat = Sout(Np+Np+1:end, Np+Np+n+1:Np+Np+n+Nc);
                info_vector = Sout(Np+Np+1:end, end);
                Rchat = Rchat; % No update to Rchat
                zchat = zchat; % No update to zchat

                % Measurement Update
                prefit_residual = sorted_measurements(ind_obs_batches{ii}(jj)).residual;
                H_tilde = sorted_measurements(ind_obs_batches{ii}(jj)).partials.wrt_X;
                Hc_tilde = sorted_measurements(ind_obs_batches{ii}(jj)).partials.wrt_C;
                R_meas = diag(sorted_measurements(ind_obs_batches{ii}(jj)).covariance);
                meas_SRI = chol(inv(R_meas), 'upper');
                H_tilde_whitened = meas_SRI * H_tilde;
                Hc_tilde_whitened = meas_SRI * Hc_tilde;
                residuals_whitened = meas_SRI * prefit_residual;
                update_matrix = [Rphat Rpxhat Rpchat zphat; ...
                                 Rxphat R_sqrt Rxchat info_vector; ...
                                 zeros(length(prefit_residual), Np) H_tilde_whitened Hc_tilde_whitened residuals_whitened];
                [~, R_new] = qr(update_matrix);
                Rphat = R_new(1:Np, 1:Np);
                Rpxhat = R_new(1:Np, Np+1:Np+n);
                Rpchat = R_new(1:Np, Np+n+1:Np+n+Nc);
                zphat = R_new(1:Np, end);
                Rxphat = R_new(Np+1:Np+n, 1:Np);
                R_sqrt = R_new(Np+1:Np+n, Np+1:Np+n);
                Rxchat = R_new(Np+1:Np+n, Np+n+1:Np+n+Nc);
                info_vector = R_new(Np+1:Np+n, end);
                Hchat = R_new(Np+n+1:end, Np+n+1:Np+n+Nc);
                pf_resid_hat = R_new(Np+n+1:end, end);
                Sc = [Rchat zchat; Hchat pf_resid_hat];
                [~, Sout] = qr(Sc);
                Rchat = Sout(1:Nc, 1:Nc);
                zchat = Sout(1:Nc, end);
                postfit_residuals(1:length(prefit_residual), obs_num) = Sout(Nc+1:end, end);
                dx_upd = R_sqrt \ (info_vector - Rxchat * (Rchat \ zchat));
                state_deviation_hist(:, obs_num) = dx_upd;
                state_corrected_hist(:, obs_num) = x_ref + dx_upd;
                chat_hist(:, obs_num) = Rchat \ zchat;
                P_xx_hist(:, :, obs_num) = inv(R_sqrt' * R_sqrt);
                obs_num = obs_num + 1;

                % Double time update at batch end
                if jj == length(ind_obs_batches{ii}) && ii < N_batches && measurement_times(ind_obs_batches{ii}(end)) < ts_batches(ii+1)
                    Phi_xx_end = squeeze(trajectory_ref(traj_indices(ind_obs_batches{ii}(end))).STM);
                    Phi_xp_end = squeeze(trajectory_ref(traj_indices(ind_obs_batches{ii}(end))).Phi_xp);
                    Vp = Phi_xx_end \ Phi_xp_end;
                    Shat = [-Rw * M Rw zeros(Np, n) zeros(Np, Nc) zeros(Np, 1); ...
                            (Rphat - Rpxhat * Vp) zeros(Np, Np) Rpxhat Rpchat zphat; ...
                            (Rxphat - R_sqrt * Vp) zeros(n, Np) R_sqrt Rxchat info_vector];
                    [~, Sout] = qr(Shat);
                    Rphat = Sout(Np+1:Np+Np, Np+1:Np+Np);
                    Rpxhat = Sout(Np+1:Np+Np, Np+Np+1:Np+Np+n);
                    Rpchat = Sout(Np+1:Np+Np, Np+Np+n+1:Np+Np+n+Nc);
                    zphat = Sout(Np+1:Np+Np, end);
                    Rxphat = Sout(Np+Np+1:end, Np+1:Np+Np);
                    R_sqrt = Sout(Np+Np+1:end, Np+Np+1:Np+Np+n);
                    Rxchat = Sout(Np+Np+1:end, Np+Np+n+1:Np+Np+n+Nc);
                    info_vector = Sout(Np+Np+1:end, end);
                end
            end
        end

        % Solve for estimates via back substitution
        chat = Rchat \ zchat;
        xhat = R_sqrt \ (info_vector - Rxchat * chat);
        phat_end = Rphat \ (zphat - Rpxhat * xhat - Rpchat * chat);
        if ii == N_batches
            tlast = tend;
        else
            tlast = ts_batches(ii+1);
        end
        Phi_pp_ts = exp(-(tlast - ts_batches(ii)) / tau_p) * eye(Np);
        p_all(:, ii) = Phi_pp_ts \ phat_end;

        % Store covariances
        Rhat_end = [Rphat Rpxhat Rpchat; zeros(n, Np) R_sqrt Rxchat; zeros(Nc, Np) zeros(Nc, n) Rchat];
        P_batch = inv(Rhat_end' * Rhat_end);
        P_pp_hist(:, :, ii) = P_batch(1:Np, 1:Np);
        P_px_hist(:, :, ii) = P_batch(1:Np, Np+1:Np+n);
        P_pc_hist(:, :, ii) = P_batch(1:Np, Np+n+1:Np+n+Nc);
        P_xx_hist(:, :, ii) = P_batch(Np+1:Np+n, Np+1:Np+n);
        P_xc_hist(:, :, ii) = P_batch(Np+1:Np+n, Np+n+1:Np+n+Nc);
        P_cc_hist(:, :, ii) = P_batch(Np+n+1:Np+n+Nc, Np+n+1:Np+n+Nc);

        % Absorb stochastic effect
        info_vector = info_vector - Rxphat * p_all(:, ii);
        zchat = zchat - Rpchat * p_all(:, ii);

        % Reset stochastic information
        Rphat = Rphat0;
        Rpxhat = Rpxhat0;
        Rpchat = Rpchat0;
        Rxphat = Rxphat0;
        zphat = zphat0;

        % Print progress
        fprintf('\b\b\b\b%3d%%', round((obs_num / T) * 100));
    end

    fprintf('\nPESF Completed!\n');

    % Store results in a struct
    results = struct( ...
        'state_corrected_hist', state_corrected_hist, ...
        'state_deviation_hist', state_deviation_hist, ...
        'chat_hist', chat_hist, ...
        'P_xx_hist', P_xx_hist, ...
        'P_xc_hist', P_xc_hist, ...
        'P_cc_hist', P_cc_hist, ...
        'P_pp_hist', P_pp_hist, ...
        'P_px_hist', P_px_hist, ...
        'P_pc_hist', P_pc_hist, ...
        'p_all', p_all, ...
        'postfit_residuals', postfit_residuals ...
    );
end
