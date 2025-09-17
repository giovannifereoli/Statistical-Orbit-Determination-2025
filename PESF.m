function results = pesf(trajectory_ref, sorted_measurements, Px0, settings_PN)
    % Pseudo-epoch State Filter (PESF) implementation
    % Inputs:
    %   trajectory_ref: Reference trajectory data (struct array with time, state, STM, etc.)
    %   sorted_measurements: Measurement data (struct array with time, residual, partials, covariance)
    %   Px0: Initial state covariance matrix
    %   settings_PN: Process noise settings (optional, struct with Q_cont, tau_p, etc.)
    % Outputs:
    %   results: Struct containing corrected state history, deviations, covariances, and residuals

    % Initialization
    T = length(sorted_measurements); % Number of measurements
    Nx = size(Px0, 1); % State dimension
    Nc = 2; % Number of parameters (e.g., MU, J2)
    Np = 3; % Number of stochastic parameters
    max_m = max(arrayfun(@(s) length(s.residual), sorted_measurements));

    % A priori parameter covariance
    Pc0 = [0.01 0; 0 1e-8]; % For MU and J2
    Pp0 = 1e-20 * eye(Np); % Initial stochastics covariance

    % Initialize histories
    xhat_batch = zeros(Nx, T);
    chat_batch = zeros(Nc, T);
    P_xx_batch = zeros(Nx, Nx, T);
    P_xc_batch = zeros(Nx, Nc, T);
    P_cc_batch = zeros(Nc, Nc, T);
    postfit_residuals = NaN(max_m, T); % Preallocate with NaN
    meas_RMS_hist = [];

    % Initial state deviation and parameter deviation
    dx = zeros(Nx, 1);
    dc = zeros(Nc, 1);
    dxchat_last = ones(Nx + Nc, 1);
    dxchat_current = zeros(Nx + Nc, 1);

    % Compute initial square root information matrices (upper triangular R)
    Rxhat0 = chol(inv(Px0), 'upper');
    Rchat0 = chol(inv(Pc0), 'upper');
    Rphat0 = chol(inv(Pp0), 'upper');
    Rpxhat0 = zeros(Np, Nx);
    Rpchat0 = zeros(Np, Nc);
    Rxphat0 = zeros(Nx, Np);
    zxhat0 = zeros(Nx, 1);
    zchat0 = zeros(Nc, 1);
    zphat0 = zeros(Np, 1);

    % Filter parameters
    max_iter = 10;
    state_tol = 1e-10;
    meas_tol = 1e-10 * 1.1 * sqrt(3) / 3;
    meas_RMS = 1e3;

    % Process noise settings
    if nargin > 3 && ~isempty(settings_PN)
        Rw_cov = settings_PN.Q_cont;
        Rw = chol(inv(Rw_cov), 'upper');
        tau_p = settings_PN.tau_p;
        M = (-1 / tau_p) * eye(Np);
        dt_batch = settings_PN.dt_batch;
    else
        Rw_cov = (1e-10^2) * eye(Np);
        Rw = chol(inv(Rw_cov), 'upper');
        tau_p = 100 * 3600;
        M = (-1 / tau_p) * eye(Np);
        dt_batch = 500 * 3600;
    end

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

    % PESF iteration loop
    iteration_num = 0;
    x0 = trajectory_ref(1).state(1:Nx)';
    c0 = [3.986004415e5; 1.082626925638815e-3]; % Initial MU, J2

    while (iteration_num < max_iter) && (norm(dxchat_current - dxchat_last) > state_tol) && (meas_RMS > meas_tol)
        iteration_num = iteration_num + 1;

        % Initialize filter variables
        Rphat = Rphat0;
        Rpxhat = Rpxhat0;
        Rpchat = Rpchat0;
        Rxphat = Rxphat0;
        Rxhat = Rxhat0;
        Rxchat = Rxchat0;
        Rchat = Rchat0;
        zxhat = zxhat0;
        zchat = zchat0;

        % Check for measurement at initial time
        special_first_meas = measurement_times(1) == t0;
        obs_num = 1;

        % Process data by batch
        for ii = 1:N_batches
            zphat = zphat0(:, ii);

            % Handle special first measurement
            if special_first_meas
                Rptilde = Rphat;
                Rpxtilde = Rpxhat;
                Rpctilde = Rpchat;
                zptilde = zphat;
                Rxptilde = Rxphat;
                Rxtilde = Rxhat;
                Rxctilde = Rxchat;
                zxtilde = zxhat;
                Rctilde = Rchat;
                zctilde = zchat;

                % Measurement Update
                traj_idx = traj_indices(ind_obs_batches{ii}(1));
                prefit_residual = sorted_measurements(ind_obs_batches{ii}(1)).residual;
                Hx_tilde = sorted_measurements(ind_obs_batches{ii}(1)).partials.wrt_X;
                Hc_tilde = sorted_measurements(ind_obs_batches{ii}(1)).partials.wrt_C;
                Ry = diag(sorted_measurements(ind_obs_batches{ii}(1)).covariance);
                Ry_whiten = chol(Ry, 'lower');
                Hx = Ry_whiten \ Hx_tilde;
                Hc = Ry_whiten \ Hc_tilde;
                whitened_prefit_resid = Ry_whiten \ prefit_residual;
                Stilde = [Rptilde Rpxtilde Rpctilde zptilde; ...
                          Rxptilde Rxtilde Rxctilde zxtilde; ...
                          zeros(length(prefit_residual), Np) Hx Hc whitened_prefit_resid];
                [~, Sout] = qr(Stilde);
                Rphat = Sout(1:Np, 1:Np);
                Rpxhat = Sout(1:Np, Np+1:Np+Nx);
                Rpchat = Sout(1:Np, Np+Nx+1:Np+Nx+Nc);
                zphat = Sout(1:Np, end);
                Rxhat = Sout(Np+1:Np+Nx, Np+1:Np+Nx);
                Rxchat = Sout(Np+1:Np+Nx, Np+Nx+1:Np+Nx+Nc);
                zxhat = Sout(Np+1:Np+Nx, end);
                Hchat = Sout(Np+Nx+1:end, Np+Nx+1:Np+Nx+Nc);
                pf_resid_hat = Sout(Np+Nx+1:end, end);
                Sc = [Rctilde zctilde; Hchat pf_resid_hat];
                [~, Sout] = qr(Sc);
                Rchat = Sout(1:Nc, 1:Nc);
                zchat = Sout(1:Nc, end);
                postfit_residuals(1:length(whitened_prefit_resid), obs_num) = Sout(Nc+1:end, end);
                obs_num = obs_num + 1;
            end

            % Process measurements in the batch
            for jj = 1:length(ind_obs_batches{ii})
                if special_first_meas && jj == 1
                    special_first_meas = 0;
                else
                    % Time Update
                    traj_idx = traj_indices(ind_obs_batches{ii}(jj));
                    Phi_xx = squeeze(trajectory_ref(traj_idx).STM);
                    Phi_xp = squeeze(trajectory_ref(traj_idx).Phi_xp);
                    Vp = Phi_xx \ Phi_xp;
                    Shat = [-Rw * M Rw zeros(Np, Nx) zeros(Np, Nc) zeros(Np, 1); ...
                            (Rphat - Rpxhat * Vp) zeros(Np, Np) Rpxhat Rpchat zphat; ...
                            (Rxphat - Rxhat * Vp) zeros(Nx, Np) Rxhat Rxchat zxhat];
                    [~, Sout] = qr(Shat);
                    Rptilde = Sout(Np+1:Np+Np, Np+1:Np+Np);
                    Rpxtilde = Sout(Np+1:Np+Np, Np+Np+1:Np+Np+Nx);
                    Rpctilde = Sout(Np+1:Np+Np, Np+Np+Nx+1:Np+Np+Nx+Nc);
                    zptilde = Sout(Np+1:Np+Np, end);
                    Rxptilde = Sout(Np+Np+1:end, Np+1:Np+Np);
                    Rxtilde = Sout(Np+Np+1:end, Np+Np+1:Np+Np+Nx);
                    Rxctilde = Sout(Np+Np+1:end, Np+Np+Nx+1:Np+Np+Nx+Nc);
                    zxtilde = Sout(Np+Np+1:end, end);
                    Rctilde = Rchat;
                    zctilde = zchat;

                    % Measurement Update
                    prefit_residual = sorted_measurements(ind_obs_batches{ii}(jj)).residual;
                    Hx_tilde = sorted_measurements(ind_obs_batches{ii}(jj)).partials.wrt_X;
                    Hc_tilde = sorted_measurements(ind_obs_batches{ii}(jj)).partials.wrt_C;
                    Ry = diag(sorted_measurements(ind_obs_batches{ii}(jj)).covariance);
                    Ry_whiten = chol(Ry, 'lower');
                    Hx = Ry_whiten \ Hx_tilde;
                    Hc = Ry_whiten \ Hc_tilde;
                    whitened_prefit_resid = Ry_whiten \ prefit_residual;
                    Stilde = [Rptilde Rpxtilde Rpctilde zptilde; ...
                              Rxptilde Rxtilde Rxctilde zxtilde; ...
                              zeros(length(prefit_residual), Np) Hx Hc whitened_prefit_resid];
                    [~, Sout] = qr(Stilde);
                    Rphat = Sout(1:Np, 1:Np);
                    Rpxhat = Sout(1:Np, Np+1:Np+Nx);
                    Rpchat = Sout(1:Np, Np+Nx+1:Np+Nx+Nc);
                    zphat = Sout(1:Np, end);
                    Rxhat = Sout(Np+1:Np+Nx, Np+1:Np+Nx);
                    Rxchat = Sout(Np+1:Np+Nx, Np+Nx+1:Np+Nx+Nc);
                    zxhat = Sout(Np+1:Np+Nx, end);
                    Hchat = Sout(Np+Nx+1:end, Np+Nx+1:Np+Nx+Nc);
                    pf_resid_hat = Sout(Np+Nx+1:end, end);
                    Sc = [Rctilde zctilde; Hchat pf_resid_hat];
                    [~, Sout] = qr(Sc);
                    Rchat = Sout(1:Nc, 1:Nc);
                    zchat = Sout(1:Nc, end);
                    postfit_residuals(1:length(whitened_prefit_resid), obs_num) = Sout(Nc+1:end, end);
                    obs_num = obs_num + 1;

                    % Double time update at batch end
                    if jj == length(ind_obs_batches{ii}) && ii < N_batches && measurement_times(ind_obs_batches{ii}(end)) < ts_batches(ii+1)
                        Phi_xx_end = squeeze(trajectory_ref(traj_indices(ind_obs_batches{ii}(end))).STM);
                        Phi_xp_end = squeeze(trajectory_ref(traj_indices(ind_obs_batches{ii}(end))).Phi_xp);
                        Vp = Phi_xx_end \ Phi_xp_end;
                        Shat = [-Rw * M Rw zeros(Np, Nx) zeros(Np, Nc) zeros(Np, 1); ...
                                (Rphat - Rpxhat * Vp) zeros(Np, Np) Rpxhat Rpchat zphat; ...
                                (Rxphat - Rxhat * Vp) zeros(Nx, Np) Rxhat Rxchat zxhat];
                        [~, Sout] = qr(Shat);
                        Rphat = Sout(Np+1:Np+Np, Np+1:Np+Np);
                        Rpchat = Sout(Np+1:Np+Np, Np+Np+Nx+1:Np+Np+Nx+Nc);
                        zphat = Sout(Np+1:Np+Np, end);
                        Rxphat = Sout(Np+Np+1:end, Np+1:Np+Np);
                        Rxhat = Sout(Np+Np+1:end, Np+Np+1:Np+Np+Nx);
                        Rxchat = Sout(Np+Np+1:end, Np+Np+Nx+1:Np+Np+Nx+Nc);
                        zxhat = Sout(Np+Np+1:end, end);
                    end
                end
            end

            % Solve for estimates via back substitution
            chat = Rchat \ zchat;
            xhat = Rxhat \ (zxhat - Rxchat * chat);
            phat = Rphat \ (zphat - Rpxhat * xhat - Rpchat * chat);

            % Store batch results
            xhat_batch(:, ii) = xhat;
            chat_batch(:, ii) = chat;
            p_all(:, ii) = p_all(:, ii) + phat;

            % Absorb stochastic effect into state and parameters
            zxhat = zxhat - Rxphat * phat;
            zchat = zchat - Rpchat * phat;

            % Reset stochastic information for next batch
            Rphat = Rphat0;
            Rpxhat = Rpxhat0;
            Rpchat = Rpchat0;
            Rxphat = Rxphat0;

            % Compute covariances
            R_full = [Rphat Rpxhat Rpchat; zeros(Nx, Np) Rxhat Rxchat; zeros(Nc, Np) zeros(Nc, Nx) Rchat];
            P_full = inv(R_full' * R_full);
            P_pp_batch(:, :, ii) = P_full(1:Np, 1:Np);
            P_px_batch(:, :, ii) = P_full(1:Np, Np+1:Np+Nx);
            P_pc_batch(:, :, ii) = P_full(1:Np, Np+Nx+1:end);
            P_xx_batch(:, :, ii) = P_full(Np+1:Np+Nx, Np+1:Np+Nx);
            P_xc_batch(:, :, ii) = P_full(Np+1:Np+Nx, Np+Nx+1:end);
            P_cc_batch(:, :, ii) = P_full(Np+Nx+1:end, Np+Nx+1:end);
        end

        % Update epoch state and parameters
        chat = Rchat \ zchat;
        xhat = Rxhat \ (zxhat - Rxchat * chat);
        x0 = x0 + xhat;
        c0 = c0 + chat;

        % Update convergence criteria
        dxchat_last = dxchat_current;
        dxchat_current = [xhat; chat];
        meas_RMS = rms(postfit_residuals(~isnan(postfit_residuals)));
        meas_RMS_hist = [meas_RMS_hist meas_RMS];

        % Print progress
        fprintf('\b\b\b\b%3d%%', round((iteration_num / max_iter) * 100));
    end

    fprintf('\nPESF Completed!\n');

    % Store results in a struct
    results = struct( ...
        'x0', x0, ...
        'c0', c0, ...
        'p_all', p_all, ...
        'xhat_batch', xhat_batch, ...
        'chat_batch', chat_batch, ...
        'P_xx_batch', P_xx_batch, ...
        'P_xc_batch', P_xc_batch, ...
        'P_cc_batch', P_cc_batch, ...
        'postfit_residuals', postfit_residuals, ...
        'meas_RMS_hist', meas_RMS_hist ...
    );

end
