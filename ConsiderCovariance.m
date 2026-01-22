%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         January 22, 2024
% CLASS:        ASEN 6080: StatOD
% INSTRUCTOR:   Prof. Jay W. Mcmahon
% ASSIGNMENT:   Homework 7 - J3 Mismatch, CCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close;

%% Create True and Perturbed Trajectories

% Constants
coeffs_true = [0, 1.0826269 * 1e-3, -2.5324 * 1e-6]; % Zonal harmonic coefficients 
coeffs_ref = [0, 1.0826269 * 1e-3,  -2.5324 * 1e-6]; 
R = 6378; % Earth radius 
mu = 398600.4415;  % Gravitational parameter 
l_max_true = 3; % Zonal Harmonics Terms 
l_max_ref = 3;

% Orbital parameters 
a = 1e4;              % Semi-major axis (km)
e = 0.001;            % Eccentricity (-)
i = deg2rad(40);      % Inclination (rad)
RAAN = deg2rad(80);   % Right Ascension of Ascending Node (rad)
omega = deg2rad(40);  % Argument of Perigee (rad)
nu0 = deg2rad(0);     % True Anomaly (rad)

% Initial state (reference orbit)
state0 = orbitalElementsToCartesian(mu, a, e, i, RAAN, omega, nu0);
state0_true = [state0; mu; coeffs_true'; reshape(eye(10), [], 1)]; 
state0_ref = [state0; mu; coeffs_ref'; reshape(eye(10), [], 1)]; 


% Perturbation vector
dr0 = 1e-3; % Initial position deviation (km)
dv0 = 1e-6; % Initial velocity deviation (km/sec)
state_perturbed0_ref = state0_ref + [dr0; 0; 0; 0; dv0; 0; 0; 0; 0; 0; zeros(100, 1)];

% Initialize ODEs
[acc_func_true, A_func_true] = precomputeZonalSPH(R, l_max_true);
[acc_func_ref, A_func_ref] = precomputeZonalSPH(R, l_max_ref);

% Simulation parameters
T_orbit = 2 * pi * sqrt(a^3 / mu); % Orbital period (s)
t_span = 0:10:(15 * T_orbit);        % Simulate for 15 orbits, every 10 seconds
options = odeset('RelTol', 2.22045 * 1e-14, 'AbsTol',  2.22045 * 1e-14);

% Integrate reference trajectory
[t_true, state_true] = ode45(@(t, state) ...
     zonalSPH_ODE(t, state, coeffs_true, mu, l_max_true, acc_func_true, A_func_true), ...
     t_span, state0_true, options);

% Integrate perturbed trajectory with STM
[t_ref, state_ref] = ode45(@(t, state) ...
     zonalSPH_ODE(t, state, coeffs_ref, mu, l_max_ref, acc_func_ref, A_func_ref), ...
     t_span, state_perturbed0_ref, options);

% Plot Orbit
gca1 = figure(1);
plot3(state_true(:, 1), state_true(:, 2), state_true(:, 3), 'b', 'LineWidth', 2);  
hold on;  
plot3(state_ref(:, 1), state_ref(:, 2), state_ref(:, 3), ...
    'r', 'LineWidth', 2);   
plot_earth();  
axis equal;
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('z (km)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Trajectory True', 'Trajectory Reference (Perturbed)', 'Earth', ...
    'Location', 'northwest', 'Interpreter', 'latex','FontSize', 14);
grid on;
exportgraphics(gca1, 'CCA_Orbit_trajectories.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Extract True and Reference States, and STM Matrices

% Extract the STM matrices for [r, v]
STM_true = reshape(state_true(:, 11:end), [], 10, 10);
STM_true = STM_true(:, 1:6, 1:6);
STM_ref = reshape(state_ref(:, 11:end), [], 10, 10);
STM_ref  = STM_ref(:,  [1:6, 10], [1:6, 10]);

% Preallocate structure arrays
num_steps_true = length(t_true);
num_steps_ref = length(t_ref);
trajectory_true = struct('time', cell(1, num_steps_ref), 'state', ...
    cell(1, num_steps_ref), 'STM', cell(1, num_steps_ref), ...
    'parameters', cell(1, num_steps_ref),'function_handle', cell(1, num_steps_ref));
trajectory_ref = struct('time', cell(1, num_steps_ref), 'state', ...
    cell(1, num_steps_ref), 'STM', cell(1, num_steps_ref), ...
    'parameters', cell(1, num_steps_ref),'function_handle', cell(1, num_steps_ref));

% Fill the true trajectory structure
for i = 1:num_steps_true
    trajectory_true(i).time = t_true(i);
    trajectory_true(i).state = state_true(i, 1:6); % Extract r, v
    trajectory_true(i).parameters = state_true(i, 7:10); % Extract mu, J1, J2
    trajectory_true(i).STM = STM_true(i, :, :);
    trajectory_true(i).function_handle = @(t, x) zonalSPH_ODE(t, ...
        x, coeffs_true, mu, l_max_true, acc_func_true, A_func_true);
end

% Fill the reference trajectory structure
for i = 1:num_steps_ref
    trajectory_ref(i).time = t_ref(i);
    trajectory_ref(i).state = state_ref(i, 1:6); % Extract r, v
    trajectory_ref(i).parameters = state_ref(i, 7:10); % Extract mu, J1, J2
    trajectory_ref(i).STM = STM_ref(i, :, :);
    trajectory_ref(i).function_handle = @(t, x) zonalSPH_ODE(t, ...
        x, coeffs_ref, mu, l_max_ref, acc_func_ref, A_func_ref);
end

%% Create Measurements 

% Constants
earth_radius = 6378; % Earth's mean radius in km
omega_earth = 2 * pi / (24 * 3600); % Earth's rotation rate (rad/s)
theta0 = 122; % Initial Earth rotation angle in degrees
elevation_mask = 10; % Elevation mask in degrees
fT_ref = 8.44 * 1e9; % Transmit frequency in Hz (X-Band)
c = 299792.458; % Speed of light in km/s
sigma_rho = 1e-3; % Range Std in km
sigma_rho_dot = 1e-6; % Range Rate Std in km/sec
R = [sigma_rho^2, sigma_rho_dot^2]; % Range, Range-rate covariance in km, km/s

% Ground station locations (lat, lon in degrees)
stations = [
    -35.398333, 148.981944;  % Station 1
    40.427222, 355.749444;   % Station 2
    35.247164, 243.205       % Station 3
];

% Extract variables
sc_true = state_true(:, 1:6)'; 
sc_ref = state_ref(:, 1:6)';
times = t_ref; % Time array

% Preallocate cell array for station measurements
measurements_cell = cell(size(stations, 1), 1);

% Process each ground station
for station_idx = 1:size(stations, 1)
    % Compute measurements for this station
    measurements_cell{station_idx} = radiometric_measurement(sc_true, sc_ref, ...
        stations(station_idx, :), theta0, times, elevation_mask, R);
end

% Merge and sort measurements across all stations
station_names = {'GS1', 'GS2', 'GS3'};
sorted_measurements = merge_and_sort_measurements(measurements_cell, station_names);

%% Plot Measurements Residuals

% Extract Observed Measurements
[range_residuals, range_rate_residuals, time_per_station] = extract_measurements(sorted_measurements, ...
    stations, 'residual');

% Plot range and range-rate measurements
gca2 = figure(2);
set(gcf, 'Position', [100, 100, 1400, 900]);
subplot(2, 1, 1);
hold on;
for station_idx = 1:size(stations, 1)
    if ~isempty(range_residuals{station_idx}) && ...
       length(time_per_station{station_idx}) == length(range_residuals{station_idx})
        plot(time_per_station{station_idx} / 3600, range_residuals{station_idx}, ...
             '.', 'DisplayName', sprintf('GS %d', station_idx), 'MarkerSize', 10);
    end
end
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Range Residual (km)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
grid on;
hold off;
subplot(2, 1, 2);
hold on;
for station_idx = 1:size(stations, 1)
    if ~isempty(range_rate_residuals{station_idx}) && ...
       length(time_per_station{station_idx}) == length(range_rate_residuals{station_idx})
        plot(time_per_station{station_idx} / 3600, range_rate_residuals{station_idx}, ...
             '.', 'DisplayName', sprintf('GS %d', station_idx), 'MarkerSize', 10);
    end
end
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Range-Rate Residual (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 8, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(gca2, 'CCA_Radiometric_Residuals.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Esimation - LKF

% Initialize CCA
J3 =  -2.5324 * 1e-6;
c_bar =  J3; % Remove J3 state estimate
P_cc = J3^2; % Consider Covariance

% Define initial covariance matrix P0
% sigma_pos = 1;       % Std for position states (km)
% sigma_vel = 1e-3;    % Std for velocity states (km/s)
% P0 = diag([sigma_pos^2, sigma_pos^2, sigma_pos^2, ...
%            sigma_vel^2, sigma_vel^2, sigma_vel^2]);
P0 = 1e4 * eye(6);

% Run the Filter (LKF w/ CCA)
results = lkf_consider(trajectory_ref, sorted_measurements, P0, c_bar, P_cc);
filename_suffix = 'LKF';

%% Plot filtering results

% Initialize
measurement_times = [sorted_measurements.time];
trajectory_times = [trajectory_true.time];
[~, traj_indices] = ismember(measurement_times, trajectory_times);

% Extract true states at exact measurement times
true_states_at_meas = zeros(6, length(measurement_times));
for i = 1:length(measurement_times)
    true_states_at_meas(:, i) = trajectory_true(traj_indices(i)).state(1:6);
end

% Compute state error (only at measurement times)
state_error = true_states_at_meas - results.state_corrected_hist;

% Compute 3-sig bounds from covariance matrix at measurement times
sigma_bounds_meas = zeros(size(state_error));
T_meas = length(measurement_times);
for t = 1:T_meas
    sigma_bounds_meas(:, t) = 3 * sqrt(diag(results.P_hist(1:6, 1:6, t)));
end

% Convert measurement times from seconds to hours
measurement_times_hours = measurement_times / 3600; 

% Plot state error
gca4 = figure(4);
set(gcf, 'Position', [100, 100, 1400, 900]);
colors = lines(6);
labels = {'$x$', '$y$', '$z$', '$\dot{x}$', '$\dot{y}$', '$\dot{z}$'};
y_labels = {'Position Error (km)', 'Velocity Error (km/s)'};
for i = 1:6
    subplot(2,3,i);
    hold on;
    scatter(measurement_times_hours, state_error(i, :), 8, colors(i, :),...
        'filled', 'DisplayName', labels{i}, 'MarkerEdgeColor', 'k'); 
    scatter(measurement_times_hours, 2 * sigma_bounds_meas(i, :), 5, ...
        colors(i, :), 'filled', 'DisplayName', '$\pm2\sigma$');
    scatter(measurement_times_hours, -3 * sigma_bounds_meas(i, :), 5, ...
        colors(i, :), 'filled', 'HandleVisibility', 'off');
    xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(y_labels{(i > 3) + 1}, 'Interpreter', 'latex', 'FontSize', 14);
    legend('show', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    hold off;
end
exportgraphics(gca4, sprintf('CCA_StateErrors_%s.pdf', filename_suffix), ...,
    'ContentType','image', 'Resolution', 1000);

% Plot state errors in semilogy
gca5 = figure(5);
set(gcf, 'Position', [100, 100, 1400, 900]);
for i = 1:6
    subplot(2,3,i);    
    sigma_bounds_meas(i, sigma_bounds_meas(i, :) == 0) = 1e-14; 
    semilogy(measurement_times_hours, abs(state_error(i, :)), ...
        'o', 'MarkerSize', 5, 'MarkerFaceColor', colors(i, :), ...
        'MarkerEdgeColor', 'k', 'DisplayName', labels{i});
    hold on;
    scatter(measurement_times_hours, 2 * sigma_bounds_meas(i, :), 5, ...
        colors(i, :), 'filled', 'DisplayName', '$+2\sigma$');
    xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(y_labels{(i > 3) + 1}, 'Interpreter', 'latex', 'FontSize', 14);
    legend('show', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    hold off;
end
exportgraphics(gca5, sprintf('CCA_StateErrorsLog_%s.pdf', filename_suffix), ...,
    'ContentType','image', 'Resolution', 1000);

% Plot range and range-rate post-fit residuals 
gca6 = figure(6);
set(gcf, 'Position', [100, 100, 1400, 900]);
subplot(2,1,1);
hold on;
scatter(measurement_times_hours, results.postfit_residuals(1, :), 8, 'b',...
    'filled', 'DisplayName', 'Range Residual', 'MarkerEdgeColor', 'k'); 
fill([measurement_times_hours, fliplr(measurement_times_hours)], ...
    [3 * sigma_rho * ones(size(measurement_times_hours)), ...
    fliplr(-3 * sigma_rho * ones(size(measurement_times_hours)))], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '$\pm3\sigma$');
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Range Residual (km)', 'Interpreter', 'latex', 'FontSize', 14);
legend('show', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
hold off;
ylim([-10 * sigma_rho, 10 * sigma_rho]);
subplot(2,1,2);
hold on;  
scatter(measurement_times_hours, results.postfit_residuals(2, :), 8, ...
    'r', 'filled', 'DisplayName', 'Range Rate Residual', 'MarkerEdgeColor', 'k'); 
fill([measurement_times_hours, fliplr(measurement_times_hours)], ...
    [3 * sigma_rho_dot * ones(size(measurement_times_hours)), ...
    fliplr(-3 * sigma_rho_dot * ones(size(measurement_times_hours)))], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '$\pm3\sigma$');
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14, 'FontSize', 14);
ylabel('Range-Rate Residual (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
legend('show', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
hold off;
ylim([-5 * sigma_rho_dot, 5 * sigma_rho_dot]);
exportgraphics(gca6, sprintf('CCA_PostFit_%s.pdf', filename_suffix), ...,
    'ContentType','image', 'Resolution', 1000);

%% Print Stats for HW7

% RMS Errors (Position & Velocity)
pos_err = state_error(1:3, :);
vel_err = state_error(4:6, :);
rms_pos = sqrt(mean(vecnorm(pos_err, 2, 1).^2));  % norm at each time, then RMS
rms_vel = sqrt(mean(vecnorm(vel_err, 2, 1).^2));
fprintf('RMS Position Error [km]: %.14f\n', rms_pos);
fprintf('RMS Velocity Error [km/s]: %.14f\n', rms_vel);

% Percent of States Outside 2sig Bounds 
two_sigma_bounds = 2 * sqrt(diag(results.P_hist(:,:,1)));  % for dimension check
outside_count = 0;
total_components = numel(state_error);
for t = 1:T_meas
    for i = 1:6
        if abs(state_error(i, t)) > 2 * sqrt(results.P_hist(i,i,t))
            outside_count = outside_count + 1;
        end
    end
end
percent_outside = 100 * outside_count / total_components;
fprintf('%% of state elements outside 2-sigma bounds: %.2f%%\n', percent_outside);

%% Map final state correction and covariance back to epoch

% Final filter estimate (deviation and covariance)
dx_f = results.state_deviation_hist(:, end);           % Includes consider correction
Pcf = results.P_full_hist(:, :, end);                  % Consider-inflated posterior

% Get STM(t_f, t_0) from trajectory_ref
STM_tf_t0 = squeeze(trajectory_ref(end).STM);          
STM_t0_tf = inv(STM_tf_t0);                           

% Split STM_t0_tf
Phi = STM_t0_tf(1:6, 1:6);        % Partial wrt state
Psi = STM_t0_tf;

% Map deviation back
dx0_cons = Psi * [dx_f; J3];
x0_est = trajectory_ref(1).state' - dx0_cons(1:6);

% Total posterior at t0
P0_cons = Psi * Pcf * Psi';

% Re-integrate from new estimate
coeffs_ref = [0, 1.0826269 * 1e-3,  0];
x0_aug = [x0_est; mu; coeffs_ref'; reshape(eye(10), [], 1)];
[t_updated, state_updated] = ode113(@(t, state) ...
    zonalSPH_ODE(t, state, coeffs_ref, mu, l_max_ref, acc_func_ref, A_func_ref), ...
    t_span, x0_aug, options);

% Extract STM for [r,v,J3]
STM_updated = reshape(state_updated(:, 11:end), [], 10, 10);
STM_updated = STM_updated(:, [1:6, 10], [1:6, 10]);  % Includes consider direction

% Propagate Pcc_0 forward using STM
N = length(t_updated);
state_errors_updated = zeros(6, N);
Pc_hist = zeros(7, 7, N);
for i = 1:N
    % Propagate Covariance
    Psi_i = squeeze(STM_updated(i, :, :));   
    Pc_hist(:, :, i) = Psi_i * P0_cons * Psi_i';

    % Error wrt true trajectory
    state_errors_updated(:, i) = trajectory_true(i).state' - state_updated(i, 1:6)';
end

% Plot
gca7 = figure(7);
set(gcf, 'Position', [100, 100, 1400, 900]);
colors = lines(6);
labels = {'$x$', '$y$', '$z$', '$\dot{x}$', '$\dot{y}$', '$\dot{z}$'};
y_labels = {'Position Error (km)', 'Velocity Error (km/s)'};
t_updated_hr = t_updated / 3600;  % convert to hours
for i = 1:6
    subplot(2,3,i);
    set(gca, 'YScale', 'log'); 
    hold on;
    scatter(t_updated_hr, abs(state_errors_updated(i, :)), 8, colors(i, :), ...
        'filled', 'DisplayName', labels{i}, 'MarkerEdgeColor', 'k');
    two_sigma = 2 * sqrt(squeeze(Pc_hist(i, i, :)))';
    scatter(t_updated_hr, +two_sigma, 5, colors(i, :), ...
        'filled', 'DisplayName', '$\pm2\sigma$');
    xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(y_labels{(i > 3) + 1}, 'Interpreter', 'latex', 'FontSize', 14);
    legend('show', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    hold off;
end
exportgraphics(gca7, sprintf('CCA_BackForPropLog_%s.pdf', filename_suffix), ...,
    'ContentType','image', 'Resolution', 1000);

% Compute running RMS norms over time
N = size(state_errors_updated, 2);
rms_pos_time = zeros(1, N);
rms_vel_time = zeros(1, N);
for k = 1:N
    pos_norms = vecnorm(state_errors_updated(1:3, 1:k), 2, 1);
    vel_norms = vecnorm(state_errors_updated(4:6, 1:k), 2, 1);
    rms_pos_time(k) = sqrt(mean(pos_norms.^2));
    rms_vel_time(k) = sqrt(mean(vel_norms.^2));
end

% Plot
gca9 = figure(9);
set(gcf, 'Position', [100, 100, 1400, 600]);
t_hr = t_updated / 3600;
colors = lines(2);
semilogy(t_hr, rms_pos_time, '--', ...
    'MarkerSize', 4, 'LineWidth', 5, ...
    'Color', colors(1,:), ...
    'DisplayName', 'Position RMS');
hold on;
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error (km, km/sec)', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
semilogy(t_hr, rms_vel_time, '--', ...
    'MarkerSize', 4, 'LineWidth', 5, ...
    'Color', colors(2,:), ...
    'DisplayName', 'Velocity RMS');
legend('show', 'Interpreter', 'latex', 'FontSize', 14, 'location', 'northwest');
grid on;
hold off;
exportgraphics(gca9, sprintf('CCA_BackForProp_RMSt_%s.pdf', filename_suffix), ...,
    'ContentType','image', 'Resolution', 1000);

%% TEST GOF
% 
% %  Mahalanobis Distance Squared for Full State Error (6 DOF) 
% dof = 6;
% confidence_level = 0.997;  % Equivalent to ~3σ
% chi2_threshold = chi2inv(confidence_level, dof);  % 99.7% chi-squared bound
% 
% chi2_dist = zeros(1, N);  % Mahalanobis D^2 over time
% 
% for k = 1:length(measurement_times)
%     e_k = state_errors_updated(1:6, k);          % Full 6x1 error
%     P_k = Pc_hist(1:6, 1:6, k);                  % Corresponding 6x6 covariance
%     chi2_dist(k) = e_k' / P_k * e_k;             % Mahalanobis D^2
% end
% 
% %  Perform GOF Test 
% frac_outside = mean(chi2_dist > chi2_threshold);
% fprintf('%.2f%% of state errors lie outside %.1f-DOF chi-squared 99.7%% bound (%.2f)\n', ...
%         100 * frac_outside, dof, chi2_threshold);
% 
% %  Plot Mahalanobis Distance Squared 
% gca10 = figure(10); clf;
% set(gca10, 'Position', [100, 100, 1200, 600]);
% semilogy(t_updated / 3600, chi2_dist, 'LineWidth', 2, ...
%     'DisplayName', '$D^2$ Mahalanobis');
% hold on;
% yline(chi2_threshold, 'r--', 'LineWidth', 2, ...
%     'DisplayName', '$D^2$ 6 DOF', ...
%     'Interpreter', 'latex', 'FontSize', 14);
% xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('$D^2 = \mathbf{e}^\top \mathbf{P}^{-1} \mathbf{e}$', ...
%     'Interpreter', 'latex', 'FontSize', 14);
% title('Mahalanobis Distance Squared: Consistency Test', ...
%     'Interpreter', 'latex', 'FontSize', 16);
% legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest');
% grid on;
% 
% %  Export plot 
% exportgraphics(gca10, sprintf('CCA_Mahalanobis_%s.pdf', filename_suffix), ...
%     'ContentType', 'image', 'Resolution', 1000);



%% Helper Functions Dynamics

function R = DCM_eciToEcef(t)
     % Earth's rotation rate 
    omega_earth = 2 * pi / (24 * 3600);
    
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

function [acceleration, A] = zonalSPH_Dynamics(position_eci, coeffs_vals, mu_val, l_max, t, acc_func, A_func)
    % Convert ECI to ECEF coordinates
    R_eci_to_ecef = DCM_eciToEcef(t);               % Rotation matrix ECI to ECEF
    position_ecef = R_eci_to_ecef * position_eci;   % Transform position to ECEF
 
    % Evaluate acceleration and parameter derivatives at the given inputs
    acceleration = acc_func(position_ecef(1), position_ecef(2), ...
                                     position_ecef(3), mu_val, coeffs_vals');

    % Convert acceleration back to ECI
    acceleration(1:3) = R_eci_to_ecef' * acceleration(1:3); 

    if nargout > 1
        % Evaluate Jacobian at the given position and parameters
        A_ecef = A_func(position_ecef(1), position_ecef(2), position_ecef(3), ...
                        0, 0, 0, mu_val, coeffs_vals');

        % Transform partial derivatives back to ECI
        Ar_eci = R_eci_to_ecef' * A_ecef(4:6, 1:3) * R_eci_to_ecef; 
        Amu_eci = R_eci_to_ecef' * A_ecef(4:6, 7);
        Acl_eci = R_eci_to_ecef' * A_ecef(4:6, 7:end);

        % Construct the full state transition matrix A
        A = [zeros(3, 3), eye(3), zeros(3, 1), zeros(3, l_max);            % Velocity derivatives
             Ar_eci, zeros(3, 3), Amu_eci, Acl_eci(:, 2:end);              % Acceleration derivatives
             zeros(1, 6), 0, zeros(1, l_max);                              % Mu derivatives
             zeros(l_max, 7), zeros(l_max, l_max)];                        % Harmonics derivatives
    end
end

function dstate = zonalSPH_ODE(t, state, coeffs, mu, l_max, acc_func, A_func)
    % Extract from state
    position_eci = state(1:3);
    velocity_eci = state(4:6);

    % Compute acceleration and A matrix
    [acceleration_eci, A_eci] = zonalSPH_Dynamics(position_eci, coeffs, mu, ...
         l_max, t, acc_func, A_func);

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

function measurements_struct = radiometric_measurement(sc_true, sc_ref, ...
    station, theta0, times, elevation_mask, R)
    % Compute ground station state in ECI
    [gs_pos, gs_vel] = compute_gs_eci(station, theta0, times);

    % Compute true and reference spacecraft state in ECI
    sc_true_pos = sc_true(1:3, :);
    sc_true_vel = sc_true(4:6, :);
    sc_ref_pos = sc_ref(1:3, :);
    sc_ref_vel = sc_ref(4:6, :);

    % Extract sigma
    sigma_rho = sqrt(R(1));
    sigma_rho_dot = sqrt(R(2));
    
    % Compute visibility (on true)
    [is_visible, ~] = compute_visibility_mask(sc_true_pos, gs_pos, elevation_mask);

    % Preallocate measurements struct
    N = length(times);
    measurements_struct = struct('time', cell(1, N), 'measurement', cell(1, N), ...
                                 'partials', cell(1, N), 'visibility', cell(1, N));

    % Loop through each time step
    for i = 1:N
        % True relative position and velocity
        dR_true = sc_true_pos(:, i) - gs_pos(:, i);
        dV_true = sc_true_vel(:, i) - gs_vel(:, i);

        % Reference relative position and velocity
        dR_ref = sc_ref_pos(:, i) - gs_pos(:, i);
        dV_ref = sc_ref_vel(:, i) - gs_vel(:, i);

        % Observed (on true) range and range-rate
        rho_obs = norm(dR_true) + sigma_rho * randn;
        rho_dot_obs = dot(dR_true, dV_true) / rho_obs + sigma_rho_dot * randn;

        % Computed (on reference) range and range-rate
        rho_comp = norm(dR_ref);
        rho_dot_comp = dot(dR_ref, dV_ref) / rho_comp;

        % Residuals
        rho_res = rho_obs - rho_comp;
        rho_dot_res = rho_dot_obs - rho_dot_comp;

        % Compute partial derivatives (on reference)
        d_rho_dR = dR_ref' / rho_comp; 
        d_rho_dV = zeros(1, 3);
        d_rho_dot_dR = (rho_comp * dV_ref' - rho_dot_comp * dR_ref') / ...
            rho_comp^2;
        d_rho_dot_dV = dR_ref' / rho_comp;

        % Spacecraft state Jacobian
        d_measurements_dX = [d_rho_dR, d_rho_dV;
                             d_rho_dot_dR, d_rho_dot_dV];

        % Ground station position Jacobian
        d_rho_dCCA = 0;
        d_rho_dot_dCCA = 0;

        % Ground station state Jacobian
        d_measurements_dCCA = [d_rho_dCCA;
                             d_rho_dot_dCCA];

        % Store function handle for measurement model
        measurement_func = @(x) measurement_function(x, gs_pos(:, i), gs_vel(:, i));

        % Store results in structure
        measurements_struct(i).time = times(i);
        measurements_struct(i).observed = [rho_obs; rho_dot_obs];
        measurements_struct(i).computed = [rho_comp; rho_dot_comp];
        measurements_struct(i).residual = [rho_res; rho_dot_res];
        measurements_struct(i).partials = struct('wrt_X', d_measurements_dX, ...
                                                 'wrt_CCA', d_measurements_dCCA);
        measurements_struct(i).covariance = R;
        measurements_struct(i).visibility = is_visible(i);
        measurements_struct(i).measurement_function = measurement_func;
    end
end

function [h_x, H] = measurement_function(x, gs_pos, gs_vel)
    % Compute range and range-rate
    dR = x(1:3) - gs_pos;
    dV = x(4:6) - gs_vel;
    rho = norm(dR); % Range
    rho_dot = dot(dR, dV) / rho; % Range rate
    h_x = [rho; rho_dot];

    % Compute Jacobian
    d_rho_dR = dR' / rho;
    d_rho_dV = zeros(1, 3);
    d_rho_dot_dR = (rho * dV' - (dot(dR, dV) / rho) * dR') / rho^2;
    d_rho_dot_dV = dR' / rho;

    % Full Jacobian
    H = [d_rho_dR, d_rho_dV;
         d_rho_dot_dR, d_rho_dot_dV];
end


function sorted_measurements = merge_and_sort_measurements(measurements_cell, station_names)
    % Flatten all measurement structures into a single array
    all_measurements = [];
    
    for i = 1:length(measurements_cell)
        station_measurements = measurements_cell{i};
        for j = 1:length(station_measurements)
            if station_measurements(j).visibility
                % Attach station name directly into struct
                station_measurements(j).station_name = station_names{i}; 
                all_measurements = [all_measurements, station_measurements(j)]; %#ok<AGROW>
            end
        end
    end

    % Sort by time
    [~, sorted_indices] = sort([all_measurements.time]);
    all_measurements = all_measurements(sorted_indices);

    % Group by unique time instances
    unique_times = unique([all_measurements.time]);
    N = length(unique_times);
    sorted_measurements = struct('time', cell(1, N), 'observed', cell(1, N), ...
                                 'computed', cell(1, N), 'residual', cell(1, N), ...
                                 'partials', cell(1, N), 'covariance', cell(1, N), ...
                                 'function_handle', cell(1, N), ...
                                 'station_labels', cell(1, N));

    for i = 1:N
        t = unique_times(i);
        matching_indices = find([all_measurements.time] == t);

        % Stack measurements and partials
        stacked_observed = [];
        stacked_computed = [];
        stacked_residual = [];
        stacked_partials_X = [];
        stacked_partials_CCA = [];
        stacked_covariance = [];
        stacked_measurement_func = [];
        stacked_labels = {};

        for j = matching_indices
            stacked_observed = [stacked_observed; ...
                all_measurements(j).observed];  %#ok<AGROW> 
            stacked_computed = [stacked_computed; ...
                all_measurements(j).computed];  %#ok<AGROW>
            stacked_residual = [stacked_residual; ...
                all_measurements(j).residual];  %#ok<AGROW>
            stacked_partials_X = [stacked_partials_X; ...
                all_measurements(j).partials.wrt_X];  %#ok<AGROW>
            stacked_partials_CCA = [stacked_partials_CCA; ...
                all_measurements(j).partials.wrt_CCA];  %#ok<AGROW>
            stacked_covariance = [stacked_covariance; ...
                all_measurements(j).covariance]; %#ok<AGROW>
            stacked_measurement_func = [stacked_measurement_func; ...
                all_measurements(j).measurement_function]; %#ok<AGROW>
            stacked_labels = [stacked_labels; 
                              sprintf('%s - Range', all_measurements(j).station_name);
                              sprintf('%s - Range Rate', all_measurements(j).station_name)];  %#ok<AGROW>
        end

        % Store in struct
        sorted_measurements(i).time = t;
        sorted_measurements(i).observed = stacked_observed;
        sorted_measurements(i).computed = stacked_computed;
        sorted_measurements(i).residual = stacked_residual;
        sorted_measurements(i).partials = struct('wrt_X', stacked_partials_X, ...
            'wrt_CCA', stacked_partials_CCA);
        sorted_measurements(i).covariance = stacked_covariance;
        sorted_measurements(i).function_handle = stacked_measurement_func;
        sorted_measurements(i).station_labels = stacked_labels;
    end
end


function [gs_positions_eci, gs_velocities_eci] = compute_gs_eci(gs_location, theta0, times)
    % Constants
    earth_radius = 6378; 
    omega_earth = 2 * pi / (24 * 3600); 

    % Convert lat/long to radians
    lat = deg2rad(gs_location(1));
    lon = deg2rad(gs_location(2));
    theta0_rad = deg2rad(theta0);

    % Pre-allocate outputs
    N = length(times);
    gs_positions_eci = zeros(3, N);
    gs_velocities_eci = zeros(3, N);

    % Loop through times to compute ECI state
    for i = 1:N
        t = times(i);
        % Earth's rotation angle at time t
        theta = theta0_rad + omega_earth * t;

        % Ground station position in ECI
        gs_positions_eci(:, i) = earth_radius * [
            cos(lon + theta) * cos(lat); 
            sin(lon + theta) * cos(lat); 
            sin(lat)
        ];

        % Ground station velocity in ECI
        gs_velocities_eci(:, i) = omega_earth * cross([0; 0; 1], ...
            gs_positions_eci(:, i));
    end
end

function [is_visible, elevation_angles] = compute_visibility_mask(sc_positions, ...
    gs_positions, elevation_mask)
    % Number of time steps
    num_steps = size(sc_positions, 2);
    elevation_angles = zeros(1, num_steps);

    % Pre-allocate visibility mask
    is_visible = false(1, num_steps);

    % Loop through each time step
    for i = 1:num_steps
        % Relative position vector 
        relative_position = sc_positions(:, i) - gs_positions(:, i);

        % Elevation angle computation 
        % I.E., projection of relative position onto local zenith direction
        zenith_direction = gs_positions(:, i) / norm(gs_positions(:, i));
        sin_elevation = dot(relative_position, zenith_direction) / norm(relative_position);
        elevation_angles(i) = asind(sin_elevation); % Convert to degrees

        % Check if elevation angle exceeds the elevation mask
        if elevation_angles(i) > elevation_mask
            is_visible(i) = true;
        end
    end
end

function state = orbitalElementsToCartesian(mu, a, e, i, RAAN, omega, nu)
    % Compute the distance (r) and speed (v) in the orbital plane
    p = a * (1 - e^2); % Semi-latus rectum
    r_orb = p / (1 + e * cos(nu)); % Radius in the orbital plane
    v_orb = sqrt(mu / p); % Velocity magnitude in the orbital plane

    % State in the orbital plane (PQW frame)
    r_PQW = r_orb * [cos(nu); sin(nu); 0];
    v_PQW = v_orb * [-sin(nu); e + cos(nu); 0];

    % Rotation matrices
    R_RAAN = [cos(RAAN), -sin(RAAN), 0;
              sin(RAAN),  cos(RAAN), 0;
              0,          0,         1];
    R_i = [1, 0,       0;
           0, cos(i), -sin(i);
           0, sin(i),  cos(i)];
    R_omega = [cos(omega), -sin(omega), 0;
               sin(omega),  cos(omega), 0;
               0,           0,          1];

    % Combined rotation matrix
    R = R_RAAN * R_i * R_omega;

    % Transform position and velocity to ECI frame
    state = [R * r_PQW; R * v_PQW];
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

function [range_measurements, range_rate_measurements, time_per_station] = extract_measurements(sorted_measurements, ...
    stations, measurement_type)

    % Extract range and range-rate measurements from sorted structure
    measurement_times = [sorted_measurements.time];
    num_measurements = length(sorted_measurements);
    
    % Preallocate cell arrays
    range_measurements = cell(size(stations, 1), 1);
    range_rate_measurements = cell(size(stations, 1), 1);
    time_per_station = cell(size(stations, 1), 1); % Store time for each GS
    
    % Extract measurements for each ground station
    for station_idx = 1:size(stations, 1)
        range_measurements{station_idx} = [];
        range_rate_measurements{station_idx} = [];
        time_per_station{station_idx} = [];
    
        for i = 1:num_measurements
            % Find the indices corresponding to this ground station
            gs_mask = contains(sorted_measurements(i).station_labels, sprintf('GS%d', station_idx));
    
            % Extract range and range-rate ensuring proper pairing
            if any(gs_mask) && length(sorted_measurements(i).(measurement_type)) >= sum(gs_mask)
                % Get the first occurrence of this GS's range measurement
                idx = find(gs_mask, 1, 'first');
    
                % Ensure index does not exceed array bounds
                if idx + 1 <= length(sorted_measurements(i).(measurement_type))
                    range_measurements{station_idx} = [range_measurements{station_idx}; ...
                                                       sorted_measurements(i).(measurement_type)(idx)];
                    range_rate_measurements{station_idx} = [range_rate_measurements{station_idx}; ...
                                                            sorted_measurements(i).(measurement_type)(idx + 1)];
                    time_per_station{station_idx} = [time_per_station{station_idx}; measurement_times(i)];
                end
            end
        end
    end
end

function results = lkf(trajectory_ref, sorted_measurements, P0)
    % Initialization
    T = length(sorted_measurements); % Number of measurements
    n = size(P0, 1);  % State dimension
    max_m = max(arrayfun(@(s) length(s.residual), sorted_measurements));

    state_deviation_hist = zeros(n, T);
    state_corrected_hist = zeros(n, T);
    P_hist = zeros(n, n, T);
    postfit_residuals = NaN(max_m, T);  % Preallocate with NaN

    % Initial state deviation and covariance
    dx = zeros(n, 1); % Initial deviation
    P = P0;           % Initial covariance

    % Identity matrix for stability
    I = eye(n);
    STM_tm = I;  % Initial transition matrix (identity)

    % Extract trajectory times for reference
    trajectory_times = [trajectory_ref.time];

    % Find indices of trajectory points that match measurement times
    measurement_times = [sorted_measurements.time];
    [~, traj_indices] = ismember(measurement_times, trajectory_times);

    % Progress tracking
    fprintf('Kalman Filter Progress: 0%%');

    % Kalman Filter loop over measurements
    for t = 1:T
        % Extract corresponding STM and state from trajectory
        traj_idx = traj_indices(t);
        STM_t = squeeze(trajectory_ref(traj_idx).STM);
        x_ref = trajectory_ref(traj_idx).state(1:6)';

        % Compute STM transition matrix
        STM_dt = STM_t / STM_tm; 

        % Extract measurement data
        prefit_residual = sorted_measurements(t).residual;
        H_tilde = sorted_measurements(t).partials.wrt_X;
        R = diag(sorted_measurements(t).covariance);

        % Prediction Step
        dx_pred = STM_dt * dx;
        P_pred = STM_dt * P * STM_dt';

        % Compute pre-fit residual
        postfit_res = prefit_residual - H_tilde * dx_pred;

        % Kalman Gain computation
        S = H_tilde * P_pred * H_tilde' + R;
        K = P_pred * H_tilde' / S;

        % Update Step
        dx_upd = dx_pred + K * postfit_res;
        P_upd = (I - K * H_tilde) * P_pred * (I - K * H_tilde)' + K * R * K';

        % Store results
        state_deviation_hist(:, t) = dx_upd;
        state_corrected_hist(:, t) = x_ref + dx_upd;
        P_hist(:, :, t) = P_upd;
        postfit_residuals(1:length(prefit_residual), t) = postfit_res;

        % Update for next iteration
        dx = dx_upd;
        P = P_upd;
        STM_tm = STM_t; % Store STM for next transition computation

        % Print progress
        fprintf('\b\b\b\b%3d%%', round((t / T) * 100));
    end

    fprintf('\nKalman Filter Completed!\n');

    % Store all outputs in a struct (dictionary)
    results = struct(...
        'state_corrected_hist', state_corrected_hist, ...
        'state_deviation_hist', state_deviation_hist, ...
        'P_hist', P_hist, ...
        'postfit_residuals', postfit_residuals ...
    );
end

function results = lkf_consider(trajectory_ref, sorted_measurements, P0, c_bar, P_cc)
% TODO: check theta
    % Initialize 
    T = length(sorted_measurements);
    n = size(P0, 1);
    n_consider = length(c_bar);
    max_m = max(arrayfun(@(s) length(s.residual), sorted_measurements));

    state_dev_hist = zeros(n, T);
    state_corr_hist = zeros(n, T);
    P_hist = zeros(n, n, T);
    P_full_hist = zeros(n + n_consider, n + n_consider, T);
    postfit_residuals = NaN(max_m, T);

    dx = zeros(n, 1);     % initial deviation
    P = P0;

    STM_tm = eye(n + n_consider);
    trajectory_times = [trajectory_ref.time];
    measurement_times = [sorted_measurements.time];
    [~, traj_indices] = ismember(measurement_times, trajectory_times);

    % Initialize S0
    t0 = 1;
    H_x0 = sorted_measurements(t0).partials.wrt_X;
    H_c0 = sorted_measurements(t0).partials.wrt_CCA;
    R0 = diag(sorted_measurements(t0).covariance);
    S_meas0 = H_x0 * P0 * H_x0' + R0;
    K0 = P0 * H_x0' / S_meas0;
    S = -K0 * H_c0;

    fprintf('Consider Cov Filter Progress: 0%%');

    for t = 1:T
        traj_idx = traj_indices(t);
        STM_t = squeeze(trajectory_ref(traj_idx).STM);  % size: (n x (n + n_consider))
        x_ref = trajectory_ref(traj_idx).state(1:6)';    % reference state (without consider params)

        % Transition
        STM_dt = STM_t / STM_tm;
        Phi = STM_dt(1:n, 1:n);
        Theta = STM_dt(1:n, n+1:end);

        % Time updates
        dx_pred = Phi * dx;
        P_bar = Phi * P * Phi';
        S_bar = Phi * S + Theta;
        % x_bar_cons = x_bar + S_bar * c_bar;
        % P_bar_cons = P_bar + S_bar * P_cc * S_bar';
        % P_bar_cross = ...

        % Measurement model
        res = sorted_measurements(t).residual;
        H_x = sorted_measurements(t).partials.wrt_X;
        H_c = sorted_measurements(t).partials.wrt_CCA;
        R = diag(sorted_measurements(t).covariance);

        % Kalman gain
        S_meas = H_x * P_bar * H_x' + R;
        K = P_bar * H_x' / S_meas;

        % Post-fit residual
        y_pred = H_x * dx_pred;
        r_postfit = res - y_pred;

        % Measurement update
        dx_upd = dx_pred + K * r_postfit;
        dx_upd_cons = dx_upd +  S * c_bar;
        P = (eye(n) - K * H_x) * P_bar;
        S = (eye(n) - K * H_x) * S_bar - K * H_c;
        % x_corr = x_ref + dx_upd;
        x_corr_cons = x_ref + dx_upd_cons;
        P_cons = P + S * P_cc * S';
        P_cross = S * P_cc;

        % Store
        state_dev_hist(:, t) = dx_upd_cons;
        state_corr_hist(:, t) = x_corr_cons;
        P_hist(:, :, t) = P_cons;
        P_full_hist(:, :, t) = [P_cons, P_cross; P_cross', P_cc];
        postfit_residuals(1:length(res), t) = r_postfit;

        % Prepare for next step
        dx = dx_upd;
        STM_tm = STM_t;

        fprintf('\b\b\b\b%3d%%', round((t / T) * 100));
    end

    fprintf('\nConsider Covariance Filter Completed!\n');

    results = struct(...
        'state_corrected_hist', state_corr_hist, ...
        'state_deviation_hist', state_dev_hist, ...
        'P_hist', P_hist, ...
        'P_full_hist', P_full_hist, ...
        'postfit_residuals', postfit_residuals ...
    );
end
