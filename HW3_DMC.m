%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         January 22, 2024
% CLASS:        ASEN 6080: StatOD
% INSTRUCTOR:   Prof. Jay W. Mcmahon
% ASSIGNMENT:   Homework 3 - J3 Mismatch, DMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close;

%% Create True and Perturbed Trajectories

% Constants
coeffs_true = [0, 1.0826269 * 1e-3, -2.5324 * 1e-6]; % Zonal harmonic coefficients 
coeffs_ref = [0, 1.0826269 * 1e-3]; 
R = 6378; % Earth radius 
mu = 398600.4415;  % Gravitational parameter 
l_max_true = 3; % Zonal Harmonics Terms 
l_max_ref = 2;

% Orbital parameters 
a = 1e4;              % Semi-major axis (km)
e = 0.001;            % Eccentricity (-)
i = deg2rad(40);      % Inclination (rad)
RAAN = deg2rad(80);   % Right Ascension of Ascending Node (rad)
omega = deg2rad(40);  % Argument of Perigee (rad)
nu0 = deg2rad(0);     % True Anomaly (rad)

% Initial state (reference orbit)
state0 = orbitalElementsToCartesian(mu, a, e, i, RAAN, omega, nu0);
state0_true = [state0; zeros(3, 1); reshape(eye(9), [], 1)]; 
state0_ref = [state0; zeros(3, 1); reshape(eye(9), [], 1)]; 

% Perturbation vector
dr0 = 1e-3; % Initial position deviation (km)
dv0 = 1e-6; % Initial velocity deviation (km/sec)
state_perturbed0_ref = state0_ref;
state_perturbed0_ref(1) = state_perturbed0_ref(1) + dr0;
state_perturbed0_ref(4) = state_perturbed0_ref(4) + dv0;

% Simulation parameters
T_orbit = 2 * pi * sqrt(a^3 / mu);   % Orbital period (s)
t_span = 0:10:(15 * T_orbit);        % Simulate for 15 orbits, every 10 seconds
options = odeset('RelTol', 2.22045 * 1e-14, 'AbsTol',  2.22045 * 1e-14);

% Initialize Gauss-Markov Process
B_DMC = (T_orbit)^-1 * eye(3);

% Initialize ODEs
[acc_func_true, A_func_true] = precomputeZonalSPH(R, l_max_true);
[acc_func_ref, A_func_ref] = precomputeZonalSPH(R, l_max_ref);
[acc_DMC, A_DMC] = precomputeGaussMarkov(B_DMC); 

% Integrate reference trajectory
[t_true, state_true] = ode45(@(t, state) ...
     sc_ODE(t, state, coeffs_true, mu, l_max_true, acc_func_true, A_func_true, ...
     acc_DMC, A_DMC), t_span, state0_true, options);

% Integrate perturbed trajectory with STM
[t_ref, state_ref] = ode45(@(t, state) ...
     sc_ODE(t, state, coeffs_ref, mu, l_max_ref, acc_func_ref, A_func_ref, ...
     acc_DMC, A_DMC), t_span, state_perturbed0_ref, options);

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

%% Extract True and Reference States, and STM Matrices

% Extract the STM matrices for [r, v]
STM_true = reshape(state_true(:, 10:end), [], 9, 9);  
STM_ref = reshape(state_ref(:, 10:end), [], 9, 9);

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
    trajectory_true(i).state = state_true(i, 1:9); % Extract r, v
    trajectory_true(i).parameters = state_true(i, 7:9); % Extract w_dmc
    trajectory_true(i).STM = STM_true(i, :, :);
    trajectory_true(i).function_handle = @(t, x) sc_ODE(t, ...
        x, coeffs_true, mu, l_max_true, acc_func_true, A_func_true, acc_DMC, A_DMC);
end

% Fill the reference trajectory structure
for i = 1:num_steps_ref
    trajectory_ref(i).time = t_ref(i);
    trajectory_ref(i).state = state_ref(i, 1:9); % Extract r, v
    trajectory_ref(i).parameters = state_ref(i, 7:9); % Extractw_dmc
    trajectory_ref(i).STM = STM_ref(i, :, :);
    trajectory_ref(i).function_handle = @(t, x) sc_ODE(t, ...
        x, coeffs_ref, mu, l_max_ref, acc_func_ref, A_func_ref, acc_DMC, A_DMC);
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
             'o', 'DisplayName', sprintf('GS %d', station_idx), ...
             'MarkerSize', 8);
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
             'o', 'DisplayName', sprintf('GS %d', station_idx), ...
             'MarkerSize', 8);
    end
end
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Range-Rate Residual (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 8, 'Interpreter', 'latex');
grid on;
hold off;

%% Estimation - Iteration to Find Optimal Q

% Define initial covariance matrix P0
sigma_pos = 1;       % Std for position states (km)
sigma_vel = 1e-3;    % Std for velocity states (km/s)
sigma_acc = 1e-11;    % Std for uncompensated acc states (km/s^2)

% Save results with suffix
filename_suffix = 'EKF_tau';

P0 = diag([sigma_pos^2, sigma_pos^2, sigma_pos^2, ...
           sigma_vel^2, sigma_vel^2, sigma_vel^2, ...
           sigma_acc^2, sigma_acc^2, sigma_acc^2]);

% Range of sigma_Q_cont values for testing
sigma_Q_cont_values = logspace(-18, -5, 20);  % From 1e-18 to 1e-5
rms_pos_vals = zeros(length(sigma_Q_cont_values), 1);
rms_vel_vals = zeros(length(sigma_Q_cont_values), 1);
rms_residuals_vals_1 = zeros(length(sigma_Q_cont_values), 1);  % For the first residual (sigma_rho)
rms_residuals_vals_2 = zeros(length(sigma_Q_cont_values), 1);  % For the second residual (sigma_rho_dot)

% Store residuals for later use (we will store 3D residuals)
residuals_outside_3sigma = zeros(length(sigma_Q_cont_values), 1); % Count residuals outside 3σ

% Loop over different sigma_Q_cont values
for i = 1:length(sigma_Q_cont_values)
    sigma_Q_cont = sigma_Q_cont_values(i);

    % Initialize Process Noise Settings with the current sigma_Q_cont
    settings_PN = struct();
    settings_PN.threshold = 10;                        % Maximum dt value
    settings_PN.Q_cont = sigma_Q_cont^2 * eye(3);      % Continuous-time covariance matrix (in km/s^2)
    settings_PN.frame_type = 'ECI';                    % Frame type (either 'RIC' or 'ECI')
    settings_PN.method = 'DMC';                        % Method choice ('SNC' or 'DMC')
    settings_PN.B = B_DMC;                             % Time Coefficients DMC (in 1/s)

    % Run the Filter (LKF, or EKF)
    results = ekf(trajectory_ref, sorted_measurements, P0, settings_PN);

    % Compute RMS errors for the current sigma_Q_cont
    [rms_pos_3d, rms_vel_3d, rms_residuals] = compute_rms_errors(results, ...
        trajectory_true, sorted_measurements);

    % Store the RMS values
    rms_pos_vals(i) = rms_pos_3d;
    rms_vel_vals(i) = rms_vel_3d;
    rms_residuals_vals_1(i) = rms_residuals(1);  % Store the first residual (sigma_rho)
    rms_residuals_vals_2(i) = rms_residuals(2);  % Store the second residual (sigma_rho_dot)

    % Count how many residuals are outside 3σ bounds (directly from postfit_residuals)
    out_of_bounds_1 = sum(abs(results.postfit_residuals(1, :)) > 3 * sigma_rho);
    out_of_bounds_2 = sum(abs(results.postfit_residuals(2, :)) > 3 * sigma_rho_dot);

    % Total number of residuals outside 3σ bounds
    residuals_outside_3sigma(i) = out_of_bounds_1 + out_of_bounds_2;

    % Print progress for each iteration
    fprintf('Iteration %d/%d: sigma_Q_cont = %.2e, RMS Position = %.4f, RMS Velocity = %.4f, RMS Residual 1 = %.4f, RMS Residual 2 = %.4f\n', ...
        i, length(sigma_Q_cont_values), sigma_Q_cont, rms_pos_3d, rms_vel_3d, rms_residuals(1), rms_residuals(2));
    fprintf('\n');
end

%% Plot for Optimal Q

% Create the figure
gca3 = figure(3);
set(gcf, 'Position', [100, 100, 1400, 900]);  % Set figure size
subplot(3, 1, 1);
loglog(sigma_Q_cont_values, rms_pos_vals, '-o', 'Color', 'k', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 3);  % Blue for position error
xlabel('Process Noise Covariance ($\sigma_{Q_{cont}}$)', ...
    'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Position Error (km)', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
subplot(3, 1, 2);
loglog(sigma_Q_cont_values, rms_vel_vals, '-s', 'Color', 'k', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 3);  % Red for velocity error
xlabel('Process Noise Covariance ($\sigma_{Q_{cont}}$)', 'Interpreter', ...
    'latex', 'FontSize', 14);
ylabel('RMS Velocity Error (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
subplot(3, 1, 3);
loglog(sigma_Q_cont_values, rms_residuals_vals_1, '-o', 'Color', 'k','MarkerSize', 10, ...
    'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'LineWidth', 3);  % Green for residual 1
hold on;
loglog(sigma_Q_cont_values, rms_residuals_vals_2, '-x', 'Color', 'k', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'LineWidth', 3);  % Magenta for residual 2
xlabel('Process Noise Covariance ($\sigma_{Q_{cont}}$)', 'Interpreter', ...
    'latex', 'FontSize', 14);
ylabel('RMS Residuals', 'Interpreter', 'latex', 'FontSize', 14);
legend({'RMS Residual 1', 'RMS Residual 2'}, 'Interpreter', 'latex', ...
    'FontSize', 12, 'Location', 'best');
grid on;
exportgraphics(gca3, sprintf('Optimal_sigma_DMC_%s.pdf', filename_suffix), 'ContentType','image',...
    'Resolution', 1000);

%% Find the Optimal sigma_Q_cont

% Find the index of the optimal sigma_Q_cont (minimum number of residuals outside 3sig)
[~, optimal_index] = min(residuals_outside_3sigma);
optimal_sigma_Q_cont = sigma_Q_cont_values(optimal_index);

% Run the filter with the optimal sigma_Q_cont
settings_PN.Q_cont = optimal_sigma_Q_cont^2 * eye(3);
results = ekf(trajectory_ref, sorted_measurements, P0, settings_PN);

% Compute and display RMS errors for the optimal sigma
[rms_pos_3d_optimal, rms_vel_3d_optimal, rms_residuals_optimal] = compute_rms_errors(results, ...
    trajectory_true, sorted_measurements);

fprintf('\n');
disp(['Optimal sigma_Q_cont: ', num2str(optimal_sigma_Q_cont)]);
disp(['RMS Position (optimal): ', num2str(rms_pos_3d_optimal)]);
disp(['RMS Velocity (optimal): ', num2str(rms_vel_3d_optimal)]);
disp(['RMS Residuals (optimal 1): ', num2str(rms_residuals_optimal(1))]);
disp(['RMS Residuals (optimal 2): ', num2str(rms_residuals_optimal(2))]);

%% Plot filtering results

% Initialize
measurement_times = [sorted_measurements.time];
trajectory_times = [trajectory_true.time];
[~, traj_indices] = ismember(measurement_times, trajectory_times);

% Extract true states at exact measurement times
true_states_at_meas = zeros(9, length(measurement_times));
for i = 1:length(measurement_times)
    true_states_at_meas(:, i) = trajectory_true(traj_indices(i)).state(:);
end

% Compute state error (only at measurement times)
state_error = true_states_at_meas - results.state_corrected_hist;

% Compute 3-sig bounds from covariance matrix at measurement times
sigma_bounds_meas = zeros(size(state_error));
T_meas = length(measurement_times);
for t = 1:T_meas
    sigma_bounds_meas(:, t) = 3 * sqrt(diag(results.P_hist(:, :, t)));
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
    scatter(measurement_times_hours, 3 * sigma_bounds_meas(i, :), 5, ...
        colors(i, :), 'filled', 'DisplayName', '$\pm3\sigma$');
    scatter(measurement_times_hours, -3 * sigma_bounds_meas(i, :), 5, ...
        colors(i, :), 'filled', 'HandleVisibility', 'off');
    xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(y_labels{(i > 3) + 1}, 'Interpreter', 'latex', 'FontSize', 14);
    legend('show', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    hold off;
end
exportgraphics(gca4, sprintf('StateErrorsDMC_%s.pdf', filename_suffix), ...,
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
    scatter(measurement_times_hours, 3 * sigma_bounds_meas(i, :), 5, ...
        colors(i, :), 'filled', 'DisplayName', '$+3\sigma$');
    xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(y_labels{(i > 3) + 1}, 'Interpreter', 'latex', 'FontSize', 14);
    legend('show', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    hold off;
end
exportgraphics(gca5, sprintf('StateErrorsLogDMC_%s.pdf', filename_suffix), ...,
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
ylim([-10 * sigma_rho_dot, 10 * sigma_rho_dot]);
exportgraphics(gca6, sprintf('PostFitDMC_%s.pdf', filename_suffix), ...,
    'ContentType','image', 'Resolution', 1000);

%% Plot Estimated Uncompensated Acceleration 

% Plot
colors = lines(6);
gca14 = figure(14);
set(gcf, 'Position', [100, 100, 1400, 900]);  % Set figure size
subplot(3,1,1);
scatter(measurement_times_hours, results.state_corrected_hist(7, :), 10, ...
    colors(1, :), 'filled', 'MarkerEdgeColor', 'k');
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$w_x\cdot 10^{-9}$ (km/s$^2$)', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
subplot(3,1,2);
scatter(measurement_times_hours, results.state_corrected_hist(8, :), 10, ...
    colors(2, :), 'filled', 'MarkerEdgeColor', 'k');
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$w_y\cdot 10^{-9}$ (km/s$^2$)', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
subplot(3,1,3);
scatter(measurement_times_hours, results.state_corrected_hist(9, :), 10, ...
    colors(3, :), 'filled', 'MarkerEdgeColor', 'k');
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$w_z\cdot 10^{-9}$ (km/s$^2$)', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
exportgraphics(gca14, sprintf('DMC_estimates_%s.pdf', filename_suffix), ...,
    'ContentType','image', 'Resolution', 1000);

%% Difference in Dynamics

% Compute the accelerations for true, reference, and estimated states at the measurement times
acc_true = zeros(length(measurement_times), 3);  % True accelerations at measurement times
acc_ref = zeros(length(measurement_times), 3);   % Reference accelerations at measurement times
acc_estimated = zeros(length(measurement_times), 3);  % Estimated accelerations (from results)

for i = 1:length(measurement_times)
    % Extract the estimated state at the measurement time
    state_estimated_at_meas = [results.state_corrected_hist(1:6, i); 0; 0; 0; reshape(eye(9), [], 1)];

    % Calculate acceleration for true state at this time
    acc_true_total = sc_ODE(measurement_times(i), state_estimated_at_meas, ...
                            coeffs_true, mu, l_max_true, acc_func_true, A_func_true, acc_DMC, A_DMC);
    acc_true(i, :) = acc_true_total(4:6);  % Extract true acceleration components

    % Calculate acceleration for reference state at this time
    acc_ref_total = sc_ODE(measurement_times(i), state_estimated_at_meas, ...
                           coeffs_ref, mu, l_max_ref, acc_func_ref, A_func_ref, acc_DMC, A_DMC);
    acc_ref(i, :) = acc_ref_total(4:6);  % Extract reference acceleration components

    % Calculate acceleration for the estimated state at this time
    acc_estimated(i, :) = results.state_corrected_hist(7:9, i);  % Extract estimated acceleration components
end

% Compute the difference between the true and reference accelerations in x, y, z
acc_diff_x = acc_true(:, 1) - acc_ref(:, 1) - acc_estimated(:, 1);
acc_diff_y = acc_true(:, 2) - acc_ref(:, 2) - acc_estimated(:, 2);
acc_diff_z = acc_true(:, 3) - acc_ref(:, 3) - acc_estimated(:, 3);

% Compute 3-sig bounds from covariance matrix at measurement times
sigma_bounds = zeros(3, length(measurement_times));
T_meas = length(measurement_times);
for t = 1:T_meas
    sigma_bounds(:, t) = 3 * sqrt(diag(results.P_hist(7:9, 7:9, t)));  % Compute 3-sigma bounds for accelerations
end

% Plot acceleration errors (both positive and negative) in linear scale
gca24 = figure(24);
set(gcf, 'Position', [100, 100, 1400, 900]);

% Labels and y-axis titles for acceleration
labels = {'$\delta a_x$', '$\delta a_y$', '$\delta a_z$'};
y_labels = {'Acceleration Error (km/s$^2$)'};

% Colors for plotting
colors = lines(3);

% Plot for each acceleration component (x, y, z)
for i = 1:3
    subplot(3,1,i);

    % Select the correct acceleration difference based on the loop index
    if i == 1
        acc_diff = acc_diff_x;
    elseif i == 2
        acc_diff = acc_diff_y;
    else
        acc_diff = acc_diff_z;
    end

    % Plot the absolute value of the error for the corresponding component
    semilogy(measurement_times / 3600, abs(acc_diff), ...
         'o', 'MarkerSize', 5, 'MarkerFaceColor', colors(i, :), ...
         'MarkerEdgeColor', 'k', 'DisplayName', labels{i});
    hold on;

    % Plot the 3-sigma bounds for both positive and negative bounds
    scatter(measurement_times / 3600, 3 * sigma_bounds(i, :), 5, ...
            colors(i, :), 'filled', 'DisplayName', '$+3\sigma$');

    % Set labels, legend, and grid
    xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(y_labels{1}, 'Interpreter', 'latex', 'FontSize', 14);
    legend('show', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    hold off;
end
exportgraphics(gca24, sprintf('dmcErrors_%s.pdf', filename_suffix), ...,
    'ContentType','image', 'Resolution', 1000);


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

function [acceleration_func, A_func] = precomputeGaussMarkov(B)
    % Define the system's state-space model: dw/dt = f(w) = -B*w
    % This is the function representing the dynamics of the Gauss-Markov process
    acceleration_func = @(w) - B * w;

    % Convert Jacobian matrix A to a numeric function
    A_func = @(w) - B ;
end

function [acceleration, A] = sc_Dynamics(position_eci, w_eci, coeffs_vals, mu_val, ...
    ~, t, acc_SPH, A_SPH, acc_DMC, A_DMC)

    % Convert ECI to ECEF coordinates
    R_eci_to_ecef = DCM_eciToEcef(t);               % Rotation matrix ECI to ECEF
    position_ecef = R_eci_to_ecef * position_eci;   % Transform position to ECEF
 
    % Evaluate acceleration and parameter derivatives at the given inputs
    acceleration_SPH = acc_SPH(position_ecef(1), position_ecef(2), ...
                                     position_ecef(3), mu_val, coeffs_vals');
    acceleration_SPH(1:3) = R_eci_to_ecef' * acceleration_SPH(1:3);
    acceleration_DMC = acc_DMC(w_eci);

    % Sum acceleration contributions
    acceleration = [acceleration_SPH(1:3) + w_eci; acceleration_DMC(1:3)];
    
   if nargout > 1
    % Evaluate Jacobian of zonal harmonics (SPH) at the given position and parameters
    A_SPH_ecef = A_SPH(position_ecef(1), position_ecef(2), position_ecef(3), ...
                    0, 0, 0, mu_val, coeffs_vals');

    % Evaluate Jacobian of DMC in ECI
    A_DMC_eci = A_DMC(w_eci);

    % Transform zonal harmonics Jacobian back to ECI
    Ar_eci = R_eci_to_ecef' * A_SPH_ecef(4:6, 1:3) * R_eci_to_ecef; 

    % Construct the full state transition matrix A
    A = [zeros(3, 3), eye(3), zeros(3, 3);          % Velocity derivatives
         Ar_eci, zeros(3, 3), eye(3);               % Acceleration derivatives
         zeros(3, 6), A_DMC_eci];                   % DMC derivatives     
   end
end


function dstate = sc_ODE(t, state, coeffs, mu, l_max, acc_SPH, A_SPH, acc_DMC, A_DMC)
    % Extract from state
    position_eci = state(1:3);
    velocity_eci = state(4:6);
    w_eci = state (7:9);

    % Compute acceleration and A matrix
    [acceleration_eci, A_eci] = sc_Dynamics(position_eci, w_eci, coeffs, mu, ...
         l_max, t, acc_SPH, A_SPH, acc_DMC, A_DMC);

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
        d_measurements_dX = [d_rho_dR, d_rho_dV, zeros(1, 3);
                             d_rho_dot_dR, d_rho_dot_dV, zeros(1, 3)];

        % Ground station position Jacobian
        d_rho_dRs = -d_rho_dR;
        d_rho_dot_dRs = -d_rho_dot_dR;

        % Ground station state Jacobian
        d_measurements_dRs = [d_rho_dRs;
                             d_rho_dot_dRs];

        % Store function handle for measurement model
        measurement_func = @(x) measurement_function(x, gs_pos(:, i), gs_vel(:, i));

        % Store results in structure
        measurements_struct(i).time = times(i);
        measurements_struct(i).observed = [rho_obs; rho_dot_obs];
        measurements_struct(i).computed = [rho_comp; rho_dot_comp];
        measurements_struct(i).residual = [rho_res; rho_dot_res];
        measurements_struct(i).partials = struct('wrt_X', d_measurements_dX, ...
                                                 'wrt_Rs', d_measurements_dRs);
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
    H = [d_rho_dR, d_rho_dV, zeros(1, 3);
         d_rho_dot_dR, d_rho_dot_dV, zeros(1, 3)];
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
        stacked_partials_Rs = [];
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
            stacked_partials_Rs = [stacked_partials_Rs; ...
                all_measurements(j).partials.wrt_Rs];  %#ok<AGROW>
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
            'wrt_Rs', stacked_partials_Rs);
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

function results = lkf(trajectory_ref, sorted_measurements, P0, settings_PN)
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
        x_ref = trajectory_ref(traj_idx).state(1:n)';

        % Compute STM transition matrix
        % OSS: if ill-conditioned, either:
        % 1. use pinv(..., lowest possible tol) 2. do integration
        % by reset STM0 = eye(n) each time. STM_dt should be avoided!
        % STM_dt = STM_t / STM_tm;
        STM_dt = STM_t * pinv(STM_tm, 1e-16); 

        % Extract measurement data
        prefit_residual = sorted_measurements(t).residual;
        H_tilde = sorted_measurements(t).partials.wrt_X;
        R = diag(sorted_measurements(t).covariance);

        % Compute Process Noise
        if nargin > 3
            dt = 0;  % Default value for the first iteration
            if t > 1
                dt = measurement_times(t) - measurement_times(t-1); % Time step
            end
            Q_disc = calculate_Q_discrete(dt, x_ref, settings_PN);  % Store Q_disc
        else 
            Q_disc = zeros(6);
        end

        % Prediction Step
        dx_pred = STM_dt * dx;
        P_pred = STM_dt * P * STM_dt' + Q_disc;
        

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

function results = ekf(trajectory_ref, sorted_measurements, P0, settings_PN)
    % Extract state and parameter dimensions dynamically
    n = size(trajectory_ref(1).state, 2);                   % Number of state variables 
    % n_par = ...
    n_full = n;                                             % Number of total variables in f(x)

    % TODO: initialize with LKF!

    % Initialization
    T = length(sorted_measurements); % Number of measurements
    max_m = max(arrayfun(@(s) length(s.residual), sorted_measurements));

    state_corrected_hist = zeros(n, T);  % Corrected state estimates
    state_deviation_hist = zeros(n, T);  % Cumulative state deviations
    P_hist = zeros(n, n, T);
    postfit_residuals = NaN(max_m, T);  % Preallocate with NaN

    % Extract trajectory times for reference
    trajectory_times = [trajectory_ref.time];

    % Find indices of trajectory points that match measurement times
    measurement_times = [sorted_measurements.time];
    [~, traj_indices] = ismember(measurement_times, trajectory_times);

    % Initial state estimate (from first trajectory reference)
    x0 = trajectory_ref(traj_indices(1)).state(:);
    STM0 = eye(n_full);  % Initial STM
    x_aug = [x0; STM0(:)];  % Augmented state with STM
    P = P0;  % Initial covariance

    % ODE options for propagation
    ode_options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

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
            STM_t = STM_t(1:n, 1:n);

            % Compute Process Noise
            if nargin > 3
                dt = 0;  % Default value for the first iteration
                if t > 1
                    dt = measurement_times(t) - measurement_times(t-1); % Time step
                end
                Q_disc = calculate_Q_discrete(dt, x_propagated(end, 1:6), settings_PN);  % Store Q_disc
            else 
                Q_disc = zeros(6);
            end
    
            % Do actual Prediction Step
            x_pred = x_propagated(end, 1:n)';  
            P_pred = STM_t * P * STM_t' + Q_disc;
        else
            % No propagation needed, use current state
            x_pred = x_aug(1:n);  
            P_pred = P;
        end

        % Extract Measurement Model Function Handle and Covariance
        measurement_function = sorted_measurements(t).function_handle;
        R = diag(sorted_measurements(t).covariance);

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
        P_upd = (eye(n) - K * H_tilde) * P_pred * (eye(n) - K * H_tilde)' + K * R * K';

        % Store results
        state_corrected_hist(:, t) = x_upd;
        state_deviation_hist(:, t) = dx_upd; 
        P_hist(:, :, t) = P_upd;
        postfit_residuals(1:length(observed_measurement), t) = postfit_res;

        % Update for next iteration
        x_aug = [x_upd; STM0(:)];  % Update augmented state for next step
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

function [rms_pos_3d, rms_vel_3d, rms_residuals] = compute_rms_errors(results, trajectory_true, sorted_measurements)
    % Extract measurement times
    measurement_times = [sorted_measurements.time];
    trajectory_times = [trajectory_true.time];
    
    % Find indices of true trajectory corresponding to measurement times
    [~, traj_indices] = ismember(measurement_times, trajectory_times);
    
    % Extract true states at measurement times
    true_states_at_meas = zeros(9, length(measurement_times));
    for i = 1:length(measurement_times)
        true_states_at_meas(:, i) = trajectory_true(traj_indices(i)).state(:);
    end

    % Compute state errors (only at measurement times)
    state_errors = true_states_at_meas - results.state_corrected_hist; % (6 x T)
    postfit_res = results.postfit_residuals; % (m x T)

    % Compute 3D RMS Position and Velocity Errors
    rms_pos_3d = sqrt(mean(sum(state_errors(1:3, :).^2, 1)));
    rms_vel_3d = sqrt(mean(sum(state_errors(4:6, :).^2, 1)));

    % Compute RMS for post-fit residuals
    rms_residuals = sqrt(mean(postfit_res.^2, 2));
end
