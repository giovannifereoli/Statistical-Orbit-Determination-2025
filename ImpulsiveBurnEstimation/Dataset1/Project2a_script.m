%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         January 22, 2024
% CLASS:        ASEN 6080: StatOD
% INSTRUCTOR:   Prof. Jay W. McMahon
% ASSIGNMENT:   Project 2 - Part 2a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% Parse Data

% Define the filename
filename = 'Project2a_obs.txt';

% Read the data
data = readmatrix(filename);

% Extract time and measurement columns
timeVec     = data(:,1);
rangeDSS34  = data(:,2);
rangeDSS65  = data(:,3);
rangeDSS13  = data(:,4);
rrDSS34     = data(:,5);
rrDSS65     = data(:,6);
rrDSS13     = data(:,7);

% Initialize storage
time_all = [];
stationID_all = [];
range_all = [];
rangeRate_all = [];

% Go row by row (time is already sorted in file)
for i = 1:length(timeVec)
    t = timeVec(i);

    if ~isnan(rangeDSS34(i)) && ~isnan(rrDSS34(i))
        time_all      = [time_all; t];
        stationID_all = [stationID_all; 34];
        range_all     = [range_all; rangeDSS34(i)];
        rangeRate_all = [rangeRate_all; rrDSS34(i)];
    end
    if ~isnan(rangeDSS65(i)) && ~isnan(rrDSS65(i))
        time_all      = [time_all; t];
        stationID_all = [stationID_all; 65];
        range_all     = [range_all; rangeDSS65(i)];
        rangeRate_all = [rangeRate_all; rrDSS65(i)];
    end
    if ~isnan(rangeDSS13(i)) && ~isnan(rrDSS13(i))
        time_all      = [time_all; t];
        stationID_all = [stationID_all; 13];
        range_all     = [range_all; rangeDSS13(i)];
        rangeRate_all = [rangeRate_all; rrDSS13(i)];
    end
end

% Set the percentage of data to keep
percentage_to_use = 1;

% Total number of data points
N_total = length(time_all);

% Compute number of points to keep
N_keep = floor(percentage_to_use * N_total);

% Apply selection
trackingData.time       = time_all(1:N_keep);
trackingData.stationID  = stationID_all(1:N_keep);
trackingData.range      = range_all(1:N_keep);
trackingData.rangeRate  = rangeRate_all(1:N_keep);

%% Constants and Parameters

% Parameters
JD0         = 2456296.25;                     % Initial Julian Date
mu_E        = 398600.432896939;               % Earth grav. param [km^3/s^2]
mu_S        = 132712440017.987;               % Sun grav. param [km^3/s^2]
AMR         = 0.01 / 1e6;                     % Area-to-mass ratio [km^2/kg]
Cr          = 1.2;                            % SRP reflectivity coefficient
c = 299792458;                                % Speed of light [m/s]
Pphi = (1357 / c) * 1e3;                      % SRP pressure in N/km^2
AU = 149597870.7;                             % Astro Units in km

% Initial Conditions (pos/vel in ECI)
r0 = [-274096790.0; -92859240.0; -40199490.0];  % km
v0 = [32.67; -8.94; -3.88];                     % km/s
par0 = Cr; 

% Initial full state vector: [x; v; Cr; STM]
state0_true = [r0; v0; par0; reshape(eye(7), [], 1)];

%% Comapare Instructor Truth Data

% Initialize
load('Project2_Prob2_truth_traj_50days.mat'); 
t_span = Tt_50;
flag_STM = 1;
n = 7;
options = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-16);

% Repropagate Your Trajectory with STM
[~, state_true] = ode113(@(t, state) ...
    flyby_ODE(t, state, JD0, mu_E, mu_S, AMR, Pphi, AU, flag_STM), ...
    t_span, Xt_50(1,:), options);

% State Error Plotting
state_error = state_true(:,1:6) - Xt_50(:,1:6);  % position/velocity only
colors = turbo(6);  % One color per state component
component_labels = {'$x$', '$y$', '$z$', '$\dot{x}$', '$\dot{y}$', '$\dot{z}$'};
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.75, 0.6]);
time_days = t_span / 86400;
plot_handles = gobjects(6,1);
for i = 1:6
    plot_handles(i) = semilogy(time_days, abs(state_error(:,i)), ...
        'LineWidth', 2.0, 'Color', colors(i,:));
    hold on;
end
grid on; box on;
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('State Error (km, km/sec)', 'Interpreter', 'latex', 'FontSize', 16);
legend(plot_handles, component_labels, ...
    'Interpreter', 'latex', ...
    'FontSize', 14, ...
    'Location', 'northeastoutside');
set(gca, 'FontSize', 14);
exportgraphics(gca, 'state_error.pdf', 'ContentType','vector');


% STM Comparison
N = length(t_span);
STM_true = zeros(n, n, N);
STM_instr = zeros(n, n, N);
STM_diff = zeros(n, n, N);
for k = 1:N
    STM_true(:,:,k) = reshape(state_true(k, 8:end), n, n);        % Your propagated STM
    STM_instr(:,:,k) = reshape(Xt_50(k, 8:end), n, n);            % Instructor’s STM
    STM_diff(:,:,k) = STM_true(:,:,k) - STM_instr(:,:,k);         % Difference
end

% All STM Element Errors Over Time 
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.85, 0.85]);
time_days = t_span / 86400;
plot_handles = gobjects(n^2,1);
legend_entries = strings(n^2,1);
colors = turbo(n^2);
idx = 1;
for i = 1:n
    for j = 1:n
        element_error = squeeze(STM_diff(i,j,:));
        plot_handles(idx) = semilogy(time_days, abs(element_error), ...
            'LineWidth', 2.0, ...
            'Color', colors(idx,:));
        hold on;
        legend_entries(idx) = sprintf('$\\Phi_{%d%d}$', i, j);
        idx = idx + 1;
    end
end
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('STM Error $\Delta\Phi_{ij}$', 'Interpreter', 'latex', 'FontSize', 16);
grid on; box on;
legend(plot_handles, legend_entries, ...
    'Interpreter', 'latex', ...
    'FontSize', 15, ...
    'NumColumns', 10, ...
    'Location', 'southeast');
set(gca, 'FontSize', 14);
exportgraphics(gca, 'stm_error.pdf', 'ContentType','vector');

%% Simulate Reference Trajectories

% Time span (fixed duration, e.g., 50 days)        
t_span = trackingData.time;             % [s]
options = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-16);
flag_STM = 1;

% Integrate true trajectory
[t_ref, state_ref] = ode113(@(t, state) ...
    flyby_ODE(t, state, JD0, mu_E, mu_S, AMR, Pphi, AU, flag_STM), ...
    t_span, state0_true, options);


%% Plot Orbits in Astronomical Units (AU), with Sun and Earth

% Preallocate Earth (Sun-to-Earth) position array
JD_vec = JD0 + t_ref / 86400;
r_earth = zeros(length(t_ref), 3);

% Loop to compute Sun-to-Earth vector at each epoch
for k = 1:length(t_ref)
    [r_earth(k,:), ~] = Ephemeride(JD_vec(k), 3, mu_S);  % Earth wrt Sun
end

% Compute Sun-relative trajectory
r_true_sun = state_ref(:,1:3) + r_earth;  % Spacecraft wrt Sun
r_true_sun_AU = r_true_sun / AU;
r_earth_AU = r_earth / AU;

% Plot
figure;
plot3(r_true_sun_AU(:,1), r_true_sun_AU(:,2), r_true_sun_AU(:,3), 'b', 'LineWidth', 2); hold on;
plot3(r_earth_AU(:,1), r_earth_AU(:,2), r_earth_AU(:,3), 'g--', 'LineWidth', 1.5);  % Earth orbit
plot3(0, 0, 0, 'yo', 'MarkerFaceColor', 'k', 'MarkerSize', 8);  % Sun at origin
axis equal;
xlabel('$x$ (AU)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$ (AU)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('$z$ (AU)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Spacecraft Trajectory', 'Earth Orbit', 'Sun (Origin)', ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest');
grid on;
view(2);
exportgraphics(gca, 'orbit.pdf', 'ContentType','vector');

%% Extract STM and Fill Trajectory Structs

% STM reshape
STM_ref = reshape(state_ref(:, 8:end), [], 7, 7);

% Preallocate structs
trajectory_ref = struct('time', {}, 'state', {}, 'STM', {}, ...
    'parameters', {}, 'function_handle', {});
for i = 1:length(t_ref)
    trajectory_ref(i).time = t_ref(i);
    trajectory_ref(i).state = state_ref(i,1:6);
    trajectory_ref(i).parameters = state_ref(i,7);
    trajectory_ref(i).STM = STM_ref(i,:,:);
    trajectory_ref(i).function_handle = @(t, x) ...
        flyby_ODE(t, x, JD0, mu_E, mu_S, AMR, Pphi, AU, flag_STM);
end

%% Create Measurement Dictionary 

% Define Measurement Noise Covariance Matrix R
sigma_rho = 5 * 1e-3; % Standard deviation of range measurement (km)
sigma_rho_dot = 0.5 * 1e-6; % Standard deviation of range rate (km/s)
R = diag([sigma_rho^2, sigma_rho_dot^2]); % Covariance matrix

% DSN Station geodetic data (lat, lon in deg, alt in km)
stations = struct();
stations.DSS34 = struct('lat', -35.398333, 'lon', 148.981944,   'alt', 0.691750);     % Canberra
stations.DSS65 = struct('lat',  40.427222, 'lon', -355.749444,   'alt', 0.834539);    % Madrid (positive longitude for Part 3)
stations.DSS13 = struct('lat',  35.247164, 'lon', 243.205000,   'alt', 1.07114904);   % Goldstone
station_ecef0 = latlon_to_ecef(stations);

% Compute Measurements Using the Radiometric Function
measurements_struct = radiometric_measurement(trackingData, state_ref', ...
    station_ecef0, R);

%% Plot Residuals (first time)

% Extract Residuals Per Station for Plotting
stations = unique(trackingData.stationID); % Unique station IDs
num_stations = length(stations);

% Initialize cell arrays to store residuals and times per station
range_residuals = cell(num_stations, 1);
range_rate_residuals = cell(num_stations, 1);
time_per_station = cell(num_stations, 1);

% Extract residuals for each station
for s = 1:num_stations
    station_id = stations(s);

    % Find indices corresponding to this station
    station_indices = find(trackingData.stationID == station_id);

    % Store time and residuals
    time_per_station{s} = trackingData.time(station_indices); % Time in seconds
    range_residuals{s} = cellfun(@(x) x(1), {measurements_struct(station_indices).residual});
    range_rate_residuals{s} = cellfun(@(x) x(2), {measurements_struct(station_indices).residual});
end

% Plot Residuals 
gca123 = figure('Position', [100, 100, 1400, 900]);
subplot(2, 1, 1);
hold on;
for s = 1:num_stations
    if ~isempty(range_residuals{s}) && length(time_per_station{s}) == length(range_residuals{s})
        plot(time_per_station{s} / 3600, range_residuals{s}, 'o', ...
             'DisplayName', sprintf('DSS%d', stations(s)), 'MarkerSize', 8, 'LineWidth', 1.2);
    end
end
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Range Residual $\cdot 10^{4}$(km)', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
subplot(2, 1, 2);
hold on;
for s = 1:num_stations
    if ~isempty(range_rate_residuals{s}) && length(time_per_station{s}) == length(range_rate_residuals{s})
        plot(time_per_station{s} / 3600, range_rate_residuals{s}, 'o', ...
             'DisplayName', sprintf('DSS%d', stations(s)), 'MarkerSize', 8, 'LineWidth', 1.2);
    end
end
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Range-Rate Residual (km/s)', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(gca123, 'pre_fit_a.pdf', 'ContentType','vector');

%% Orbit Determination

% Initialize
sigma_pos = 100;       % Std for position states (km)
sigma_vel = 1e-1;      % Std for velocity states (km/s)
sigma_Cr = 0.1;        % Std for parameter states (-)
P0 = diag([sigma_pos^2, sigma_pos^2, sigma_pos^2, ...
           sigma_vel^2, sigma_vel^2, sigma_vel^2, sigma_Cr^2]);

% Run the filter
% filename_suffix = 'BatchSrif';
results = batch_srif_filter(trajectory_ref, measurements_struct, P0);

%% Iterate (2)

% Time span (fixed duration, e.g., 50 days)        
t_span = trackingData.time;             % [s]
options = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-16);
flag_STM = 1;
eyeI = eye(7);

% Integrate trajectory
[t_ref2, state_ref2] = ode113(@(t, state) ...
    flyby_ODE(t, state, JD0, mu_E, mu_S, AMR, Pphi, AU, flag_STM), ...
    t_span, [results.state_corrected_hist(:, 1); eyeI(:)], options);

% STM reshape
STM_ref2 = reshape(state_ref2(:, 8:end), [], 7, 7);

% Preallocate structs
trajectory_ref2 = struct('time', {}, 'state', {}, 'STM', {}, ...
    'parameters', {}, 'function_handle', {});
for i = 1:length(t_ref)
    trajectory_ref2(i).time = t_ref2(i);
    trajectory_ref2(i).state = state_ref2(i,1:6);
    trajectory_ref2(i).parameters = state_ref2(i,7);
    trajectory_ref2(i).STM = STM_ref2(i,:,:);
    trajectory_ref2(i).function_handle = @(t, x) ...
        flyby_ODE(t, x, JD0, mu_E, mu_S, AMR, Pphi, AU, flag_STM);
end

% Compute Measurements Using the Radiometric Function
measurements_struct2 = radiometric_measurement(trackingData, state_ref2', ...
    station_ecef0, R);

% Run the filter
% filename_suffix = 'BatchSrif2';
results = batch_srif_filter(trajectory_ref2, measurements_struct2, P0);

%% Iterate (3)

% Time span (fixed duration, e.g., 50 days)        
t_span = trackingData.time;             % [s]
options = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-16);
flag_STM = 1;
eyeI = eye(7);

% Integrate trajectory
[t_ref3, state_ref3] = ode113(@(t, state) ...
    flyby_ODE(t, state, JD0, mu_E, mu_S, AMR, Pphi, AU, flag_STM), ...
    t_span, [results.state_corrected_hist(:, 1); eyeI(:)], options);

% STM reshape
STM_ref3 = reshape(state_ref3(:, 8:end), [], 7, 7);

% Preallocate structs
trajectory_ref3 = struct('time', {}, 'state', {}, 'STM', {}, ...
    'parameters', {}, 'function_handle', {});
for i = 1:length(t_ref)
    trajectory_ref3(i).time = t_ref3(i);
    trajectory_ref3(i).state = state_ref3(i,1:6);
    trajectory_ref3(i).parameters = state_ref3(i,7);
    trajectory_ref3(i).STM = STM_ref3(i,:,:);
    trajectory_ref3(i).function_handle = @(t, x) ...
        flyby_ODE(t, x, JD0, mu_E, mu_S, AMR, Pphi, AU, flag_STM);
end

% Compute Measurements Using the Radiometric Function
measurements_struct3 = radiometric_measurement(trackingData, state_ref3', ...
    station_ecef0, R);

% Run the filter
% filename_suffix = 'BatchSrif3';
results = batch_srif_filter(trajectory_ref3, measurements_struct3, P0);

%% Plot filtering results

% Time vector
measurement_times_hours = [measurements_struct.time] / 3600;

% Compute ±3σ bounds
n_states = size(results.P_hist, 1);
n_times = length(measurement_times_hours);
sigma_bounds = zeros(n_states, n_times);
for t = 1:n_times
    sigma_bounds(:, t) = 3 * sqrt(diag(results.P_hist(:,:,t)));
end

% Extract state history for overlay (for w and dv)
state_hist = results.state_corrected_hist;

% Labels for 13 states
labels = {'$x\cdot 10^{-3}$', '$y\cdot 10^{-3}$', '$z\cdot 10^{-3}$', ...
          '$\dot{x}\cdot 10^{-10}$', '$\dot{y}\cdot 10^{-7}$', '$\dot{z}\cdot 10^{-7}$', ...
          '$C_R$'};

% Group axis labels
y_labels = {'Position Error (km)', ...
            'Velocity Error (km/s)', ...
            'SRP Coeff. Error'};

% Colors
cmap = lines(n_states);

% Plot
gca234 = figure('Position', [100, 100, 1600, 1000]);
for i = 1:n_states
    subplot(3, 3, i);

    % Plot ±3σ bounds
    semilogy(measurement_times_hours, sigma_bounds(i, :), ...
             '.', 'Color', cmap(i,:), 'MarkerSize', 10, ...
             'DisplayName', ['$3\sigma$: ', labels{i}]);
    hold on;

    xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);

    % Y-axis label selection
    if i <= 3
        ylabel(y_labels{1}, 'Interpreter', 'latex', 'FontSize', 14);
    elseif i <= 6
        ylabel(y_labels{2}, 'Interpreter', 'latex', 'FontSize', 14);
    elseif i == 7
        ylabel(y_labels{3}, 'Interpreter', 'latex', 'FontSize', 14);
    end

    title(labels{i}, 'Interpreter', 'latex', 'FontSize', 16);
    grid on;
    legend('Interpreter', 'latex', 'FontSize', 10, 'Location', 'northeast');
end
exportgraphics(gca234, 'covariances_1a.pdf', 'ContentType','vector');

% Prep Info
range_res     = results.postfit_residuals(1, :);
rangerate_res = results.postfit_residuals(2, :);

range_mean = mean(range_res);
range_std  = std(range_res);

rangerate_mean = mean(rangerate_res);
rangerate_std  = std(rangerate_res);

% Post-fits Residuals Plot
gca456 = figure('Position', [100, 100, 1400, 800]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Range Plot
ylims_range = [-1, 1] * 10 * sigma_rho;
% Left: Range residuals over time
nexttile(1);
scatter(measurement_times_hours, range_res, 10, 'b', 'filled');
hold on;
fill([measurement_times_hours, fliplr(measurement_times_hours)], ...
     [3*sigma_rho*ones(size(measurement_times_hours)), ...
      fliplr(-3*sigma_rho*ones(size(measurement_times_hours)))], ...
     'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Range Residual (km)', 'Interpreter', 'latex', 'FontSize', 14);
title('Range Residuals vs. Time', 'Interpreter', 'latex', 'FontSize', 16);
ylim(ylims_range);
grid on;
% Right: Histogram (horizontal, aligned)
nexttile(2);
histogram(range_res, 50, 'Orientation', 'horizontal', ...
    'FaceColor', [0.2 0.2 0.8], 'FaceAlpha', 0.7);
hold on;
yline(range_mean, 'k-', 'LineWidth', 1.5);
yline(range_mean + 3*range_std, 'r--', 'LineWidth', 1.2);
yline(range_mean - 3*range_std, 'r--', 'LineWidth', 1.2);
xlabel('Count', 'Interpreter', 'latex', 'FontSize', 14);
title('Range Residual Histogram', 'Interpreter', 'latex', 'FontSize', 16);
ylim(ylims_range); % MATCH left plot
set(gca, 'YTickLabel', []); % Hide duplicate label
grid on;
legend({...
    sprintf('Mean = %.2e', range_mean), ...
    sprintf('Std = %.2e', range_std)}, ...
    'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast');

% Range-Rate Plots
ylims_rangerate = [-1, 1] * 10 * sigma_rho_dot;
% Left: Range-rate residuals over time
nexttile(3);
scatter(measurement_times_hours, rangerate_res, 10, 'r', 'filled');
hold on;
fill([measurement_times_hours, fliplr(measurement_times_hours)], ...
     [3*sigma_rho_dot*ones(size(measurement_times_hours)), ...
      fliplr(-3*sigma_rho_dot*ones(size(measurement_times_hours)))], ...
     'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Range-Rate Residual $\cdot 10^{-6}$(km/s)', 'Interpreter', 'latex', 'FontSize', 14);
title('Range-Rate Residuals vs. Time', 'Interpreter', 'latex', 'FontSize', 16);
ylim(ylims_rangerate);
grid on;
% Right: Histogram
nexttile(4);
histogram(rangerate_res, 50, 'Orientation', 'horizontal', ...
    'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.7);
hold on;
yline(rangerate_mean, 'k-', 'LineWidth', 1.5);
yline(rangerate_mean + 3*rangerate_std, 'r--', 'LineWidth', 1.2);
yline(rangerate_mean - 3*rangerate_std, 'r--', 'LineWidth', 1.2);
xlabel('Count', 'Interpreter', 'latex', 'FontSize', 14);
title('Range-Rate Residual Histogram', 'Interpreter', 'latex', 'FontSize', 16);
ylim(ylims_rangerate); % MATCH left plot
set(gca, 'YTickLabel', []); % Hide duplicate label
grid on;
legend({...
    sprintf('Mean = %.2e', rangerate_mean), ...
    sprintf('Std = %.2e', rangerate_std)}, ...
    'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast');
exportgraphics(gca456, 'postfits_a.pdf', 'ContentType','vector');

%% B-Plane Mapping and Plots

% Constants
R_SOI    = 924000;               % Earth's SOI [km]
r_target = 3 * R_SOI;            % Target radius for event

% ODE options (1)
opts_event = odeset( ...
    'Events',    @(t,x) reach_r_event(t, x, r_target), ...
    'RelTol',    1e-12, ...
    'AbsTol',    1e-14);

% Initial time and state for DCO propagation
t0_DCO    = measurements_struct(end).time;
t1_DCO    = 10 * t0_DCO;
x0_DCO    = [results.state_corrected_hist(:, end); reshape(eye(7), [], 1) ];

% Propagate until reaching 3 * RSOI
[t_rSOI, state_rSOI] = ode113( ...
    @(t, x) flyby_ODE(t, x, JD0, mu_E, mu_S, AMR, Pphi, AU, 1), ...
    [t0_DCO, t1_DCO], ...
    x0_DCO, ...
    opts_event);

% Extract position & velocity at event
r1       = state_rSOI(end, 1:3).';
v1       = state_rSOI(end, 4:6).';
v1_norm  = norm(v1);

% Build B‑plane frame
S_hat    = v1 / v1_norm;
N_hat    = [0; 0; 1];
T_hat    = cross(S_hat, N_hat);
T_hat    = T_hat / norm(T_hat);
R_hat    = cross(S_hat, T_hat);
DCM_ECI_to_B = [T_hat, R_hat, S_hat]';

% ODE Options (2)
opts_prop  = odeset( ...
    'Events', @(t, x) bplane_cross_event(t, x, S_hat), ...
    'RelTol',    1e-12, ...
    'AbsTol',    1e-14);

% Compute hyperbolic orbit parameters
r1_norm = norm(r1);
h_vec = cross(r1, v1);
h = norm(h_vec);
e_vec = (1 / mu_E) * ((v1_norm^2 - mu_E / r1_norm) * r1 - dot(r1, v1) * v1);
e = norm(e_vec);
a = -mu_E / (v1_norm^2 - 2 * mu_E / r1_norm);
cos_th = dot(r1, e_vec) / (r1_norm * e);
theta = acos(cos_th);
f = acosh( 1 + (v1_norm^2) / mu_E * (a * (1 - e^2)) / (1 + e * cos(theta)));
LTOF = (mu_E / v1_norm^3) * (sinh(f) - f);

% Propagate from 3*RSOI to B‑plane crossing
t_start_B = t_rSOI(end);
t_end_B   = t_start_B + 10 * LTOF;
x0_B      = state_rSOI(end, :)';
[t_bplane, state_bplane] = ode113( ...
    @(t, x) flyby_ODE(t, x, JD0, mu_E, mu_S, AMR, Pphi, AU, 1), ...
    [t_start_B, t_end_B], ...
    x0_B, ...
    opts_prop);

% Extract final position and STM
r_bplane          = state_bplane(end, 1:3)';
Phi_DCO_to_bplane = reshape(state_bplane(end, 8:end), 7, 7);

% Map the covariance into the B‑plane frame
P_DCO    = results.P_hist(:, :, end);
P_bplane = Phi_DCO_to_bplane * P_DCO * Phi_DCO_to_bplane';
P_r      = P_bplane(1:3, 1:3);
P_r_B    = DCM_ECI_to_B * P_r * DCM_ECI_to_B';

% Compute BdotT, BdotR, BdotS
B_eci      = r1 - dot(r1, v1) * v1;
r_bplane_B = DCM_ECI_to_B * B_eci;
BdotT      = r_bplane_B(1);
BdotR      = r_bplane_B(2);
BdotS      = r_bplane_B(3);

% Save results
totalDays = 200;  
spanDays  = percentage_to_use * totalDays; 
Rstruct = struct( ...
    'spanDays', spanDays, ...
    'P_r_B',    P_r_B, ...
    'BdotT',    BdotT, ...
    'BdotR',    BdotR, ...
    'BdotS',    BdotS ...
);
fname = sprintf('bplane_%ddays.mat', spanDays);
save(fname, 'Rstruct');


% Plot
figure('Position', [100, 100, 1400, 900]); hold on;
plot(BdotT, BdotR, 'x', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', [0 0 1]);
error_ellipse(P_r_B(1:2,1:2), [BdotT; BdotR], 'C', [0 0 1]);
BdotR_true  = 14970.824; BdotT_true  = 9796.737;
% plot( BdotT_true, BdotR_true, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'k' );
% text( BdotT_true + 200, BdotR_true + 200, 'Target', 'FontSize', 12 );
xlabel( '$B \cdot \hat{T}$ (km)', 'Interpreter', 'latex', 'FontSize', 14 );
ylabel( '$B \cdot \hat{R}$ (km)', 'Interpreter', 'latex', 'FontSize', 14 );
grid on; box on;
legend( 'DCO Mean', 'DCO Covariance', 'Location', 'best' );
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xtickformat( '%.0f' );
ytickformat( '%.0f' );
set(gca,'YDir','reverse');


% 8) DEBUG PRINTS
fprintf( '[DEBUG] t_DCO (days): %.6f\n', t_rSOI(1)/86400 );
fprintf( '[DEBUG] r - 3 * R_SOI (should be small): %.14f km\n', ...
    norm(state_rSOI(end,1:3)) - r_target );
fprintf( '[DEBUG] t_3rsoi (days): %.6f\n', t_rSOI(end)/86400 );
fprintf( '[DEBUG] BdotT = %.3f (true: %.3f) -> DELTA = %.3f\n', ...
    BdotT, BdotT_true, BdotT - BdotT_true );
fprintf( '[DEBUG] BdotR = %.3f (true: %.3f) -> DELTA = %.3f\n', ...
    BdotR, BdotR_true, BdotR - BdotR_true );
fprintf( '[DEBUG] t_bplane (days): %.6f\n', t_bplane(end)/86400 );
fprintf( '[DEBUG] BdotS (should be small): %.14f km\n',dot(r_bplane, S_hat));

%% Plot B‑plane Estimates for 50, 100, 150 & 200 day Runs 

% Initialize
spans = [50, 100, 150, 200];
cols  = lines(numel(spans));
figure('Position', [100, 100, 1400, 900]); hold on;
legendEntries = cell(1, 2*numel(spans) + 1);
for i = 1:numel(spans)
    % Load result struct
    fname = sprintf('bplane_%ddays.mat', spans(i));
    S     = load(fname);
    fn    = fieldnames(S);
    R     = S.(fn{1});

    % Extract mean and covariance
    B_T   = R.BdotT;
    B_R   = R.BdotR;
    P_TR  = R.P_r_B;

    % Plot covariance ellipse
    error_ellipse(3 * P_TR(1:2, 1:2), [B_T; B_R], 'C', cols(i,:), 'style', '-');
    legendEntries{2*i-1} = sprintf('DCO %d-day 3-sig Error Ellipse', spans(i));

    % Plot mean cross
    plot(B_T, B_R, 'x', ...
         'Color',    cols(i,:), ...
         'LineWidth', 2, ...
         'MarkerSize',12);
    legendEntries{2*i} = sprintf('DCO %d-day Target', spans(i));
end
% plot(9796.737, 14970.824, 'kp', 'MarkerFaceColor','k','MarkerSize',8);
legendEntries{end} = 'TCM Design';
xlabel('$B\!\cdot\!\hat T$ (km)', 'Interpreter','latex','FontSize',14);
ylabel('$B\!\cdot\!\hat R$ (km)', 'Interpreter','latex','FontSize',14);
grid on; axis equal;
legend(legendEntries, 'Location','northwest', 'FontSize',16, ...
       'Box','on', 'Interpreter','latex')
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
set(gca,'YDir','reverse');
exportgraphics(gca, 'bplanes_comparisons.pdf', 'ContentType','vector');

%% Compute Reduced Chi-Squared

% Flatten postfit residuals
postfit_residuals_flat = results.postfit_residuals(:);

% Build weight vector
n_meas = length(postfit_residuals_flat);
weights = repmat([sigma_rho; sigma_rho_dot], n_meas/2, 1);

% Normalize residuals
normalized_residuals = postfit_residuals_flat ./ weights;

% Number of estimated parameters
n_estimated = 7; % [x,y,z, vx,vy,vz, Cr]

% Degrees of Freedom (DOF)
nu = n_meas - n_estimated;

% Compute raw chi-squared
chi2_raw = sum(normalized_residuals.^2);

% Compute reduced chi-squared
chi2_reduced = chi2_raw / nu;

% Display
fprintf('\n========= Chi-Squared Evaluation =========\n');
fprintf('Raw Chi-Squared: %.6f\n', chi2_raw);
fprintf('Degrees of Freedom: %d\n', nu);
fprintf('Reduced Chi-Squared: %.6f\n', chi2_reduced);
fprintf('Ideal reduced value should be close to 1.0\n');
fprintf('===========================================\n');


%% Helper Functions Dynamics

function [value, isterminal, direction] = bplane_cross_event(~, state, S_hat)
    r = state(1:3);               % position in ECI
    value = dot(r, S_hat);        % event when r perp to S_hat
    isterminal = 1;               % stop when crossing plane
    direction = 1;                % trigger when going from negative to positive
end

function [value, isterminal, direction] = reach_r_event(~, state, r_target)
    r_vec = state(1:3);              % Extract position vector
    r_norm = norm(r_vec);            % Compute distance from origin
    value = r_norm - r_target;       % Zero crossing at r = r_target
    isterminal = 1;                  % Stop integration when event fires
    direction = -1;                   % Detect only when r is increasing
end

function error_ellipse(C, mu, varargin)
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'style', '-', @ischar);
    addParameter(p, 'C', [0 0 1]);  % default color
    parse(p, varargin{:});
    style = p.Results.style;
    color = p.Results.C;

    % Number of points to draw the ellipse
    theta = linspace(0, 2*pi, 100);

    % Eigen-decomposition
    [V, D] = eig(C);

    % Generate unit circle points
    circle = [cos(theta); sin(theta)];

    % Transform circle to ellipse
    ellipse = V * sqrt(D) * circle;

    % Translate ellipse to mean
    ellipse = ellipse + mu;

    % Plot
    fill(ellipse(1,:), ellipse(2,:), color, ...    % pale blue
      'EdgeColor', color, ...
      'LineWidth', 2, ...
      'FaceAlpha', 0.3);
end

function gs_ecef = latlon_to_ecef(stations)
    % Constants
    Re = 6378.1363;  % Earth radius [km]

    % Init output map
    gs_ecef = containers.Map();

    % Station names
    station_names = fieldnames(stations);

    for i = 1:length(station_names)
        name = station_names{i};

        % Convert lat/lon/alt to ECEF
        lat = deg2rad(stations.(name).lat);
        lon = deg2rad(stations.(name).lon);
        alt = stations.(name).alt;

        % ECEF position
        r_ecef = (Re + alt) * [cos(lat) * cos(lon);
                               cos(lat) * sin(lon);
                               sin(lat)];

        % Store in output map
        gs_ecef(name) = r_ecef;
    end
end


function R = DCM_eciToEcef(t)
     % Earth's rotation rate 
    omega_earth =7.29211585275553e-5; % Earth's rotation rate [rad/s]
    
    % Rotation angle
    theta = omega_earth * t;
    
    % Rotation matrix for ECI to ECEF
    R = [cos(theta), sin(theta), 0;
        - sin(theta), cos(theta), 0;
         0,          0,          1];
end

function dstate = flyby_ODE(t, state, JD0, mu_E, mu_S, AM_ratio, Pphi, AU, flag_STM)
    % Determine state size
    n = 7;  % 3 position + 3 velocity + 1 parameter (Cr)
    if flag_STM
        dstate = zeros(n + n^2,1);
    else
        dstate = zeros(n,1);
    end

    % Extract states
    r     = state(1:3);    % Position [km]
    v     = state(4:6);    % Velocity [km/s]
    Cr    = state(7);      % Reflectivity coefficient [-]

    % Sun position in ECI frame at current time
    JD = JD0 + t / 86400;
    [r_earth_to_sun, ~] = Ephemeride(JD, 3, mu_S);  % Earth wrt Sun
    r_sun = -r_earth_to_sun;                        % Sun wrt Earth

    % Compute relative vectors
    r_norm     = norm(r);
    rS_norm    = norm(r_sun);
    r_rel      = r_sun - r;
    r_rel_norm = norm(r_rel);
    r_rel_dir  = r_rel / r_rel_norm;

    % Accelerations from masses
    a_earth = -mu_E * r / r_norm^3;
    a_sun   = mu_S * (r_rel / r_rel_norm^3 - r_sun / rS_norm^3);

    % SRP acceleration
    a_srp = - AM_ratio * Cr * Pphi * (AU^2 / r_rel_norm^2) * r_rel_dir;  

    % Total acceleration
    a_total = a_earth + a_sun + a_srp;

    % Populate derivative vector
    dstate(1:3) = v;
    dstate(4:6) = a_total;
    dstate(7)   = 0;  % Cr constant

    % STM propagation (if enabled)
    if flag_STM == 1
        % Reshape STM from state
        STM = reshape(state(n+1:end), n, n);
        A = zeros(n);  

        % d r_dot / dv 
        A(1:3,4:6) = eye(3);

        % d a_earth / dr
        dadr_earth = mu_E * (3 * (r * r.') / r_norm^5 - eye(3) / r_norm^3);

        %  d a_sun/ dr
        dadr_sun = mu_S * (3 * (r_rel * r_rel.') / r_rel_norm^5 - eye(3) / r_rel_norm^3);

        % d a_srp / dr
        scale = AM_ratio * Cr * Pphi * AU^2;
        dadr_srp = - scale * (3 * (r_rel * r_rel.') / r_rel_norm^5 - eye(3) / r_rel_norm^3);

        % d a_total / dr, sum of all
        A(4:6,1:3) = dadr_earth + dadr_sun + dadr_srp;

        % d a_total / d C_r 
        A(4:6,7) = a_srp / Cr;

        % STM derivative
        Phi_dot = A * STM;
        dstate(n+1:end) = Phi_dot(:);
    end
end

function [gs_positions_eci, gs_velocities_eci] = compute_gs_eci(gs_location_ecef, times)
    % Constants
    omega_E = 7.29211585275553e-5;  % Earth rotation rate [rad/s]

    % Pre-allocate outputs
    N = length(times);
    gs_positions_eci = zeros(3, N);
    gs_velocities_eci = zeros(3, N);

    % Loop through times to compute ECI state
    for i = 1:N
        t = times(i);
        % Earth's rotation angle at time t
        R_eci2ecef = DCM_eciToEcef(t);

        % Ground station position in ECI
        gs_positions_eci(:, i) = R_eci2ecef' * gs_location_ecef;

        % Ground station velocity in ECI
        gs_velocities_eci(:, i) = omega_E * cross([0; 0; 1], ...
            gs_positions_eci(:, i));
    end
end

function measurements_struct = radiometric_measurement(trackingData, sc_ref, ...
    station_ecef0, R)
    
    % Compute reference spacecraft state in ECI
    sc_ref_pos = sc_ref(1:3, :);
    sc_ref_vel = sc_ref(4:6, :);

    % Preallocate measurements struct
    N = length(trackingData.time);
    measurements_struct = struct('time', cell(1, N), 'stationID', cell(1, N), ...
                                 'observed', cell(1, N), 'computed', cell(1, N), ...
                                 'residual', cell(1, N), 'partials', cell(1, N), ...
                                 'covariance', cell(1, N), 'measurement_function', cell(1, N));

    % Loop through each time step
    for i = 1:N
        % Identify the ground station ID for the current time
        station_id = trackingData.stationID(i);

        % Extract precomputed GS ECI position and velocity
        [gs_pos, gs_vel] = compute_gs_eci(station_ecef0("DSS" + string(station_id)), ...
            trackingData.time(i));

        % Reference relative position and velocity
        dR_ref = sc_ref_pos(:, i) - gs_pos;
        dV_ref = sc_ref_vel(:, i) - gs_vel;

        % Computed (on reference) range and range-rate
        rho_comp = norm(dR_ref);
        rho_dot_comp = dot(dR_ref, dV_ref) / rho_comp;

        % Observed (true) range and range-rate from trackingData
        rho_obs = trackingData.range(i);
        rho_dot_obs = trackingData.rangeRate(i);

        % Compute residuals
        rho_res = rho_obs - rho_comp;
        rho_dot_res = rho_dot_obs - rho_dot_comp;

        % Compute partial derivatives (on reference)
        d_rho_dR = dR_ref' / rho_comp; 
        d_rho_dV = zeros(1, 3);
        d_rho_dot_dR = (rho_comp * dV_ref' - rho_dot_comp * dR_ref') / rho_comp^2;
        d_rho_dot_dV = dR_ref' / rho_comp;

        % Spacecraft state Jacobian in ECI
        d_measurements_dX = [d_rho_dR, d_rho_dV, 0;
                             d_rho_dot_dR, d_rho_dot_dV, 0];

        % Store function handle for measurement model
        measurement_func = @(x) measurement_function(x, gs_pos, gs_vel);

        % Store results in structure
        measurements_struct(i).time = trackingData.time(i);
        measurements_struct(i).stationID = station_id;
        measurements_struct(i).observed = [rho_obs; rho_dot_obs];
        measurements_struct(i).computed = [rho_comp; rho_dot_comp];
        measurements_struct(i).residual = [rho_res; rho_dot_res];
        measurements_struct(i).partials = struct('wrt_X', d_measurements_dX);
        measurements_struct(i).covariance = R;
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

    % Full Jacobian Spacecraft
    H = [d_rho_dR, d_rho_dV, 0;
         d_rho_dot_dR, d_rho_dot_dV, 0];
end

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