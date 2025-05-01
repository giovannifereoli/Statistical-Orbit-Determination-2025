%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         January 22, 2024
% CLASS:        ASEN 6080: StatOD
% INSTRUCTOR:   Prof. Jay W. McMahon
% ASSIGNMENT:   Project 2 - Part 2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;


%% Parse Data

% Define the filename
filename = 'Project2b_obs.txt';

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
percentage_to_use = 0.725;

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
Cr          = 1.0;                            % SRP reflectivity coefficient
c = 299792458;                                % Speed of light [m/s]
Pphi = (1357 / c) * 1e3;                      % SRP pressure in N/km^2
AU = 149597870.7;                             % Astro Units in km
B = (60 * 24 * 3600)^-1 * eye(3);

% Initial Conditions (pos/vel in ECI)
r0 = [-274096770.76544; -92859266.44990661; -40199493.6677441];                    % km
v0 = [32.6704564599943; -8.93838913761049; -3.87881914050316];                     % km/s
[pos_sun, vel_sun] = planetEphemeris(JD0, 'Earth', 'Sun');
par0 = [Cr; pos_sun'; vel_sun']; 

% Initial full state vector: [x; v; Cr; STM]
state0_true = [r0; v0; par0; reshape(eye(13), [], 1)];

%% Simulate Reference Trajectories

% Time span (fixed duration, e.g., 50 days)        
t_span = trackingData.time;             % [s]
options = odeset('RelTol', 2.22045e-10, 'AbsTol', 2.22045e-12);
flag_STM = 1;

% Integrate true trajectory
[t_ref, state_ref] = ode113(@(t, state) ...
    flyby_ODE(t, state, JD0, B, mu_E, mu_S, AMR, Pphi, AU, flag_STM), ...
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

%% Extract STM and Fill Trajectory Structs

% STM reshape
STM_ref = reshape(state_ref(:, 14:end), [], 13, 13);

% Preallocate structs
trajectory_ref = struct('time', {}, 'state', {}, 'STM', {}, ...
    'parameters', {}, 'function_handle', {});
for i = 1:length(t_ref)
    trajectory_ref(i).time = t_ref(i);
    trajectory_ref(i).state = state_ref(i,1:6);
    trajectory_ref(i).parameters = state_ref(i,7:13);
    trajectory_ref(i).STM = STM_ref(i,:,:);
    trajectory_ref(i).function_handle = @(t, x) ...
        flyby_ODE(t, x, JD0, B, mu_E, mu_S, AMR, Pphi, AU, flag_STM);
end

%% Create Measurement Dictionary 

% Define Measurement Noise Covariance Matrix R
sigma_rho = 5 * 1e-3; % Standard deviation of range measurement (km)
sigma_rho_dot = 0.5 * 1e-6; % Standard deviation of range rate (km/s)
% sigma_rho = 10 * 1e-3; % Standard deviation of range measurement (km), 5
% sigma_rho_dot = 1 * 1e-6; % Standard deviation of range rate (km/s), 0.5
R = diag([sigma_rho^2, sigma_rho_dot^2]); % Covariance matrix

% DSN Station geodetic data (lat, lon in deg, alt in km)
stations = struct();
stations.DSS34 = struct('lat', -35.398333, 'lon', 148.981944,   'alt', 0.691750);     % Canberra
stations.DSS65 = struct('lat',  40.427222, 'lon', 355.749444,   'alt', 0.834539);     % Madrid (positive longitude for Part 3)
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
ylabel('Range Residual (km)', 'Interpreter', 'latex', 'FontSize', 16);
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

%% Orbit Determination

% Initialize
sigma_pos = 100;       % Std for position states (km)
sigma_vel = 1e-1;      % Std for velocity states (km/s)
sigma_Cr = 0.1;        % Std for parameter states (-)
sigma_sun_pos = 1e6;
sigma_vel_pos = 1e2;
P0 = diag([sigma_pos^2, sigma_pos^2, sigma_pos^2, ...
           sigma_vel^2, sigma_vel^2, sigma_vel^2, sigma_Cr^2, ...
           sigma_sun_pos^2, sigma_sun_pos^2, sigma_sun_pos^2, ...
           sigma_vel_pos^2, sigma_vel_pos^2, sigma_vel_pos^2]);

% Initialize Process Noise Settings with the current sigma_Q_cont
settings_PN = struct();
settings_PN.threshold = 10;                        % Maximum dt value
settings_PN.frame_type = 'ECI';                    % Frame type (either 'RIC' or 'ECI')
settings_PN.method = 'SNC';                        % Method choice ('SNC' or 'DMC')
settings_PN.B = B;
optimal_sigma_Q_cont = 1e-8;  % Continuous-time covariance matrix (in km/s^2)
settings_PN.Q_cont = optimal_sigma_Q_cont^2 * eye(3);

% Run the filter
% filename_suffix = 'BatchSrif';
results = lkf(trajectory_ref, measurements_struct, P0, settings_PN);

%% Iterate (2)

% Time span (fixed duration, e.g., 50 days)        
t_span = trackingData.time;             % [s]
options = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-16);
flag_STM = 1;
eyeI = eye(13);

% Integrate trajectory
[t_ref2, state_ref2] = ode113(@(t, state) ...
    flyby_ODE(t, state, JD0, B, mu_E, mu_S, AMR, Pphi, AU, flag_STM), ...
    t_span, [results.state_corrected_hist(:, 1); eyeI(:)], options);

% STM reshape
STM_ref2 = reshape(state_ref(:, 14:end), [], 13, 13);

% Preallocate structs
trajectory_ref2 = struct('time', {}, 'state', {}, 'STM', {}, ...
    'parameters', {}, 'function_handle', {});
for i = 1:length(t_ref)
    trajectory_ref2(i).time = t_ref2(i);
    trajectory_ref2(i).state = state_ref2(i,1:6);
    trajectory_ref2(i).parameters = state_ref2(i,7:13);
    trajectory_ref2(i).STM = STM_ref2(i,:,:);
    trajectory_ref2(i).function_handle = @(t, x) ...
        flyby_ODE(t, x, JD0, B, mu_E, mu_S, AMR, Pphi, AU, flag_STM);
end

% Compute Measurements Using the Radiometric Function
measurements_struct2 = radiometric_measurement(trackingData, state_ref2', ...
    station_ecef0, R);

% Run the filter
% filename_suffix = 'BatchSrif2';
results = lkf(trajectory_ref2, measurements_struct2, P0, settings_PN);
results = rts_smoother(results, trajectory_ref2, measurements_struct2);

%% Iterate (3)
% 
% % Time span (fixed duration, e.g., 50 days)        
% t_span = trackingData.time;             % [s]
% options = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-16);
% flag_STM = 1;
% eyeI = eye(10);
% 
% % Integrate trajectory
% [t_ref3, state_ref3] = ode113(@(t, state) ...
%     flyby_ODE(t, state, JD0, B, mu_E, mu_S, AMR, Pphi, AU, flag_STM), ...
%     t_span, [results.state_corrected_hist(:, 1); eyeI(:)], options);
% 
% % STM reshape
% STM_ref3 = reshape(state_ref3(:, 11:end), [], 10, 10);
% 
% % Preallocate structs
% trajectory_ref3 = struct('time', {}, 'state', {}, 'STM', {}, ...
%     'parameters', {}, 'function_handle', {});
% for i = 1:length(t_ref)
%     trajectory_ref3(i).time = t_ref3(i);
%     trajectory_ref3(i).state = state_ref3(i,1:6);
%     trajectory_ref3(i).parameters = state_ref3(i,7:10);
%     trajectory_ref3(i).STM = STM_ref3(i,:,:);
%     trajectory_ref3(i).function_handle = @(t, x) ...
%         flyby_ODE(t, x, JD0, B, mu_E, mu_S, AMR, Pphi, AU, flag_STM);
% end
% 
% % Compute Measurements Using the Radiometric Function
% measurements_struct3 = radiometric_measurement(trackingData, state_ref3', ...
%     station_ecef0, R);
% 
% % Run the filter
% % filename_suffix = 'BatchSrif3';
% results = lkf(trajectory_ref3, measurements_struct3, P0);
% results = rts_smoother(results, trajectory_ref3, measurements_struct3);

%% Plot filtering results

% Time vector in hours
measurement_times_hours = [measurements_struct.time] / 3600;

% Initial a priori state
x_apriori = [r0; v0; par0];  % 7x1

% Compute ±3σ bounds from EKF covariance at each step
n_states = size(results.P_hist, 1);
n_times = length(measurement_times_hours);
sigma_bounds = zeros(n_states, n_times);

for t = 1:n_times
    sigma_bounds(:, t) = 3 * sqrt(diag(results.P_hist(:, :, t)));
end

% Plot: semilogy with different colors and improved visuals
figure('Position', [100, 100, 1400, 900]);

% State labels
labels = {'$x$', '$y$', '$z$', ...
          '$\dot{x}$', '$\dot{y}$', '$\dot{z}$', ...
          '$C_R$', '$w_x$', '$w_y$', '$w_z$'};

% Y-axis group labels
y_labels = {'Position Error (km)', ...
            'Velocity Error (km/s)', ...
            'CR Error', ...
            'Empirical Accel. Error (km/s$^2$)'};

% Distinct colors
cmap = lines(10);  % nice set of distinct colors

for i = 1:10
    subplot(4, 3, i);
    
    semilogy(measurement_times_hours, sigma_bounds(i, :), ...
        '.', 'Color', cmap(i,:), 'MarkerSize', 10, ...
        'DisplayName', ['$3\sigma$: ', labels{i}]);
    
    xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
    
    % Choose correct y-axis label
    if i <= 3
        ylabel(y_labels{1}, 'Interpreter', 'latex', 'FontSize', 14);
    elseif i <= 6
        ylabel(y_labels{2}, 'Interpreter', 'latex', 'FontSize', 14);
    elseif i == 7
        ylabel(y_labels{3}, 'Interpreter', 'latex', 'FontSize', 14);
    else
        ylabel(y_labels{4}, 'Interpreter', 'latex', 'FontSize', 14);
    end

    title(labels{i}, 'Interpreter', 'latex', 'FontSize', 16);
    grid on;
    legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
end

% Prep Info
range_res     = results.postfit_residuals(1, :);
rangerate_res = results.postfit_residuals(2, :);

range_mean = mean(range_res);
range_std  = std(range_res);

rangerate_mean = mean(rangerate_res);
rangerate_std  = std(rangerate_res);

% Post-fits Residuals Plot
figure('Position', [100, 100, 1400, 800]);
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
histogram(range_res, 40, 'Orientation', 'horizontal', ...
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
ylabel('Range-Rate Residual (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
title('Range-Rate Residuals vs. Time', 'Interpreter', 'latex', 'FontSize', 16);
ylim(ylims_rangerate);
grid on;
% Right: Histogram
nexttile(4);
histogram(rangerate_res, 40, 'Orientation', 'horizontal', ...
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


%% Final Solutions
% 
% % Time span (fixed duration, e.g., 50 days)        
% t_span = trackingData.time;             % [s]
% options = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-16);
% flag_STM = 1;
% eyeI = eye(10);
% 
% % Integrate trajectory
% [t_ref4, state_ref4] = ode113(@(t, state) ...
%     flyby_ODE(t, state, JD0, B, mu_E, mu_S, AMR, Pphi, AU, flag_STM), ...
%     t_span, [results.state_corrected_hist(:, 1); eyeI(:)], options);
% 
% 
% % Compute Measurements Using the Radiometric Function
% measurements_struct4 = radiometric_measurement(trackingData, state_ref4', ...
%     station_ecef0, R);
% 
% 
% % Extract Residuals Per Station for Plotting
% stations = unique(trackingData.stationID); % Unique station IDs
% num_stations = length(stations);
% 
% % Initialize cell arrays to store residuals and times per station
% range_residuals = cell(num_stations, 1);
% range_rate_residuals = cell(num_stations, 1);
% time_per_station = cell(num_stations, 1);
% 
% % Extract residuals for each station
% for s = 1:num_stations
%     station_id = stations(s);
% 
%     % Find indices corresponding to this station
%     station_indices = find(trackingData.stationID == station_id);
% 
%     % Store time and residuals
%     time_per_station{s} = trackingData.time(station_indices); % Time in seconds
%     range_residuals{s} = cellfun(@(x) x(1), {measurements_struct4(station_indices).residual});
%     range_rate_residuals{s} = cellfun(@(x) x(2), {measurements_struct4(station_indices).residual});
% end
% 
% % Plot Residuals 
% gca123 = figure('Position', [100, 100, 1400, 900]);
% subplot(2, 1, 1);
% hold on;
% for s = 1:num_stations
%     if ~isempty(range_residuals{s}) && length(time_per_station{s}) == length(range_residuals{s})
%         plot(time_per_station{s} / 3600, range_residuals{s}, 'o', ...
%              'DisplayName', sprintf('DSS%d', stations(s)), 'MarkerSize', 8, 'LineWidth', 1.2);
%     end
% end
% xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 16);
% ylabel('Range Residual (km)', 'Interpreter', 'latex', 'FontSize', 16);
% legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');
% grid on;
% hold off;
% subplot(2, 1, 2);
% hold on;
% for s = 1:num_stations
%     if ~isempty(range_rate_residuals{s}) && length(time_per_station{s}) == length(range_rate_residuals{s})
%         plot(time_per_station{s} / 3600, range_rate_residuals{s}, 'o', ...
%              'DisplayName', sprintf('DSS%d', stations(s)), 'MarkerSize', 8, 'LineWidth', 1.2);
%     end
% end
% xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 16);
% ylabel('Range-Rate Residual (km/s)', 'Interpreter', 'latex', 'FontSize', 16);
% legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');
% grid on;
% hold off;
% 


%% Helper Functions Dynamics

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

function dstate = flyby_ODE(~, state, ~, ~, mu_E, mu_S, AM_ratio, Pphi, AU, flag_STM)
    % state = [r_sc(3); v_sc(3); Cr; r_sun(3); v_sun(3);
    %          Phi(:) if flag_STM]
    
    %--- dimensions ---
    n_sc  = 3+3+1;   % [r_sc; v_sc; Cr]
    n_sun = 3+3;     % [r_sun; v_sun]
    n     = n_sc + n_sun;
    
    if flag_STM
        dstate = zeros(n + n^2,1);
    else
        dstate = zeros(n,1);
    end

    %--- unpack state ---
    r      = state(1:3);       % sc pos wrt Earth
    v      = state(4:6);       % sc vel wrt Earth
    Cr     = state(7);         % reflectivity coeff
    r_sun  = state(8:10);      % Sun pos wrt Earth
    v_sun  = state(11:13);     % Sun vel wrt Earth

    %--- geometry ---
    r_norm      = norm(r);
    r_sun_norm  = norm(r_sun);
    r_rel       = r_sun - r;
    r_rel_norm  = norm(r_rel);

    %--- accelerations on spacecraft ---
    a_earth = -mu_E * r / r_norm^3;
    a_sun   = mu_S * (r_rel/r_rel_norm^3 - r_sun/r_sun_norm^3);
    a_srp   = - AM_ratio * Cr * Pphi * (AU^2) * (r_rel / r_rel_norm^3);
    a_tot   = a_earth + a_sun + a_srp;

    %--- fill derivatives ---
    dstate(1:3)   = v;           % dr_sc/dt
    dstate(4:6)   = a_tot;       % dv_sc/dt
    dstate(7)     = 0;           % Cr constant
    dstate(8:10)  = v_sun;       % dr_sun/dt = Earth vel (neg Sun vel)
    dstate(11:13) = -mu_S * r_sun / r_sun_norm^3;  % dv_sun/dt

    %--- STM propagation if requested ---
    if flag_STM
        % reshape STM
        Phi = reshape(state(n+1:end), n, n);
        A   = zeros(n);

        % 1) Spacecraft blocks
        A(1:3,4:6) = eye(3);

        % d a_earth / dr
        dadr_E = mu_E*(3*(r*r.')/r_norm^5 - eye(3)/r_norm^3);

        % d a_sun / dr_sc  (r_rel = r_sun - r => ∂r_rel/∂r = -I)
        dadr_S = -mu_S*(eye(3)/r_rel_norm^3 - 3*(r_rel*r_rel.')/r_rel_norm^5);

        % d a_srp / dr_sc  (same form as sun but scaled)
        scale2   = AM_ratio * Cr * Pphi * AU^2;
        dadr_SRP =  scale2*(eye(3)/r_rel_norm^3 - 3*(r_rel*r_rel.')/r_rel_norm^5);

        A(4:6,1:3) = dadr_E + dadr_S + dadr_SRP;

        % ∂a_tot/∂Cr
        A(4:6,7) = - AM_ratio * Pphi * AU^2 * (r_rel / r_rel_norm^3);

        % ∂a_tot/∂r_sun
        dadw_sun     = mu_S*( eye(3)/r_rel_norm^3 ...
                             - 3*(r_rel*r_rel.')/r_rel_norm^5 ...
                             - ( eye(3)/r_sun_norm^3 ...
                                 - 3*(r_sun*r_sun.')/r_sun_norm^5 ) );
        dadw_srp =  scale2*( eye(3)/r_rel_norm^3 ...
                    - 3*(r_rel*r_rel.')/r_rel_norm^5 );
        A(4:6,8:10)  = dadw_sun + dadw_srp;

        % 2) Sun‐ephemeris blocks
        A(8:10,11:13)   = eye(3);  % dr_sun/dt = v_sun
        A(11:13,8:10)   = -mu_S*( ...
                             eye(3)/r_sun_norm^3 ...
                           - 3*(r_sun*r_sun.')/r_sun_norm^5 );

        % propagate STM
        dPhi = A * Phi;
        dstate(n+1:end) = dPhi(:);
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
        d_measurements_dX = [d_rho_dR, d_rho_dV, 0, 0, 0, 0, 0, 0, 0;
                             d_rho_dot_dR, d_rho_dot_dV, 0, 0, 0, 0, 0, 0, 0];

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

function results = lkf(trajectory_ref, sorted_measurements, P0, settings_PN)
    % Initialization
    T = length(sorted_measurements); % Number of measurements
    n = size(P0, 1);  % State dimension
    max_m = max(arrayfun(@(s) length(s.residual), sorted_measurements));

    state_deviation_hist = zeros(n, T);
    state_corrected_hist = zeros(n, T);
    Pbar_hist = zeros(n, n, T);
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
        x_ref = [trajectory_ref(traj_idx).state(1:6)'; trajectory_ref(traj_idx).parameters(:)];

        % Compute STM transition matrix
        STM_dt = STM_t / STM_tm; 

        % Extract measurement data
        prefit_residual = sorted_measurements(t).residual;
        H_tilde = sorted_measurements(t).partials.wrt_X;
        R = sorted_measurements(t).covariance;

        % Compute Process Noise
        if nargin > 3
            dt = 0;  % Default value for the first iteration
            if t > 1
                dt = measurement_times(t) - measurement_times(t-1); % Time step
            end
            Q_9x9 = calculate_Q_discrete(dt, x_ref(1:6), settings_PN);  % 9x9
            Q_disc = zeros(n);  % Final 10x10 output
            
            % Fill in Q(1:6,1:6) – position and velocity
            Q_disc(1:6, 1:6) = Q_9x9(1:6, 1:6);

        else 
            Q_disc = zeros(n);
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
        Pbar_hist(:, :, t) = P_pred;
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
        'Pbar_hist', Pbar_hist, ...
        'postfit_residuals', postfit_residuals ...
    );
end


function smoothed_results = rts_smoother(lkf_results, trajectory_ref, sorted_measurements)
    % Extract dimensions
    T = size(lkf_results.state_corrected_hist, 2);
    n = size(lkf_results.state_corrected_hist, 1);
    
    % Initialize smoothed state and covariance
    state_smoothed_hist = zeros(n, T);
    state_dev_smoothed_hist = zeros(n, T);
    P_smoothed_hist = zeros(n, n, T);
    postfit_residuals_smoothed = NaN(size(lkf_results.postfit_residuals));
    
    % Final smoothed states are the same as filtered
    state_smoothed_hist(:, end) = lkf_results.state_corrected_hist(:, end);
    state_dev_smoothed_hist(:, end) = lkf_results.state_deviation_hist(:, end);
    P_smoothed_hist(:, :, end) = lkf_results.P_hist(:, :, end);
    postfit_residuals_smoothed(:, end) = lkf_results.postfit_residuals(:, end);
    
    % Extract trajectory times for reference
    trajectory_times = [trajectory_ref.time];
    measurement_times = [sorted_measurements.time];
    [~, traj_indices] = ismember(measurement_times, trajectory_times);

    % Progress tracking
    fprintf('RTS Smoother Progress: 0%%');
    
    % Backward smoothing
    for t = T-1:-1:1
        traj_idx_next = traj_indices(t+1);
        traj_idx = traj_indices(t);
        STM_tk_tkp1 = squeeze(trajectory_ref(traj_idx_next).STM) / ...
            squeeze(trajectory_ref(traj_idx).STM);
        
        % Retrieve filter estimates
        x_filt = lkf_results.state_deviation_hist(:, t);
        P_filt = lkf_results.P_hist(:, :, t);
        P_bar_next = lkf_results.Pbar_hist(:, :, t+1);
        
        % Compute smoothing gain
        S = P_filt * STM_tk_tkp1' / P_bar_next;
        
        % Smoothed estimate
        x_smoothed = x_filt + S * (state_dev_smoothed_hist(:, t+1) - STM_tk_tkp1 * x_filt);
        P_smoothed = P_filt + S * (P_smoothed_hist(:, :, t+1) - P_bar_next) * S';
        
        % Compute smoothed postfit residuals
        postfit_residuals_smoothed(:, t) = sorted_measurements(t).residual - sorted_measurements(t).partials.wrt_X * x_smoothed;
        
        % Store results
        state_dev_smoothed_hist(:, t) = x_smoothed;
        state_smoothed_hist(:, t) = [trajectory_ref(traj_idx).state(1:6)'; ...
            trajectory_ref(traj_idx).parameters(:)] + x_smoothed;
        P_smoothed_hist(:, :, t) = P_smoothed;

        % Print progress
        fprintf('\b\b\b\b%3d%%', round(((T - t) / T) * 100));
    end

    fprintf('\nRTS Smoother Completed!\n');
    
    % Store all outputs in a struct
    smoothed_results = struct(...
        'state_corrected_hist', state_smoothed_hist, ...
        'state_deviation_hist', state_dev_smoothed_hist, ...
        'P_hist', P_smoothed_hist, ...
        'postfit_residuals', postfit_residuals_smoothed ...
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

