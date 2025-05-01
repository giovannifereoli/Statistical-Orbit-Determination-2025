%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         January 22, 2024
% CLASS:        ASEN 6080: StatOD
% INSTRUCTOR:   Prof. Jay W. McMahon
% ASSIGNMENT:   Project 2 - Part 2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% TODO: try mu, try bias, try GS, try maneuvers
% TODO: why 8 months and 8 oscillation like increasing toward sun?

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
Cr          = 1.0;                            % SRP reflectivity coefficient
c = 299792458;                                % Speed of light [m/s]
Pphi = (1357 / c) * 1e3;                      % SRP pressure in N/km^2
AU = 149597870.7;                             % Astro Units in km

% Initial Conditions (pos/vel in ECI)
r0 = [-274096770.76544; -92859266.44990661; -40199493.6677441];                    % km
v0 = [32.6704564599943; -8.93838913761049; -3.87881914050316];                     % km/s
par0 = Cr; 

% Initial full state vector: [x; v; Cr; STM]
state0_true = [r0; v0; par0; reshape(eye(7), [], 1)];

% Define Measurement Noise Covariance Matrix R
sigma_rho = 5 * 1e-3; % Standard deviation of range measurement (km)
sigma_rho_dot = 0.5 * 1e-6; % Standard deviation of range rate (km/s)
R = diag([sigma_rho^2, sigma_rho_dot^2]); % Covariance matrix

% DSN Station geodetic data (lat, lon in deg, alt in km)
stations = struct();
stations.DSS34 = struct('lat', -35.398333, 'lon', 148.981944,   'alt', 0.691750);     % Canberra
stations.DSS65 = struct('lat',  40.427222, 'lon', 355.749444,   'alt', 0.834539);     % Madrid (positive longitude for Part 3)
stations.DSS13 = struct('lat',  35.247164, 'lon', 243.205000,   'alt', 1.07114904);   % Goldstone
station_ecef0 = latlon_to_ecef(stations);

%% Run LSQNONLIN - 1

% % Initialize
% theta0 = [r0; v0; 1.0; mu_E; mu_S];
% resfun = @(theta) residual_vector_weighted(theta, trackingData, station_ecef0, AMR, Pphi, AU, JD0);
% opts = optimoptions('lsqnonlin', 'Display', 'iter', 'Algorithm','interior-point', ...
%     'MaxIterations', 100);
% [theta_est, resnorm, residuals] = lsqnonlin(resfun, theta0, [], [], opts);
% 
% % Function to be Minimized
% function r_weighted = residual_vector_weighted(theta, trackingData, station_ecef0, AMR, Pphi, AU, JD0)
% 
%     % Unpack parameters
%     r0   = theta(1:3);
%     v0   = theta(4:6);
%     Cr   = theta(7);
%     mu_E = theta(8);
%     mu_S = theta(9);
% 
%     % Initial state
%     STM0 = eye(7);
%     x0 = [r0; v0; Cr; STM0(:)];
% 
%     % Propagate trajectory
%     t_span = trackingData.time;
%     options = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-16);
%     [~, state_ref] = ode113(@(t, state) ...
%         flyby_ODE(t, state, JD0, mu_E, mu_S, AMR, Pphi, AU, 1), ...
%         t_span, x0, options);
% 
%     % Measurement model
%     measurements_struct = radiometric_measurement(trackingData, state_ref', station_ecef0, eye(2));
% 
%     % Stack weighted residuals
%     r_weighted = [];
%     for i = 1:length(measurements_struct)
%         W = inv(measurements_struct(i).covariance);
%         r = measurements_struct(i).residual;
%         r_weighted = [r_weighted; sqrtm(W) * r / 1000];  % decorrelate
%     end
% end
% 
% % Save Results
% r0   = theta_est(1:3);
% v0   = theta_est(4:6);
% Cr   = theta_est(7);
% mu_E = theta_est(8);
% mu_S = theta_est(9);
% state0_true = [r0; v0; Cr; reshape(eye(7), [], 1)];

%% 

% Initial ECEF positions (3 stations × 3 coords)
ecef_init = [station_ecef0("DSS34"); station_ecef0("DSS65"); station_ecef0("DSS13")];
theta0 = [r0; v0; 1.0; mu_E; mu_S; ecef_init(:)];
resfun = @(theta) residual_vector_weighted2(theta, trackingData, AMR, Pphi, AU, JD0);
opts = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxIterations', 10);
[theta_est, resnorm, residuals] = lsqnonlin(resfun, theta0, [], [], opts);

function r_weighted = residual_vector_weighted2(theta, trackingData, AMR, Pphi, AU, JD0)
    % Unpack trajectory + params
    r0   = theta(1:3);
    v0   = theta(4:6);
    Cr   = theta(7);
    mu_E = theta(8);
    mu_S = theta(9);
    
    % Unpack station positions
    ecef_flat = theta(10:end);
    station_ecef0 = containers.Map();
    station_ecef0("DSS34") = ecef_flat(1:3);
    station_ecef0("DSS65") = ecef_flat(4:6);
    station_ecef0("DSS13") = ecef_flat(7:9);

    % Initial state
    STM0 = eye(7);
    x0 = [r0; v0; Cr; STM0(:)];

    % Propagate trajectory
    t_span = trackingData.time;
    options = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-16);
    [~, state_ref] = ode113(@(t, state) ...
        flyby_ODE(t, state, JD0, mu_E, mu_S, AMR, Pphi, AU, 1), ...
        t_span, x0, options);

    % Measurement model
    measurements_struct = radiometric_measurement(trackingData, state_ref', station_ecef0, eye(2));

    % Stack weighted residuals
    r_weighted = [];
    for i = 1:length(measurements_struct)
        W = inv(measurements_struct(i).covariance);
        r = measurements_struct(i).residual;
        r_weighted = [r_weighted; sqrtm(W) * r];  % decorrelate
    end
end

%% 
% Re-initialize
r0      = theta_est(1:3);
v0      = theta_est(4:6);
Cr      = theta_est(7);
mu_E    = theta_est(8);
mu_S    = theta_est(9);
ecef_34 = theta_est(10:12);
ecef_65 = theta_est(13:15);
ecef_13 = theta_est(16:18);
station_ecef0 = containers.Map();
station_ecef0("DSS34") = ecef_34;
station_ecef0("DSS65") = ecef_65 ;
station_ecef0("DSS13") = ecef_13;
state0_true = [r0; v0; Cr; reshape(eye(7), [], 1)];

%%
% 
% ecef_init = [station_ecef0("DSS34"); station_ecef0("DSS65"); station_ecef0("DSS13")];
% DV_guess = [0; 0; 0];                   % Δv [km/s]
% tb_guess = trackingData.time(round(end/2));  % Guess burn time in the middle of arc
% theta0 = [r0; v0; 1.0; mu_E; mu_S; DV_guess; tb_guess; ecef_init(:)];
% resfun = @(theta) residual_vector_weighted3(theta, trackingData, AMR, Pphi, AU, JD0);
% opts = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxIterations', 100);
% [theta_est, resnorm, residuals] = lsqnonlin(resfun, theta0, [], [], opts);
% 
% function r_weighted = residual_vector_weighted3(theta, trackingData, AMR, Pphi, AU, JD0)
% 
%     % === Unpack parameters ===
%     r0   = theta(1:3);
%     v0   = theta(4:6);
%     Cr   = theta(7);
%     mu_E = theta(8);
%     mu_S = theta(9);
%     DV   = theta(10:12);      % Δv [km/s]
%     t_burn = theta(13);       % Burn time [s]
% 
%     ecef_flat = theta(14:end);
%     station_ecef0 = containers.Map();
%     station_ecef0("DSS34") = ecef_flat(1:3);
%     station_ecef0("DSS65") = ecef_flat(4:6);
%     station_ecef0("DSS13") = ecef_flat(7:9);
% 
%     % === Initial state ===
%     STM0 = eye(7);
%     x0 = [r0; v0; Cr; STM0(:)];
% 
%     t_full = trackingData.time;
%     options = odeset('RelTol', 2.22045e-14, 'AbsTol', 2.22045e-16);
% 
%     % === Pre-burn propagation ===
%     t_pre = t_full(t_full <= t_burn);
%     if isempty(t_pre), t_pre = [t_full(1); t_burn]; end  % at least two points
% 
%     [t_pre_out, state_pre] = ode113(@(t, state) ...
%         flyby_ODE(t, state, JD0, mu_E, mu_S, AMR, Pphi, AU, 1), ...
%         [t_pre; t_burn], x0, options);
% 
%     % Interpolate to burn time
%     state_burn = interp1(t_pre_out, state_pre, t_burn, 'linear')';
%     state_burn(4:6) = state_burn(4:6) + DV;  % Apply impulse
% 
%     % === Post-burn propagation ===
%     t_post = t_full(t_full > t_burn);
%     if isempty(t_post), t_post = t_burn + [1]; end
% 
%     JD_burn = JD0 + t_burn / 86400;
%     [~, state_post] = ode113(@(t, state) ...
%         flyby_ODE(t, state, JD_burn, mu_E, mu_S, AMR, Pphi, AU, 1), ...
%         t_post, state_burn, options);
% 
%     % === Combine ===
%     state_all = [state_pre(1:end-1,:); state_burn'; state_post];
%     t_all = [t_pre; t_post];
% 
%     % === Measurement model ===
%     measurements_struct = radiometric_measurement(trackingData, state_all', station_ecef0, eye(2));
% 
%     % === Stack weighted residuals ===
%     r_weighted = [];
%     for i = 1:length(measurements_struct)
%         W = inv(measurements_struct(i).covariance);
%         r = measurements_struct(i).residual;
%         r_weighted = [r_weighted; sqrtm(W) * r];
%     end
% end
% 

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

% % Initialize
% sigma_pos = 100;       % Std for position states (km)
% sigma_vel = 1e-1;      % Std for velocity states (km/s)
% sigma_Cr = 0.1;        % Std for parameter states (-)
% P0 = diag([sigma_pos^2, sigma_pos^2, sigma_pos^2, ...
%            sigma_vel^2, sigma_vel^2, sigma_vel^2, sigma_Cr^2]);
% 
% % Run the filter
% % filename_suffix = 'BatchSrif';
% results = batch_srif_filter(trajectory_ref, measurements_struct, P0);

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
