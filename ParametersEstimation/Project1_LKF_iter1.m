%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         January 22, 2024
% CLASS:        ASEN 6080: StatOD
% INSTRUCTOR:   Prof. Jay W. Mcmahon
% ASSIGNMENT:   Project 1 - LKF - Iter1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close;

%% Parse Data

% Define the filename (update with the actual path if needed)
filename = 'project.txt';  % Replace with the actual filename

% Read the data from the file
data = readmatrix(filename);

% Extract columns into variables
time = data(:, 1);          % Time (seconds)
stationID = data(:, 2);     % Station ID
range = data(:, 3) / 1e3;         % Range (km)
rangeRate = data(:, 4) / 1e3;     % Range rate (km/s)

% Optionally, store in a struct for easy access
trackingData.time = time;
trackingData.stationID = stationID;
trackingData.range = range;
trackingData.rangeRate = rangeRate;

%% Create Reference Trajectory

% Constants
coeffs = [0, 1.082626925638815 * 1e-3]; % Zonal harmonic coefficients (i.e., [J1, J2])
R = 6378.1363; % Earth radius 
mu = 398600.4415;  % Gravitational parameter 
l_max = 2; % Zonal Harmonics Terms 
A = 3 * 1e-6; % Cross-Sectional Area
m = 970; % Mass
CD = 2.0; % Drag Coefficient
param_meas = 9;

% Define the station positions in ECEF at t0
station_ecef0 = containers.Map(...
    {'101', '337', '394'}, ...
    {[-5127.5100; -3794.1600; 0.0], ...
     [3860.9100; 3238.4900; 3898.0940], ...
     [549.5050; -1380.8720; 6182.1970]});

% Initial state (reference orbit)
state0 = [757.7; 5222.607; 4851.5; 2.21321; 4.67834; -5.3713]; % (km, km/s)
state0 = [state0; mu; coeffs'; CD;  station_ecef0('101'); ...
     station_ecef0('337'); station_ecef0('394'); ...
     reshape(eye(19), [], 1)]; % Adding mu, J1, J2, CD

% Initialize ODEs
[acc_SPH, A_SPH] = precomputeZonalSPH(R, l_max);
[acc_Drag, A_Drag] = precomputeDrag(A, m, R);

% Simulation parameters
T_orbit = 18340; % Orbital period (s)
t_span = trackingData.time;      
options = odeset('RelTol', 2.22045 * 1e-14, 'AbsTol',  2.22045 * 1e-14);

% Integrate reference trajectory
[t_ref, state_ref] = ode45(@(t, state) ...
     sc_ODE(t, state, coeffs, mu, l_max, CD, acc_SPH, A_SPH, acc_Drag, ...
     A_Drag, param_meas), t_span, state0, options);

% Plot Orbit
gca1 = figure(1);
hold on;  
plot3(state_ref(:, 1), state_ref(:, 2), state_ref(:, 3), ...
    'r.', 'MarkerSize', 10);   
plot_earth();  
axis equal;
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('z (km)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Reference Trajectory at $t_{obs}$', 'Earth', ...
    'Location', 'northwest', 'Interpreter', 'latex','FontSize', 14);
grid on;
view(3);

%% Construct Reference Trajectory Dictionary

% Extract the STM matrices for [r, v]
size = (-1 + sqrt(1 + 4 * size(state0, 1))) / 2;
STM_ref = reshape(state_ref(:, (size + 1):end), [], size, size);

% Preallocate structure arrays
num_steps_ref = length(t_ref);
trajectory_ref = struct('time', cell(1, num_steps_ref), 'state', ...
    cell(1, num_steps_ref), 'STM', cell(1, num_steps_ref), ...
    'parameters', cell(1, num_steps_ref),'function_handle', cell(1, num_steps_ref));

% Fill the reference trajectory structure
for i = 1:num_steps_ref
    trajectory_ref(i).time = t_ref(i);
    trajectory_ref(i).state = state_ref(i, 1:size); % Extract everything
    trajectory_ref(i).parameters = state_ref(i, 7:size); % Extract mu, J1, J2, CD
    trajectory_ref(i).STM = STM_ref(i, :, :);
    trajectory_ref(i).function_handle = @(t, x)  sc_ODE(t, x, coeffs, ...
        mu, l_max, CD, acc_SPH, A_SPH, acc_Drag, A_Drag);
end

%% Create Measurement Dictionary 

% Define Measurement Noise Covariance Matrix R
sigma_rho = 1e-5; % Standard deviation of range measurement (km)
sigma_rho_dot = 1e-6; % Standard deviation of range rate (km/s)
R = 1e-12 * diag([sigma_rho^2, sigma_rho_dot^2]); % Covariance matrix

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
             'DisplayName', sprintf('GS %d', stations(s)), 'MarkerSize', 8, 'LineWidth', 1.2);
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
             'DisplayName', sprintf('GS %d', stations(s)), 'MarkerSize', 8, 'LineWidth', 1.2);
    end
end
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Range-Rate Residual $\cdot 1e-3$ (km/s)', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(gca123, 'Prefit_LKF_Iter1.pdf', 'ContentType',...
 'vector', 'Resolution', 1000);


%% Create Filter Data Set (aka, improved Measurement Dictionary)

% Create Filter Data Set
params_dyn = 4;
params_meas = 9;
filter_dataset = create_filter_dataset(measurements_struct, ...
    params_dyn, params_meas);

%% Estimation - LKF

% Define initial covariance matrix P0
P0 = diag([1, 1, 1, 1, 1, 1, 1e11, 1e-16, 1e6, 1e6, ...
           1e-16, 1e-16, 1e-16, 1, 1, 1, 1, 1, 1]);

% Run the Filter (LKF, Batch, or EKF)
results = lkf(trajectory_ref, filter_dataset, P0);
filename_suffix = 'LKFIter1';

% Print Results
compute_rms_errors(results, filter_dataset);

%% Plot filtering results

% Initialize
measurement_times = [filter_dataset.time];
measurement_times_hours = measurement_times / 3600; 

% Deviation from Reference in Time for 19 states
% Define state labels with units
state_labels = {'$r_x$', '$r_y$ ', '$r_z$ ', ...
                '$v_x$', '$v_y$', '$v_z$', ...
                '$\mu$', '$J_1$', '$J_2$', '$C_D$', ...
                '$r_{x,GS1}$', '$r_{y,GS1}$', '$r_{z,GS1}$', ...
                '$r_{x,GS2}$', '$r_{y,GS2}$', '$r_{z,GS2}$', ...
                '$r_{x,GS3}$', '$r_{y,GS3}$', '$r_{z,GS3}$'};
position_indices = 1:3; % rx, ry, rz
velocity_indices = 4:6; % vx, vy, vz
dynamic_indices = 7:10; % mu, J1, J2, CD
measurement_indices = 11:19; % Ground station parameters
gca10 = figure(10);
set(gcf, 'Position', [100, 100, 1400, 900]); 
colors = lines(19);
subplot(2,2,1);
hold on;
for i = position_indices
    scatter(measurement_times_hours, results.state_deviation_hist(i, :), ...
        20, colors(i, :), 'filled', 'DisplayName', state_labels{i});
end
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Position Deviation (km)', 'Interpreter', 'latex', 'FontSize', 14);
legend('show', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
grid on;
hold off;
subplot(2,2,2);
hold on;
for i = velocity_indices
    scatter(measurement_times_hours, results.state_deviation_hist(i, :), ...
        20, colors(i, :), 'filled', 'DisplayName', state_labels{i});
end
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Velocity Deviation (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
legend('show', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
grid on;
hold off;
subplot(2,2,3);
hold on;
for i = dynamic_indices
    scatter(measurement_times_hours, results.state_deviation_hist(i, :), ...
        20, colors(i, :), 'filled', 'DisplayName', state_labels{i});
end
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Deviation $\cdot 1e4$ (-)', 'Interpreter', 'latex', 'FontSize', 14);
legend('show', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
grid on;
hold off;
subplot(2,2,4);
hold on;
for i = measurement_indices
    scatter(measurement_times_hours, results.state_deviation_hist(i, :), ...
        10, colors(i, :), 'filled', 'DisplayName', state_labels{i});
end
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Deviation (km)', 'Interpreter', 'latex', 'FontSize', 14);
legend('show', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
grid on;
hold off;
exportgraphics(gca10, sprintf('StateDeviationsSubplots_%s.pdf', filename_suffix), ...
     'ContentType','image', 'Resolution', 1000);


% Plot range and range-rate post-fit residuals 
gca6 = figure(6);
set(gcf, 'Position', [100, 100, 1400, 900]);
subplot(2,1,1);
hold on;
scatter(measurement_times_hours, results.postfit_residuals(1, :), 35, 'b',...
    'filled', 'DisplayName', 'Range Residual', 'MarkerEdgeColor', 'k'); 
fill([measurement_times_hours, fliplr(measurement_times_hours)], ...
    [3 * sigma_rho * ones(1, length(measurement_times_hours)), ...
    fliplr(-3 * sigma_rho * ones(1, length(measurement_times_hours)))], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '$\pm3\sigma$');
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Range Residual $\cdot 1e-4$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
legend('show', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
hold off;
ylim([-100 * sigma_rho, 100 * sigma_rho]);
xlim([measurement_times_hours(1), measurement_times_hours(end)]);
subplot(2,1,2);
hold on;  
scatter(measurement_times_hours, results.postfit_residuals(2, :), 35, ...
    'r', 'filled', 'DisplayName', 'Range Rate Residual', 'MarkerEdgeColor', 'k'); 
fill([measurement_times_hours, fliplr(measurement_times_hours)], ...
    [3 * sigma_rho_dot * ones(1, length(measurement_times_hours)), ...
    fliplr(-3 * sigma_rho_dot * ones(1, length(measurement_times_hours)))], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '$\pm3\sigma$');
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14, 'FontSize', 14);
ylabel('Range-Rate Residual $\cdot 1e-5$ (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
legend('show', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
hold off;
ylim([-100 * sigma_rho_dot, 100 * sigma_rho_dot]);
xlim([measurement_times_hours(1), measurement_times_hours(end)]);
exportgraphics(gca6, sprintf('PostFit_%s.pdf', filename_suffix), ...,
     'ContentType','image', 'Resolution', 1000);

%% Analysis Numerical Issues - Plot Trace for [r, v]

% Initialize matrices
trace_P_pos = zeros(1, length(measurement_times));  % Trace of position covariance
trace_P_vel = zeros(1, length(measurement_times));  % Trace of velocity covariance

T_meas = length(measurement_times);
for t = 1:T_meas
    % Extract covariance matrices for position and velocity
    P_pos = results.P_hist(1:3, 1:3, t);  % Position covariance (rx, ry, rz)
    P_vel = results.P_hist(4:6, 4:6, t);  % Velocity covariance (vx, vy, vz)

    % Compute traces
    trace_P_pos(t) = trace(P_pos);
    trace_P_vel(t) = trace(P_vel);
end

%Plot
gca11 = figure(11);
set(gca11, 'Position', [100, 100, 1400, 600]); 
subplot(2,1,1);
scatter(measurement_times_hours, trace_P_pos, 30, 'b', 'Filled', 'MarkerEdgeColor', 'k');
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Trace of Covariance for $\mathbf{r}$ (-)', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
subplot(2,1,2);
scatter(measurement_times_hours, trace_P_vel, 30, 'r', 'Filled', 'MarkerEdgeColor', 'k');
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Trace of Covariance for $\mathbf{v}$ (-)', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
exportgraphics(gca11, sprintf('Trace_P_%s.pdf', filename_suffix), ...
    'ContentType', 'vector', 'Resolution', 1000);

%% Re-iterate Again

% Save correction
state0_deviation= results.state0_deviation;
save('state0_deviation_iter1.mat', 'state0_deviation');

% Save solution
results_LKF_Iter1 = results;
save('results_LKF_Iter1.mat', 'results_LKF_Iter1');

% Save Time
save('meas_times.mat', 'time');

%% Print out results

% Extract the final covariance matrix
P_final_lkf = results.P_hist(:, :, end);

% Extract diagonal elements (variances)
variances = diag(P_final_lkf);

% Compute standard deviations
std_devs = sqrt(variances);

% Extract the final state values
state_final_lkf = results.state_corrected_hist(:, end);

% Define state labels
state_labels = {'$r_x$', '$r_y$', '$r_z$', ...
                '$v_x$', '$v_y$', '$v_z$', ...
                '$\mu$', '$J_1$', '$J_2$', '$C_D$', ...
                '$r_{x,GS1}$', '$r_{y,GS1}$', '$r_{z,GS1}$', ...
                '$r_{x,GS2}$', '$r_{y,GS2}$', '$r_{z,GS2}$', ...
                '$r_{x,GS3}$', '$r_{y,GS3}$', '$r_{z,GS3}$'};

% Print results in a structured way
fprintf('\n%-15s %-20s %-20s\n', 'State', 'Mean Value', 'Final Std Dev');
fprintf('------------------------------------------------------------\n');
for i = 1:length(state_labels)
    fprintf('%-15s %-25.14f %-25.14f\n', state_labels{i}, state_final_lkf(i), std_devs(i));
end


%% Helper Functions

function R = DCM_eciToEcef(t)
    % Earth's rotation rate 
    % omega_earth = 2 * pi / (24 * 3600);
    omega_earth = 7.2921158553 * 1e-5; % OSS: Given my assignment.
    
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

function [acceleration_func, A_func] = precomputeDrag(A, m, REarth)
    % Constants
    omega_E = 7.2921158553 * 1e-5; % Earth's angular velocity (rad/s)
    rho0 = 3.614 * 1e-4; % Reference Density(kg/km^3)
    r0 = (700.00 + REarth); % Reference Radius (km)
    H = 88.6670; % Scale Height(km)

    % Symbolic definitions
    syms x y z vx vy vz CD real

    % Define position and velocity vectors
    r_vec = [x; y; z];
    v_vec = [vx; vy; vz];

    % Earth's rotation vector (assumed along z-axis)
    omega_vec = [0; 0; omega_E];

    % Compute velocity relative to atmosphere
    v_rel = v_vec - cross(omega_vec, r_vec);
    v_rel_mag = norm(v_rel);

    % Compute atmospheric density at current altitude
    r = norm(r_vec);
    rho = rho0 * exp(-(r - r0) / H);

    % Drag acceleration using relative velocity
    a_drag = -0.5 * (CD * A / m) * rho * v_rel_mag * v_rel;
    zero_acc_CD = sym(0); % Placeholder for CD derivative
    augmented_acceleration = [a_drag; zero_acc_CD];

    % Convert to a numeric function
    acceleration_func = matlabFunction(augmented_acceleration, ...
        'Vars', {x, y, z, vx, vy, vz, CD});

    % Compute the augmented acceleration and Jacobian of the system
    augmented_acceleration_full = [vx; vy; vz; augmented_acceleration];
    A_sym = jacobian(augmented_acceleration_full, [x, y, z, vx, vy, vz, CD]);

    % Convert to a numeric function
    A_func = matlabFunction(A_sym, 'Vars', {x, y, z, vx, vy, vz, CD});
end

function [acceleration, A] = sc_Dynamics(position_eci, velocity_eci, coeffs_vals, mu_val, ...
    l_max, CD, t, acc_SPH, A_SPH, acc_Drag, A_Drag, params_meas)

    % Convert ECI to ECEF coordinates
    R_eci_to_ecef = DCM_eciToEcef(t);               % Rotation matrix ECI to ECEF
    position_ecef = R_eci_to_ecef * position_eci;   % Transform position to ECEF
 
    % Evaluate acceleration and parameter derivatives at the given inputs
    acceleration_SPH = acc_SPH(position_ecef(1), position_ecef(2), ...
                                     position_ecef(3), mu_val, coeffs_vals');
    acceleration_SPH(1:3) = R_eci_to_ecef' * acceleration_SPH(1:3);
    acceleration_Drag = acc_Drag(position_eci(1), position_eci(2), ...
                                     position_eci(3), velocity_eci(1), velocity_eci(2), ...
                                     velocity_eci(3), CD);

    % Sum acceleration contributions
    acceleration = [acceleration_SPH(1:3) + acceleration_Drag(1:3); ...
        acceleration_SPH(4:end); acceleration_Drag(4:end); ...
        zeros(params_meas, 1)];
    
   if nargout > 1
    % Evaluate Jacobian of zonal harmonics (SPH) at the given position and parameters
    A_SPH_ecef = A_SPH(position_ecef(1), position_ecef(2), position_ecef(3), ...
                    0, 0, 0, mu_val, coeffs_vals');

    % Evaluate Jacobian of drag in ECI
    A_Drag_eci = A_Drag(position_eci(1), position_eci(2), ...
                                 position_eci(3), velocity_eci(1), velocity_eci(2), ...
                                 velocity_eci(3), CD);

    % Transform zonal harmonics Jacobian back to ECI
    Ar_eci = R_eci_to_ecef' * A_SPH_ecef(4:6, 1:3) * R_eci_to_ecef; 
    Amu_eci = R_eci_to_ecef' * A_SPH_ecef(4:6, 7);
    Acl_eci = R_eci_to_ecef' * A_SPH_ecef(4:6, 7:end);

    % Extract drag contributions in ECI
    Ar_drag = A_Drag_eci(4:6, 1:3);  % Partial derivatives wrt position
    Av_drag = A_Drag_eci(4:6, 4:6);  % Partial derivatives wrt velocity
    A_CD = A_Drag_eci(4:6, 7);       % Partial derivative wrt C_D

    % Sum the contributions from both models
    Ar_total = Ar_eci + Ar_drag;
    Av_total = Av_drag; % Only drag affects velocity derivatives

    % Construct the full state transition matrix A 
    A = [zeros(3, 3), eye(3), zeros(3, 1), zeros(3, l_max + 1), zeros(3, params_meas);      % Velocity derivatives
         Ar_total, Av_total, Amu_eci, Acl_eci(:, 2:end), A_CD, zeros(3, params_meas);       % Acceleration derivatives
         zeros(1, 6), 0, zeros(1, l_max + 1), zeros(1, params_meas);                        % Mu derivatives
         zeros(l_max, 7), zeros(l_max, l_max + 1), zeros(l_max, params_meas);                   % Harmonics derivatives
         zeros(1, 7), zeros(1, l_max + 1), zeros(1, params_meas);
         zeros(params_meas, 8+l_max+params_meas)];                          % CD derivatives
   end
end

function dstate = sc_ODE(t, state, coeffs, mu, l_max, CD, acc_SPH, A_SPH, ...
    acc_Drag, A_Drag, params_meas)
    % Extract from state
    position_eci = state(1:3);
    velocity_eci = state(4:6);

    % Compute acceleration and A matrix
    [acceleration_eci, A_eci] = sc_Dynamics(position_eci, velocity_eci, coeffs, mu, ...
         l_max, CD, t, acc_SPH, A_SPH, acc_Drag, A_Drag, params_meas);

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

function [gs_positions_eci, gs_velocities_eci] = compute_gs_eci(gs_location_ecef, times)
    % Constants
    omega_earth = 7.2921158553 * 1e-5; 

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
        gs_velocities_eci(:, i) = omega_earth * cross([0; 0; 1], ...
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
        [gs_pos, gs_vel] = compute_gs_eci(station_ecef0(string(station_id)), ...
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
        d_measurements_dX = [d_rho_dR, d_rho_dV;
                             d_rho_dot_dR, d_rho_dot_dV];

        % Ground station position Jacobian in ECEF
        R_eci_to_ecef = DCM_eciToEcef(trackingData.time(i));
        d_rho_dRs = -d_rho_dR * R_eci_to_ecef';
        d_rho_dot_dRs = -d_rho_dot_dR * R_eci_to_ecef';
        d_measurements_dRs = [d_rho_dRs;
                             d_rho_dot_dRs];

        % Store function handle for measurement model
        measurement_func = @(x) measurement_function(x, gs_pos, gs_vel,trackingData.time(i));

        % Store results in structure
        measurements_struct(i).time = trackingData.time(i);
        measurements_struct(i).stationID = station_id;
        measurements_struct(i).observed = [rho_obs; rho_dot_obs];
        measurements_struct(i).computed = [rho_comp; rho_dot_comp];
        measurements_struct(i).residual = [rho_res; rho_dot_res];
        measurements_struct(i).partials = struct('wrt_X', d_measurements_dX, ...
                                                 'wrt_Rs', d_measurements_dRs);
        measurements_struct(i).covariance = R;
        measurements_struct(i).measurement_function = measurement_func;
    end
end

function [h_x, H, H_dRs] = measurement_function(x, gs_pos, gs_vel, time)
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
    H = [d_rho_dR, d_rho_dV;
         d_rho_dot_dR, d_rho_dot_dV];

    % Position Jacobian Spacecraft
    R_eci_to_ecef = DCM_eciToEcef(time);
    d_rho_dRs = -d_rho_dR * R_eci_to_ecef';
    d_rho_dot_dRs = -d_rho_dot_dR * R_eci_to_ecef';
    H_dRs = [d_rho_dRs; 
        d_rho_dot_dRs];
end

function filter_dataset = create_filter_dataset(measurements_struct, params_dyn, params_meas)
    % Copy the input structure to preserve all fields
    filter_dataset = measurements_struct;

    % Number of measurements
    N = length(measurements_struct);

    % Loop through measurements and modify only the partials field
    for i = 1:N
        % Extract measurement partial derivatives
        d_measurements_dX = measurements_struct(i).partials.wrt_X; % Spacecraft state partials
        d_measurements_dRs = measurements_struct(i).partials.wrt_Rs; % Ground station partials

        % Initialize new partials matrix (2 x num_state_params)
        H = zeros(2, 6 + params_dyn + params_meas); % 2 rows for range & range-rate

        % Fill in the partials w.r.t. spacecraft states (r, v, mu, J1, J2, CD)
        H(:, 1:6) = d_measurements_dX;

        % Define the mapping of station IDs to index shifts
        station_indices = containers.Map({101, 337, 394}, ...
                                         {1, 4, 7}); 
        
        % Identify the ground station and place the corresponding partials
        if isKey(station_indices, measurements_struct(i).stationID)
            idx_shift = station_indices(measurements_struct(i).stationID); % Get relative index
            H(:, 6 + params_dyn + idx_shift : 6 + params_dyn + idx_shift + 2) = d_measurements_dRs;
        end


        % Replace the entire .partials field with H (2 x num_state_params)
        filter_dataset(i).partials = H;
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
        x_ref = trajectory_ref(traj_idx).state(1:n)';

        % Compute STM transition matrix
        STM_dt = STM_t / STM_tm; 

        % Extract measurement data
        prefit_residual = sorted_measurements(t).residual;
        H_tilde = sorted_measurements(t).partials;
        R = sorted_measurements(t).covariance;

        % Prediction Step
        dx_pred = STM_dt * dx;
        P_pred = STM_dt * P * STM_dt'; % OSS:  + 1e-14 * eye(19), Q works!
        % [V, D] = eig(P_pred);
        % D(D < 1e-12) = 1e-12; % Replace negative eigenvalues with a small positive value
        % P_pred = V * D * V';

        % Compute pre-fit residual
        postfit_res = prefit_residual - H_tilde * dx_pred;

        % Kalman Gain computation
        S = H_tilde * P_pred * H_tilde' + R;
        K = P_pred * H_tilde' / S;

        % Update Step
        dx_upd = dx_pred + K * postfit_res;
        P_upd = (I - K * H_tilde) * P_pred * (I - K * H_tilde)' + K * R * K';
        % [V, D] = eig(P_upd);
        % D(D < 1e-12) = 1e-12; % Replace negative eigenvalues with a small positive value
        % P_upd = V * D * V';

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

    % Compute initial state correction
    dx0 = STM_t \ dx_upd;

    fprintf('\nKalman Filter Completed!\n');

    % Store all outputs in a struct (dictionary)
    results = struct(...
        'state_corrected_hist', state_corrected_hist, ...
        'state_deviation_hist', state_deviation_hist, ...
        'P_hist', P_hist, ...
        'postfit_residuals', postfit_residuals, ...
        'state0_deviation', dx0);
end


function compute_rms_errors(results, sorted_measurements)
    % Compute state errors (only at measurement times)
    prefit_res = [sorted_measurements.residual]; % (m x T)
    postfit_res = results.postfit_residuals; % (m x T)

    % Compute RMS for pre-fit residuals
    rms_prefit = sqrt(mean(prefit_res.^2, 2));

    % Compute RMS for post-fit residuals
    rms_postfit = sqrt(mean(postfit_res.^2, 2));

    % Print pre-fit residual RMS errors
    fprintf('\nPre-Fit Residual RMS Errors:\n');
    for i = 1:size(prefit_res, 1)
        fprintf('Residual %d: %.14e\n', i, rms_prefit(i));
    end

    % Print post-fit residual RMS errors
    fprintf('\nPost-Fit Residual RMS Errors:\n');
    for i = 1:size(postfit_res, 1)
        fprintf('Residual %d: %.14e\n', i, rms_postfit(i));
    end
end
