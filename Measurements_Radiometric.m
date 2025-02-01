%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         January 22, 2024
% CLASS:        ASEN 6080: StatOD
% INSTRUCTOR:   Prof. Jay W. Mcmahon
% ASSIGNMENT:   Homework 1 - Exercises: 3, 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read text file - 1 
clear; clc; close;

% File name
filename = 'HW1_prob3b_data.txt';

% Read file
[sc_eci_parsed, gs_eci_parsed, Htilde_parsed] = parse_input_file(filename);

%% Exercise 3 - Compare Measurement Partials wrt S/C State

% Initialize
sc_pos = sc_eci_parsed.r';
sc_vel = sc_eci_parsed.v';
gs_pos = gs_eci_parsed.Rs';
gs_vel = gs_eci_parsed.Vs';

% Compute measurements and partials
[~, partials] = radiometric_measurement(sc_pos, sc_vel, ...
                                                    gs_pos, gs_vel);

% Compute the element-by-element error
errorMatrix = abs(Htilde_parsed - partials.wrt_X) ./ (abs(Htilde_parsed) + 1e-14);

% Plot
gca1 = figure(1);
imagesc(errorMatrix); 
purpleColormap = [linspace(1, 0.5, 256)', linspace(1, 0.3, 256)', ...
    linspace(1, 0.8, 256)'];
colormap(purpleColormap);
colorbar_handle = colorbar; 
colorbar_handle.Label.String = '$|A^\ast_{ij}-A_{ij}|/|A^\ast_{ij}|$ (-)'; 
colorbar_handle.Label.Interpreter = 'latex'; 
colorbar_handle.Label.FontSize = 14; 
xlabel('Column Index (-)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Row Index (-)', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'XTick', 1:size(Htilde_parsed, 2), 'YTick', ...
    1:size(Htilde_parsed, 1), 'TickLabelInterpreter', 'latex');
grid on;
axis square;
exportgraphics(gca1, 'HW1_1b.pdf', 'ContentType','image',...
    'Resolution', 1000);


%% Read text file - 2 
clear; clc; close;

% File name
filename = 'HW1_prob3d_data.txt';

% Read file
[sc_eci_parsed, gs_eci_parsed, Htilde_parsed] = parse_input_file(filename);

%% Exercise 3 - Compare Measurement Partials wrt GS State

% Initialize
sc_pos = sc_eci_parsed.r';
sc_vel = sc_eci_parsed.v';
gs_pos = gs_eci_parsed.Rs';
gs_vel = gs_eci_parsed.Vs';

% Compute measurements and partials
[measurements, partials] = radiometric_measurement(sc_pos, sc_vel, ...
                                                    gs_pos, gs_vel);

% Compute the element-by-element error
errorMatrix = abs(Htilde_parsed - partials.wrt_Rs) ./ (abs(Htilde_parsed) + 1e-14);

% Plot
gca2 = figure(2);
imagesc(errorMatrix); 
purpleColormap = [linspace(1, 0.5, 256)', linspace(1, 0.3, 256)', ...
    linspace(1, 0.8, 256)'];
colormap(purpleColormap);
colorbar_handle = colorbar; 
colorbar_handle.Label.String = '$|A^\ast_{ij}-A_{ij}|/|A^\ast_{ij}|$ (-)'; 
colorbar_handle.Label.Interpreter = 'latex'; 
colorbar_handle.Label.FontSize = 14; 
xlabel('Column Index (-)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Row Index (-)', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'XTick', 1:size(Htilde_parsed, 2), 'YTick', ...
    1:size(Htilde_parsed, 1), 'TickLabelInterpreter', 'latex');
grid on;
axis square;
exportgraphics(gca2, 'HW1_2b.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Exercise 4 - Simulating Measurements
clear; clc; close;
 
% Constants
earth_radius = 6378; % Earth's mean radius in km
omega_earth = 2 * pi / (24 * 3600); % Earth's rotation rate (rad/s)
theta0 = 122; % Initial Earth rotation angle in degrees
elevation_mask = 10; % Elevation mask in degrees
fT_ref = 8.44 * 1e9; % Transmit frequency in Hz (X-Band)
c = 299792.458; % Speed of light in km/s

% Ground station locations (lat, lon in degrees)
stations = [
    -35.398333, 148.981944;  % Station 1
    40.427222, 355.749444;   % Station 2
    35.247164, 243.205       % Station 3
];

% Simulated spacecraft orbit
loaded_data = load('reference_trajectory.mat');

% Extract variables
state_ref = loaded_data.state_ref;
sc_pos = state_ref(:, 1:3)';
sc_vel = state_ref(:, 4:6)';
times = loaded_data.t_ref;

% Pre-allocate measurements
range_measurements = [];
range_rate_measurements = [];
station_indices = [];
measurement_times = [];
is_visibile_GSs = [];
elevation_angles_GSs = [];

% Process each station
for station_idx = 1:size(stations, 1)
    % Compute ground station ECI positions and velocities
    [gs_pos_eci, gs_vel_eci] = compute_gs_eci(stations(station_idx, :), theta0, times);
    
    % Compute visibility mask
    [is_visible, elevation_angles] = compute_visibility_mask(sc_pos, gs_pos_eci, elevation_mask);
    is_visibile_GSs = [is_visibile_GSs; is_visible];
    elevation_angles_GSs = [elevation_angles_GSs; elevation_angles];
    
    % Extract visible measurements
    for i = 1:length(times)
        if is_visible(i)
            % Compute range and range-rate
            [measurement, ~] = radiometric_measurement(sc_pos(:, i), sc_vel(:, i), ...
                                                       gs_pos_eci(:, i), gs_vel_eci(:, i));
            range_measurements = [range_measurements, measurement(1)];
            range_rate_measurements = [range_rate_measurements, measurement(2)];
            station_indices = [station_indices, station_idx];
            measurement_times = [measurement_times, times(i)];
        end
    end
end

% Plot range and range-rate measurements
gca3 = figure(3);
subplot(2, 1, 1);
hold on;
for station_idx = 1:size(stations, 1)
    mask = station_indices == station_idx;
    plot(measurement_times(mask) / 3600, range_measurements(mask), '.', ...
        'DisplayName', sprintf('GS %d', station_idx));
end
xlabel('Time (hours)','Interpreter', 'latex', 'FontSize', 14);
ylabel('$\rho$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
grid on;
subplot(2, 1, 2);
hold on;
for station_idx = 1:size(stations, 1)
    mask = station_indices == station_idx;
    plot(measurement_times(mask) / 3600, range_rate_measurements(mask), '.', ...
        'DisplayName', sprintf('GS %d', station_idx));
end
xlabel('Time (hours)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\dot{\rho}$ (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 8, 'Interpreter', 'latex');
grid on;
exportgraphics(gca3, 'HW1_3b.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Print
fprintf('First measurement time: %.2f seconds\n', min(measurement_times));
fprintf('Last measurement time: %.2f seconds\n', max(measurement_times));

% Plot visibility
gca4 = figure(4);
subplot(2, 1, 1);
hold on;
% Plot elevation angles for each station
for station_idx = 1:size(stations, 1)
    % Plot elevation angles with customized line styles and markers
    plot(times / 3600, elevation_angles_GSs(station_idx, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('GS %d', station_idx));
end
yline(elevation_mask, 'r-.', 'LineWidth', 2, 'DisplayName', 'Min Elevation Angle');
grid on;
xlabel('Time (hours)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$El$ (deg)', 'FontSize', 14, 'Interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 8, 'Interpreter', 'latex');
xlim([0, max(times) / 3600]);
subplot(2, 1, 2);
hold on;
for station_idx = 1:size(stations, 1)
    % Plot visibility as a step function
    stairs(times / 3600, is_visibile_GSs(station_idx, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('GS %d', station_idx));
end
xlabel('Time (hours)', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Visibility (Logical)', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
legend('Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
xlim([0, max(times) / 3600]);
ylim([-0.1, 1.1]);
exportgraphics(gca4, 'HW1_4b.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Convert range to Range Units (RU) and range rate to Doppler shift (Hz)
range_ru = (221 / 749 * range_measurements) / (c / fT_ref);
doppler_shift = -2 * range_rate_measurements * fT_ref / c;

% Plot range in RU and Doppler shift
gca5 = figure(5);
subplot(2, 1, 1);
hold on;
for station_idx = 1:size(stations, 1)
    mask = station_indices == station_idx;
    plot(measurement_times(mask) / 3600, range_ru(mask), '.', ...
        'DisplayName', sprintf('GS %d', station_idx));
end
xlabel('Time (hours)', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('RU $\cdot 10^{7} \, (s^2)$', 'FontSize', 12, 'Interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
grid on;
subplot(2, 1, 2);
hold on;
for station_idx = 1:size(stations, 1)
    mask = station_indices == station_idx;
    plot(measurement_times(mask) / 3600, doppler_shift(mask), '.', ...
        'DisplayName', sprintf('GS %d', station_idx));
end
xlabel('Time (hours)', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('$\Delta f \cdot 10^{5}$ (Hz)', 'FontSize', 12, 'Interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
grid on;
exportgraphics(gca5, 'HW1_5b.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Add Gaussian noise to range-rate measurements
sigma = 0.5 / 1e6; % 0.5 mm/s converted to km/s
noisy_range_rate = range_rate_measurements + ...
    sigma * randn(size(range_rate_measurements));

% Plot noisy measurements and differences
gca6 = figure(6);
subplot(2, 1, 1);
hold on;
for station_idx = 1:size(stations, 1)
    mask = station_indices == station_idx;
    plot(measurement_times(mask) / 3600, noisy_range_rate(mask), '.', ...
        'DisplayName', sprintf('GS %d', station_idx));
end
xlabel('Time (hours)', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Noisy $\dot{\rho}$ (km/s)', 'FontSize', 12, 'Interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
grid on;
subplot(2, 1, 2);
hold on;
for station_idx = 1:size(stations, 1)
    mask = station_indices == station_idx;
    plot(measurement_times(mask) / 3600, ...
        noisy_range_rate(mask) - range_rate_measurements(mask), '.', ...
        'DisplayName', sprintf('GS %d', station_idx));
end
yline(3 * sigma, 'r-.', 'LineWidth', 2, 'DisplayName', '$\pm 3 \sigma$');
yline(-3 * sigma, 'r-.', 'LineWidth', 2, 'HandleVisibility', 'off');
xlabel('Time (hours)', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Residuals $\cdot 10^{-6}$ (km/s)', 'FontSize', 12, 'Interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
grid on;
exportgraphics(gca6, 'HW1_6b.pdf', 'ContentType','image',...
    'Resolution', 1000);


%% Helper Functions

function [measurements, partials] = radiometric_measurement(sc_pos, sc_vel, gs_pos, gs_vel)
    % Relative position and velocity
    dR = sc_pos - gs_pos; 
    dV = sc_vel - gs_vel; 

    % Compute range and range-rate
    rho = norm(dR);
    rho_dot = dot(dR, dV) / rho;

    % Compute partial derivatives
    % Partial derivatives w.r.t. spacecraft state [r; v]
    d_rho_dR = dR' / rho; 
    d_rho_dV = zeros(1, 3); 
    d_rho_dot_dR = (rho * dV' - rho_dot * dR') / rho^2;
    d_rho_dot_dV = dR' / rho;

    % Combine for spacecraft state
    d_measurements_dX = [d_rho_dR, d_rho_dV; 
                         d_rho_dot_dR, d_rho_dot_dV]; 

    % Partial derivatives w.r.t. GS state [Rs; Vs]
    d_rho_dRs = -d_rho_dR;
    % d_rho_dVs = zeros(1, 3); 
    d_rho_dot_dRs = -d_rho_dot_dR; 
    % d_rho_dot_dVs = -d_rho_dot_dV; 

    % Combine for GS position
    d_measurements_dRs = [d_rho_dRs; 
                         d_rho_dot_dRs]; 

    % Prepare outputs
    measurements = [rho; rho_dot];
    partials = struct('wrt_X', d_measurements_dX, ...
                      'wrt_Rs', d_measurements_dRs);
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


function [state, station_info, Htilde] = parse_input_file(filename)
    % Open the file
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Could not open file %s.', filename);
    end

    % Initialize outputs
    state = struct('r', [], 'v', []);
    station_info = struct('Rs', [], 'Vs', []);
    Htilde = [];

    % Read the file line by line
    while ~feof(fileID)
        line = strtrim(fgets(fileID));
        
        % Parse spacecraft position (r)
        if startsWith(line, 'r(')
            parts = sscanf(line, 'r(%d) = %f');
            idx = parts(1);
            state.r(idx) = parts(2);
        end

        % Parse spacecraft velocity (v)
        if startsWith(line, 'v(')
            parts = sscanf(line, 'v(%d) = %f');
            idx = parts(1);
            state.v(idx) = parts(2);
        end

        % Parse station position (Rs)
        if startsWith(line, 'Rs(')
            parts = sscanf(line, 'Rs(%d) = %f');
            idx = parts(1);
            station_info.Rs(idx) = parts(2);
        end

        % Parse station velocity (Vs)
        if startsWith(line, 'Vs(')
            parts = sscanf(line, 'Vs(%d) = %f');
            idx = parts(1);
            station_info.Vs(idx) = parts(2);
        end

        % Parse Htilde matrix
        if startsWith(line, 'Htilde(')
            parts = sscanf(line, 'Htilde(%d,%d) = %f');
            row = parts(1);
            col = parts(2);
            value = parts(3);

            % Dynamically resize Htilde if needed
            if size(Htilde, 1) < row || size(Htilde, 2) < col
                Htilde(max(row, size(Htilde, 1)), max(col, size(Htilde, 2))) = 0;
            end
            Htilde(row, col) = value;
        end
    end

    % Close the file
    fclose(fileID);
end
