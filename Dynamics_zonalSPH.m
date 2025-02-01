%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         January 22, 2024
% CLASS:        ASEN 6080: StatOD
% INSTRUCTOR:   Prof. Jay W. Mcmahon
% ASSIGNMENT:   Homework 1 - Exercises: 1, 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read text file - 1 
clear; clc; close;

% File name
filename = 'HW1_prob1c_data.txt';

% Open the file
fileID = fopen(filename, 'r');

% Initialize variables
r_parsed = zeros(3, 1);
v_parsed = zeros(3, 1);
A_parsed = zeros(9, 9);

% Read the file line by line
while ~feof(fileID)
    line = fgetl(fileID);
    if contains(line, 'r(')
        % Parse r components
        idx = sscanf(line, 'r(%d) = %f');
        r_parsed(idx(1)) = idx(2);
    elseif contains(line, 'v(')
        % Parse v components
        idx = sscanf(line, 'v(%d) = %f');
        v_parsed(idx(1)) = idx(2);
    elseif contains(line, 'mu =')
        % Parse mu
        mu_parsed = sscanf(line, 'mu = %f');
    elseif contains(line, 'J2 =')
        % Parse J2
        J2_parsed = sscanf(line, 'J2 = %f');
    elseif contains(line, 'J3 =')
        % Parse J3
        J3_parsed = sscanf(line, 'J3 = %f');
    elseif contains(line, 'A(')
        % Parse A matrix components
        idx = sscanf(line, 'A(%d,%d) = %f');
        A_parsed(idx(1), idx(2)) = idx(3);
    end
end

% Close the file
fclose(fileID);

%% Exercise 1 - Compare SPH Partials

% Define initial conditions
position_eci = r_parsed; % Initial position 
velocity_eci = v_parsed; % Initial velocity 
coeffs = [0, J2_parsed, J3_parsed]; % Zonal harmonic coefficients (i.e., [J1, J2, J3])
R = 6378; % Earth radius 
mu = mu_parsed; % Gravitational parameter 
l_max = 3; % Zonal Harmonics Terms 
t = 0; % Initial Time

% Dynamics function
[acc_func, A_func] = precomputeZonalSPH(R, l_max);
[acc, A] = zonalSPH_Dynamics(position_eci, coeffs, mu, l_max, t, acc_func, A_func);
A(8, :) = []; A(:, 8) = []; 

% Compute the element-by-element error
errorMatrix = abs(A_parsed - A) ./ (abs(A_parsed) + 1e-14);

% Plot results
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
set(gca, 'XTick', 1:size(A, 2), 'YTick', 1:size(A, 1), ...
    'TickLabelInterpreter', 'latex');
grid on;
axis square;
exportgraphics(gca1, 'HW1_1.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Read text file - 2

% File name
filename = 'HW1_prob2b_data.txt';

% Open the file
fileID = fopen(filename, 'r');

% Initialize variables
X0_parsed = [];
Phi0_parsed = [];
Xdot_parsed = [];
Phidot_parsed = [];

% Read line by line
while ~feof(fileID)
    line = strtrim(fgets(fileID));

    % Parse X0
    if startsWith(line, 'X0(')
        % Extract index and value
        parts = sscanf(line, 'X0(%d) = %f');
        X0_parsed(parts(1)) = parts(2);
    end

    % Parse Phi0
    if startsWith(line, 'Phi0(')
        % Extract row, column, and value
        parts = sscanf(line, 'Phi0(%d,%d) = %f');
        Phi0_parsed(parts(1), parts(2)) = parts(3);
    end

    % Parse Xdot
    if startsWith(line, 'Xdot(')
        % Extract index and value
        parts = sscanf(line, 'Xdot(%d) = %f');
        Xdot_parsed(parts(1)) = parts(2);
    end

    % Parse Phidot
    if startsWith(line, 'Phidot(')
        % Extract row, column, and value
        parts = sscanf(line, 'Phidot(%d,%d) = %f');
        Phidot_parsed(parts(1), parts(2)) = parts(3);
    end
end

% Re-organize
X0_parsed = X0_parsed';
Xdot_parsed = Xdot_parsed';

% Close the file
fclose(fileID);

%% Exercise 2a - Compare ODE

% Define initial conditions
coeffs = [0, X0_parsed(end)]; % Zonal harmonic coefficients (i.e., [J1, J2])
R = 6378; % Earth radius 
mu = 398600.4418; % Gravitational parameter 
l_max = 2; % Zonal Harmonics Terms 
t = 0; % Initial Time

% Modify parsed Phi0 (adding mu and J0)
Phi0_parsed_mod = zeros(9);
Phi0_parsed_mod(1:6, 1:6) = Phi0_parsed(1:6, 1:6);
Phi0_parsed_mod(7, 7) = 1;
Phi0_parsed_mod(8, 8) = 1;
Phi0_parsed_mod(1:6, 9) = Phi0_parsed(1:6, 7); 
Phi0_parsed_mod(9, 1:6) = Phi0_parsed(7, 1:6); 
Phi0_parsed_mod(9, 9) = Phi0_parsed(7, 7);    

% Modified state parsed
X0_parsed_mod = [X0_parsed(1:end-1); mu; 0; X0_parsed(end); ...
    reshape(Phi0_parsed_mod, [], 1)];

% Initialize dynamics
[acc_func, A_func] = precomputeZonalSPH(R, l_max);

% Compute ODE
dstate = zonalSPH_ODE(t, X0_parsed_mod, coeffs, mu, l_max, acc_func, A_func); 

% Extract Xdot
Xdot = [dstate(1:6); dstate(9)];

% Extract Phidot
Phidot = reshape(dstate(10:end), 9, 9);
Phidot([7, 8], :) = []; % Remove rows
Phidot(:, [7, 8]) = []; % Remove columns

% Compute relative differences
relative_diff_Xdot = abs(Xdot_parsed - Xdot) ./ abs(Xdot_parsed);         % Compare Xdot
relative_diff_Xdot(isnan(relative_diff_Xdot)) = 0;
relative_diff_Phidot = abs(Phidot_parsed - Phidot) ./ abs(Phidot_parsed); % Compare Phidot
relative_diff_Phidot(isnan(relative_diff_Phidot)) = 0;

% Plots
gca2 = figure(2);
plot(1:length(relative_diff_Xdot), relative_diff_Xdot, 'b*','LineWidth', 2);
grid on;
xlabel('Index (-)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$|a^\ast_{i} - a_{i}|/|a^\ast_{i}|\cdot 10^{-10}$ (-)', 'Interpreter', 'latex', ...
    'FontSize', 14);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
exportgraphics(gca2, 'HW1_2.pdf', 'ContentType','image',...
    'Resolution', 1000);

gca3 = figure(3);
imagesc(relative_diff_Phidot);
purpleColormap = [linspace(1, 0.5, 256)', linspace(1, 0.3, 256)', ...
    linspace(1, 0.8, 256)'];
colormap(purpleColormap);
colorbar_handle = colorbar;
colorbar_handle.Label.String = '$|A^\ast_{ij}-A_{ij}|/|A^\ast_{ij}|$ (-)'; 
colorbar_handle.Label.Interpreter = 'latex';
colorbar_handle.Label.FontSize = 14;
xlabel('Column Index (-)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Row Index (-)', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
axis square;
grid on;
exportgraphics(gca3, 'HW1_3.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Exercise 2b - Compare STM Deviations

% Constants
coeffs = [0, 1.0826269 * 1e-3]; % Zonal harmonic coefficients (i.e., [J1, J2])
R = 6378; % Earth radius 
mu = 398600.4415;  % Gravitational parameter 
l_max = 2; % Zonal Harmonics Terms 
t = 0; % Initial Time

% Orbital parameters (GEO)
% a = 42163;              % Semi-major axis (km)
% e = 1e-5;            % Eccentricity
% i = deg2rad(i);      % Inclination (rad)
% RAAN = deg2rad(70);   % Right Ascension of Ascending Node (rad)
% omega = deg2rad(0);  % Argument of Perigee (rad)
% nu0 = deg2rad(0);     % True Anomaly (rad)

% Orbital parameters (Hyperbolic)
% state0 =[-7737.559071593195; -43881.87809094457; 0.0; 3.347424567061589; ...
%         3.828541915617483; 0.0];

% Orbital parameters 
a = 1e4;              % Semi-major axis (km)
e = 0.001;            % Eccentricity (-)
i = deg2rad(40);      % Inclination (rad)
RAAN = deg2rad(80);   % Right Ascension of Ascending Node (rad)
omega = deg2rad(40);  % Argument of Perigee (rad)
nu0 = deg2rad(0);     % True Anomaly (rad)

% Initial state (reference orbit)
state0 = orbitalElementsToCartesian(mu, a, e, i, RAAN, omega, nu0);
state0 = [state0; mu; 0; coeffs(end); reshape(eye(9), [], 1)]; % Adding mu, J1, J2

% Perturbation vector
delta_x0 = [1; 0; 0; 0; 1e-2; 0; 0; 0; 0]; % Initial perturbation 
state_perturbed0 = state0 + [delta_x0; zeros(81, 1)];

% Initialize ODEs
[acc_func, A_func] = precomputeZonalSPH(R, l_max);

% Simulation parameters
T_orbit = 2 * pi * sqrt(a^3 / mu); % Orbital period (s)
t_span = 0:10:(15 * T_orbit);        % Simulate for 15 orbits, every 10 seconds
options = odeset('RelTol', 2.22045 * 1e-14, 'AbsTol',  2.22045 * 1e-14);

% Integrate reference trajectory
[t_ref, state_ref] = ode45(@(t, state) ...
     zonalSPH_ODE(t, state, coeffs, mu, l_max, acc_func, A_func), ...
     t_span, state0, options);
save('reference_trajectory.mat', 't_ref', 'state_ref');

% Integrate perturbed trajectory with STM
[t_perturbed, state_perturbed] = ode45(@(t, state) ...
     zonalSPH_ODE(t, state, coeffs,  mu, l_max, acc_func, A_func), ...
     t_span, state_perturbed0, options);

% Extract true deviation
delta_x_true = state_perturbed(:, 1:6) - state_ref(:, 1:6);

% STM-based deviation propagation
STM = reshape(state_perturbed(:, 10:end), [], 9, 9);
STM = STM(: , 1:6, 1:6);
delta_x_stm = zeros(size(delta_x_true));
for k = 1:length(t_perturbed)
    delta_x_stm(k, :) = (squeeze(STM(k, :, :)) * delta_x0(1:6))';
end

% Plot
gca4 = figure(4);
subplot(2,1,1);
plot(t_perturbed / 3600, delta_x_true(:, 1:3), 'LineWidth', 2); hold on;
plot(t_perturbed / 3600, delta_x_stm(:, 1:3), '--', 'LineWidth', 2);
xlabel('Time (hours)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\delta \mathbf{r}$ (km)', 'FontSize', 14, 'Interpreter', 'latex');
legend({'$\delta x_{True}$', '$\delta y_{True}$', '$\delta z_{True}$', ...
        '$\delta x_{STM}$', '$\delta y_{STM}$', '$\delta z_{STM}$'}, ...
        'Location', 'northwest', 'FontSize', 8, 'Interpreter', 'latex');
grid on; 
set(gca, 'FontSize', 12, 'LineWidth', 1); 
xlim([min(t_perturbed / 3600), max(t_perturbed / 3600)]);
subplot(2,1,2);
plot(t_perturbed / 3600, delta_x_true(:, 4:6), 'LineWidth', 2); hold on;
plot(t_perturbed / 3600, delta_x_stm(:, 4:6), '--', 'LineWidth', 2);
xlabel('Time (hours)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\delta \mathbf{v}$ (km/s)', 'FontSize', 14, 'Interpreter', 'latex');
legend({'$\delta v_{x,True}$', '$\delta v_{y,True}$', '$\delta v_{z,True}$', ...
        '$\delta v_{x,STM}$', '$\delta v_{y,STM}$', '$\delta v_{z,STM}$'}, ...
        'Location', 'northwest', 'FontSize', 8, 'Interpreter', 'latex');
grid on; 
set(gca, 'FontSize', 12, 'LineWidth', 1); 
xlim([min(t_perturbed / 3600), max(t_perturbed / 3600)]);
exportgraphics(gca4, 'HW1_4.pdf', 'ContentType','image',...
    'Resolution', 1000);

gca5 = figure(5);
subplot(2,1,1);
plot(t_perturbed / 3600, delta_x_true(:, 1) - delta_x_stm(:, 1), ...
    'LineWidth', 2, 'DisplayName', 'Position $\delta x$'); hold on;
plot(t_perturbed / 3600, delta_x_true(:, 2) - delta_x_stm(:, 2), ...
    'LineWidth', 2, 'DisplayName', 'Position $\delta y$');
plot(t_perturbed / 3600, delta_x_true(:, 3) - delta_x_stm(:, 3), ...
    'LineWidth', 2, 'DisplayName', 'Position $\delta z$');
xlabel('Time (hours)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\delta \mathbf{r}_{true} - \delta \mathbf{r}_{STM}$(km)', 'FontSize', ...
    14, 'Interpreter', 'latex');
legend('Location', 'northwest', 'FontSize', 8, 'Interpreter', 'latex');
grid on; 
set(gca, 'FontSize', 12, 'LineWidth', 1); 
xlim([min(t_perturbed / 3600), max(t_perturbed / 3600)]); 
subplot(2,1,2);
plot(t_perturbed / 3600, delta_x_true(:, 4) - delta_x_stm(:, 4), ...
    'LineWidth', 2, 'DisplayName', 'Velocity $\delta v_x$'); hold on;
plot(t_perturbed / 3600, delta_x_true(:, 5) - delta_x_stm(:, 5), ...
    'LineWidth', 2, 'DisplayName', 'Velocity $\delta v_y$');
plot(t_perturbed / 3600, delta_x_true(:, 6) - delta_x_stm(:, 6), ...
    'LineWidth', 2, 'DisplayName', 'Velocity $\delta v_z$');
xlabel('Time (hours)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\delta \mathbf{v}_{true} - \delta \mathbf{v}_{STM}$(km/s)', 'FontSize', ...
    14, 'Interpreter', 'latex');
legend('Location', 'northwest', 'FontSize', 8, 'Interpreter', 'latex');
grid on; 
set(gca, 'FontSize', 12, 'LineWidth', 1); 
xlim([min(t_perturbed / 3600), max(t_perturbed / 3600)]); 
exportgraphics(gca5, 'HW1_5.pdf', 'ContentType','image',...
    'Resolution', 1000);

gca6 = figure(6);
plot3(state_ref(:, 1), state_ref(:, 2), state_ref(:, 3), 'b', 'LineWidth', 2);  
hold on;  
plot3(state_perturbed(:, 1), state_perturbed(:, 2), state_perturbed(:, 3), ...
    'r', 'LineWidth', 2);   
plot_earth();  
axis equal;
xlabel('x (km)', 'Interpreter', 'latex');
ylabel('y (km)', 'Interpreter', 'latex');
zlabel('z (km)', 'Interpreter', 'latex');
legend('Trajectory True', 'Trajectory Perturbed', 'Earth', 'Location', 'northwest',...
    'Interpreter', 'latex');
grid on;
exportgraphics(gca6, 'HW1_6.pdf', 'ContentType','image',...
    'Resolution', 1000);


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
