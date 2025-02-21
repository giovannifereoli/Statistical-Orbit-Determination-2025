%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         January 22, 2024
% CLASS:        ASEN 6080: StatOD
% INSTRUCTOR:   Prof. Jay W. Mcmahon
% ASSIGNMENT:   Project 1 - Comparisons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close;

%% Parse Data

% Define folder paths
load("Batch/results_Batch_Iter3.mat");
load("LKF/results_LKF_Iter1.mat");
load("LKF/meas_times.mat");

%% Compare Batch and CKF Covariance Ellipsoids in 3D

P_final_batch = results_Batch_Iter3.P_hist(1:6, 1:6, end);
P_final_LKF = results_LKF_Iter1.P_hist(1:6, 1:6, end);

state_final_batch = results_Batch_Iter3.state_corrected_hist(:, end);
state_final_LKF = results_LKF_Iter1.state_corrected_hist(:, end);
state_corrected_history_batch = results_Batch_Iter3.state_corrected_hist';
state_corrected_history_LKF = results_LKF_Iter1.state_corrected_hist';

% Covariance Ellipsoids for Position and Velocity
gca11 = figure(11);
set(gca11, 'Position', [100, 100, 1400, 600]); 
subplot(1,2,1);
hold on;
plot_covariance_ellipsoid(P_final_batch(1:3, 1:3), state_final_batch(1:3), 'b', 0.1, 'Batch', 'b');
plot_covariance_ellipsoid(P_final_LKF(1:3, 1:3), state_final_batch(1:3), 'r', 0.1, 'LKF', 'r');
xlabel('$x$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('$z$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
title('Covariance Ellipsoid $\mathbf{P}_{\mathbf{r}}$', 'Interpreter', 'latex', 'FontSize', 16);
legend('show', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
grid on; axis equal;
hold off;
view(3);
subplot(1,2,2);
hold on;
plot_covariance_ellipsoid(P_final_batch(4:6, 4:6), state_final_batch(4:6), 'b', 0.1, 'Batch', 'b');
plot_covariance_ellipsoid(P_final_LKF(4:6, 4:6), state_final_batch(4:6), 'r', 0.1, 'LKF', 'r');
xlabel('$v_x$ (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$v_y$ (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('$v_z$ (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
title('Covariance Ellipsoid $\mathbf{P}_{\mathbf{v}}$', 'Interpreter', 'latex', 'FontSize', 16);
legend('show', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
grid on; axis equal;
hold off;
view(3);
exportgraphics(gca11, 'Covariance_Ellipsoids.pdf', 'ContentType', 'vector', 'Resolution', 1000);


%% State Corrected History with Uncertainty Bounds

% Plot
colors = {'b', 'r', 'g', 'm', 'c', [1 0.5 0]}; 
gca15 = figure(15);
set(gca15, 'Position', [100, 100, 1400, 600]); 
for i = 1:3
    subplot(3,2,2*i-1); hold on;
    semilogy(time, abs(state_corrected_history_batch(:, i+1) - state_corrected_history_LKF(:, i+1)), ...
             'o', 'MarkerSize', 4, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k');
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(['$|\delta ' char(88+i-1) '|$ (km)'], 'Interpreter', 'latex', 'FontSize', 14);
    grid on; set(gca, 'YScale', 'log'); % Log scale
end
for i = 4:6
    subplot(3,2,2*i-6); hold on;
    semilogy(time, abs(state_corrected_history_batch(:, i+1) - state_corrected_history_LKF(:, i+1)), ...
             'o', 'MarkerSize', 4, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k');
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(['$|\delta v_' char(120+i-4) '|$ (km/s)'], 'Interpreter', 'latex', 'FontSize', 14);
    grid on; set(gca, 'YScale', 'log'); % Log scale
end
exportgraphics(gca15, 'State_Corrected_History_Difference_Log.pdf', ...
    'ContentType', 'vector', 'Resolution', 1000);

%% Function to Plot 3D Covariance Ellipsoids

function plot_covariance_ellipsoid(P, mu, color, alpha, label, edgeC)
    [V, D] = eig(P);
    [X, Y, Z] = sphere(30);
    unit_sphere = [X(:), Y(:), Z(:)]';
    ellipsoid_points = V * sqrt(D) * unit_sphere;
    
    X_ellip = reshape(ellipsoid_points(1, :) + mu(1), size(X));
    Y_ellip = reshape(ellipsoid_points(2, :) + mu(2), size(Y));
    Z_ellip = reshape(ellipsoid_points(3, :) + mu(3), size(Z));
    
    surf(X_ellip, Y_ellip, Z_ellip, 'FaceColor', color, 'EdgeColor', ...
        edgeC, 'FaceAlpha', alpha, 'DisplayName', label);
    camlight;
    lighting phong;
end
