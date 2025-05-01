% B-PLANE OVERLAY PLOT FOR MULTIPLE OD SOLUTIONS
clear; clc; close all;

% USER SETTINGS
excluded_solutions = {'OD05'};

% Constants
R_earth = 6378.1363;  % Earth's mean radius [km]

% Load Files
files = dir('bplane_*.mat');
gca123 = figure('Position', [100, 100, 1400, 900]); hold on;
colors = lines(length(files));

% Initialize for legend
legend_handles = gobjects(0);
legend_labels = {};

% Loop Over Solutions
plot_idx = 0;
for k = 1:length(files)
    % Load file
    data = load(files(k).name);
    Rstruct = data.Rstruct;

    % Extract solution name
    tokens = regexp(files(k).name, 'bplane_(.+)\.mat', 'tokens');
    solution_label = tokens{1}{1};

    % Check if excluded
    if ismember(solution_label, excluded_solutions)
        continue;
    end

    plot_idx = plot_idx + 1;

    % Extract B-plane parameters
    BdotT = Rstruct.BdotT;
    BdotR = Rstruct.BdotR;
    P_r_B = Rstruct.P_r_B;

    % Plot mean point and store handle
    h = plot(BdotT, BdotR, 'x', 'Color', colors(plot_idx,:), ...
        'MarkerSize', 10, 'LineWidth', 2);
    legend_handles(end+1) = h;

    % Compute standard deviations
    sigma_BdotT = sqrt(P_r_B(1,1));
    sigma_BdotR = sqrt(P_r_B(2,2));
    
    % Construct detailed legend label
    legend_labels{end+1} = sprintf(...
    '$\\texttt{%s:}\\ B \\cdot \\hat{T} = %.1f \\pm %.1f\\ \\mathrm{km},\\ B \\cdot \\hat{R} = %.1f \\pm %.1f\\ \\mathrm{km}$', ...
    solution_label, BdotT, sigma_BdotT, BdotR, sigma_BdotR);

    % Plot covariance ellipse
    error_ellipse(P_r_B(1:2,1:2), [BdotT; BdotR], ...
        'C', colors(plot_idx,:));

    % Annotate with jitter
    text(BdotT + 500 * randn(), BdotR + 500 * randn(), solution_label, ...
        'FontSize', 12, 'FontWeight', 'bold', ...
        'Color', colors(plot_idx,:));

    % Print to terminal
    rho_TR = P_r_B(1,2) / (sigma_BdotT * sigma_BdotR);
    fprintf('\n===== %s =====\n', solution_label);
    fprintf('B · T       = %.14f km\n', BdotT);
    fprintf('B · R       = %.14f km\n', BdotR);
    fprintf('σ(B · T)    = %.14f km\n', sigma_BdotT);
    fprintf('σ(B · R)    = %.14f km\n', sigma_BdotR);
    fprintf('ρ(BT, BR)   = %.14f\n', rho_TR);
    fprintf('===========================\n');
end

% Plot Earth Disk 
theta = linspace(0, 2*pi, 300);
x_earth = R_earth * cos(theta);
y_earth = R_earth * sin(theta);
plot(x_earth, y_earth, 'k-', 'LineWidth', 2);
text(0, 0, 'Earth', 'HorizontalAlignment', 'center', ...
    'FontSize', 12, 'FontWeight', 'bold');

% Final Plot Formatting
xlabel('$B \cdot \hat{T}$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$B \cdot \hat{R}$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
grid on; box on; axis equal;
set(gca, 'YDir', 'reverse');
legend(legend_handles, legend_labels, ...
    'Location', 'southwest', ...
    'Interpreter', 'latex', ...
    'FontSize', 16);
exportgraphics(gca123, 'Bplane_2.pdf', 'ContentType','vector');


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

