function [Tref, Xref, Phi_xx, Phi_xp, Phi_xc] = integrate_reference(x0in, c0, p_all, N_batches, ts_batches, tobs, dyn_function_handle, Nx, Nc, Np, tau_p, integ_args)
    % Integrate reference trajectory with stochastics in a single pass
    % Inputs:
    %   x0in: Initial state vector (Nx x 1)
    %   c0: Initial parameter vector (Nc x 1)
    %   p_all: Stochastic parameters per batch (Np x N_batches)
    %   N_batches: Number of batches
    %   ts_batches: Batch start times (1 x N_batches)
    %   tobs: Observation times (1 x Nobs)
    %   dyn_function_handle: Dynamics function handle
    %   Nx, Nc, Np: Dimensions of state, parameters, stochastics
    %   tau_p: Stochastic time constant
    %   integ_args: Additional arguments for integrator (e.g., {REarth})
    % Outputs:
    %   Tref: Times (1 x Ntimes)
    %   Xref: States (Nx x Ntimes)
    %   Phi_xx: State transition matrices (Nx x Nx x Ntimes)
    %   Phi_xp: State-to-stochastic partials (Nx x (Np*N_batches) x Ntimes)
    %   Phi_xc: State-to-parameter partials (Nx x Nc x Ntimes)

    % Collect all unique times for integration
    t0 = 0;
    tend = tobs(end);
    t_all = unique(sort([t0 ts_batches tobs tend]));

    % Flatten stochastic parameters
    p_total = reshape(p_all, [], 1); % Np*N_batches x 1

    % Build integrator initial state
    xint0 = [x0in; reshape(eye(Nx), [], 1); zeros(Nx * (Np * N_batches), 1); zeros(Nx * Nc, 1)];

    % Prepare arguments for ODE
    all_args = {c0, p_total, ts_batches, integ_args};

    % Set ODE options
    abs_tol = 1e-10;
    ode_options = odeset('RelTol', 1e-10, 'AbsTol', abs_tol);

    % Integrate
    [Tref, Xint] = ode45(dyn_function_handle, t_all, xint0, ode_options, all_args{:});

    % Post-process
    Tref = Tref';
    Xref = Xint(:, 1:Nx)';
    Phi_xx = reshape(Xint(:, Nx+1:Nx+Nx^2)', [Nx Nx length(t_all)]);
    Phi_xp = reshape(Xint(:, Nx+Nx^2+1:Nx+Nx^2 + Nx*(Np*N_batches))', [Nx (Np*N_batches) length(t_all)]);
    Phi_xc = reshape(Xint(:, Nx+Nx^2 + Nx*(Np*N_batches)+1:end)', [Nx Nc length(t_all)]);
end

function xdot = two_body_wPhi_stochastics(t, X, c, p, ts_batches, integ_args)
    % Dynamics function for two-body problem with J2 and batch stochastics
    % Inputs:
    %   t: Current time
    %   X: Integration state vector [x; reshape(Phi_xx, Nx*Nx, 1); reshape(Phi_xp, Nx*(Np*N_batches), 1); reshape(Phi_xc, Nx*Nc, 1)]
    %   c: Parameter vector (Nc x 1, e.g., [MU; J2])
    %   p: Flattened stochastic vector (Np*N_batches x 1)
    %   ts_batches: Batch start times (1 x N_batches)
    %   integ_args: Additional arguments {REarth, Nx, Nc, Np, tau_p}
    % Outputs:
    %   xdot: Derivatives [xdot; reshape(Phi_xx_dot, Nx*Nx, 1); reshape(Phi_xp_dot, Nx*(Np*N_batches), 1); reshape(Phi_xc_dot, Nx*Nc, 1)]

    % Extract arguments
    REarth = integ_args{1};
    Nx = integ_args{2};
    Nc = integ_args{3};
    Np = integ_args{4};
    tau_p = integ_args{5};

    % Extract state and partials
    x = X(1:Nx); % Position and velocity
    phi_xx = reshape(X(Nx+1:Nx+Nx^2), Nx, Nx);
    phi_xp = reshape(X(Nx+Nx^2+1:Nx+Nx^2 + Nx*(Np*N_batches)), Nx, Np*N_batches);
    phi_xc = reshape(X(Nx+Nx^2 + Nx*(Np*N_batches)+1:end), Nx, Nc);

    % Extract position and velocity
    pos = x(1:3); % [x y z]
    vel = x(4:6); % [vx vy vz]

    % Extract parameters (assume c(1) = MU, c(2) = J2 if Nc >= 2)
    MU = c(1);
    J2 = c(2) if Nc >= 2 else 0;

    % Compute radius terms
    radius_2 = sum(pos.^2);
    radius = sqrt(radius_2);
    radius_3 = radius_2 * radius;
    radius_5 = radius_3 * radius_2;
    radius_7 = radius_5 * radius_2;
    radius_9 = radius_7 * radius_2;

    % Gravitational constant for J2
    GravCoeff = J2 * MU * REarth^2 if Nc >= 2 else 0;

    % Two-body + J2 acceleration
    accel = -MU * pos / radius_3;
    if Nc >= 2 && J2 ~= 0
        accel(1) = accel(1) + (15 * GravCoeff * pos(3)^2 * pos(1) / radius_7 - 3 * GravCoeff * pos(1) / radius_5) / 2;
        accel(2) = accel(2) + (15 * GravCoeff * pos(3)^2 * pos(2) / radius_7 - 3 * GravCoeff * pos(2) / radius_5) / 2;
        accel(3) = accel(3) + (15 * GravCoeff * pos(3)^3 / radius_7 - 9 * GravCoeff * pos(3) / radius_5) / 2;
    end

    % Add stochastic accelerations for current batch
    p_all = reshape(p, Np, []);
    accel_stoch = zeros(3,1);
    fp = zeros(Nx, Np*N_batches);
    if t < ts_batches(1)
        % Before first batch, no stochastics
    else
        k = find(ts_batches <= t, 1, 'last');
        if ~isempty(k) && (k == length(ts_batches) || t < ts_batches(k+1))
            Phi_pp_k = exp(-(t - ts_batches(k)) / tau_p) * eye(Np);
            accel_stoch = Phi_pp_k * p_all(:,k);
            col_start = (k-1)*Np + 1;
            col_end = k*Np;
            fp(4:6, col_start:col_end) = Phi_pp_k;
        end
    end
    accel = accel + accel_stoch;

    % State derivatives
    xdot = [vel; accel];

    % State Jacobian (A = df/dx)
    A = zeros(Nx, Nx);
    A(1:3, 4:6) = eye(3);
    A(4,1) = -MU/radius_3 + 3*MU*pos(1)^2/radius_5 + 15/2*GravCoeff*pos(3)^2/radius_7 - 105/2*GravCoeff*pos(3)^2*pos(1)^2/radius_9 - 3/2*GravCoeff/radius_5 + 15/2*GravCoeff*pos(1)^2/radius_7;
    A(4,2) = 3*MU*pos(1)*pos(2)/radius_5 - 105/2*GravCoeff*pos(3)^2*pos(1)*pos(2)/radius_9 + 15/2*GravCoeff*pos(1)*pos(2)/radius_7;
    A(4,3) = 3*MU*pos(1)*pos(3)/radius_5 + 45/2*GravCoeff*pos(3)*pos(1)/radius_7 - 105/2*GravCoeff*pos(3)^3*pos(1)/radius_9;
    A(5,1) = A(4,2); % Symmetric
    A(5,2) = -MU/radius_3 + 3*MU*pos(2)^2/radius_5 + 15/2*GravCoeff*pos(3)^2/radius_7 - 105/2*GravCoeff*pos(3)^2*pos(2)^2/radius_9 - 3/2*GravCoeff/radius_5 + 15/2*GravCoeff*pos(2)^2/radius_7;
    A(5,3) = 3*MU*pos(2)*pos(3)/radius_5 + 45/2*GravCoeff*pos(3)*pos(2)/radius_7 - 105/2*GravCoeff*pos(3)^3*pos(2)/radius_9;
    A(6,1) = A(4,3); % Symmetric
    A(6,2) = A(5,3); % Symmetric
    A(6,3) = -MU/radius_3 + 3*MU*pos(3)^2/radius_5 + 45*GravCoeff*pos(3)^2/radius_7 - 105/2*GravCoeff*pos(3)^4/radius_9 - 9/2*GravCoeff/radius_5;

    % State partials derivative
    phi_xx_dot = A * phi_xx;

    % Stochastic partials derivative
    phi_xp_dot = A * phi_xp + fp;

    % Parameter partials derivative
    fc = zeros(Nx, Nc);
    fc(4:6, 1) = -pos / radius_3; % w.r.t. MU
    if Nc >= 2
        fc(4,1) = fc(4,1) + 15/2*J2*REarth^2*pos(3)^2*pos(1)/radius_7 - 3/2*J2*REarth^2*pos(1)/radius_5;
        fc(5,1) = fc(5,1) + 15/2*J2*REarth^2*pos(3)^2*pos(2)/radius_7 - 3/2*J2*REarth^2*pos(2)/radius_5;
        fc(6,1) = fc(6,1) + 15/2*J2*REarth^2*pos(3)^3/radius_7 - 9/2*J2*REarth^2*pos(3)/radius_5;
        fc(4,2) = 15/2*MU*REarth^2*pos(3)^2*pos(1)/radius_7 - 3/2*MU*REarth^2*pos(1)/radius_5;
        fc(5,2) = 15/2*MU*REarth^2*pos(3)^2*pos(2)/radius_7 - 3/2*MU*REarth^2*pos(2)/radius_5;
        fc(6,2) = 15/2*MU*REarth^2*pos(3)^3/radius_7 - 9/2*MU*REarth^2*pos(3)/radius_5;
    end
    phi_xc_dot = A * phi_xc + fc;

    % Combine derivatives
    xdot = [xdot; reshape(phi_xx_dot, Nx*Nx, 1); reshape(phi_xp_dot, Nx*(Np*N_batches), 1); reshape(phi_xc_dot, Nx*Nc, 1)];
end
