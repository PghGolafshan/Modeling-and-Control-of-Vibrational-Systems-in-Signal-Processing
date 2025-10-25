function sinesumconv
    % Constants
    beta_values = [7.5, 18.776, 31.419, 43.982, 56.549];
    L = 0.25;
    N = 5;
    a = [3.248, 20.355, 56.994, 111.685, 184.623];
    b = [263.804, 10358.158, 81206.525, 311837.785, 852137.646];
    c = [0.1136, 0.4131, 0.9280, 1.3417, 1.512];
    w = [16.2338, 101.72344, 284.82312, 558.14128, 922.64504];

    % Time span for the impulse response
    t_impulse_span = linspace(0, 2, 100);
    q_values_impulse = zeros(N, length(t_impulse_span));
    x_range = linspace(0, L, 100);
    y_values_impulse_2d = zeros(length(x_range), length(t_impulse_span));

    % Solve for each q_i(t) with impulse input
    for i = 1:N
        sol = ode45(@(t, q) odefun_impulse(t, q, a(i), b(i), c(i)), [0, 2], [0, 0]);
        q_values_impulse(i, :) = deval(sol, t_impulse_span, 1);
        for j = 1:length(x_range)
            w_i_value = w_i(x_range(j), beta_values(i), L);
            y_values_impulse_2d(j, :) = y_values_impulse_2d(j, :) + w_i_value * q_values_impulse(i, :);
        end
    end

    % Time span for the sinusoidal response
    t_sin_span = linspace(0, 2, 500);
    q_values_sin = zeros(N, length(t_sin_span), length(w));
    y_values_sin_2d = zeros(length(x_range), length(t_sin_span), length(w));

    % Compute the response to sinusoidal inputs
    for k = 1:length(w)
        sin_input = sin(w(k) * t_sin_span);
        for i = 1:N
            conv_result = conv(q_values_impulse(i, :), sin_input, 'same');
            dt = mean(diff(t_sin_span));
            q_values_sin(i, :, k) = conv_result(1:length(t_sin_span)) * dt;
            for j = 1:length(x_range)
                y_values_sin_2d(j, :, k) = y_values_sin_2d(j, :, k) + w_i(x_range(j), beta_values(i), L) * q_values_sin(i, :, k);
            end
        end
    end

    % Plotting the sinusoidal response for each frequency
    figure;
    for k = 1:length(w)
        subplot(3, 2, k);
        plot(t_sin_span, y_values_sin_2d(round(length(x_range) / 2), :, k));
        title(sprintf('Sinusoidal Response at x = %.3f, Frequency w%d = %.2f rad/s', x_range(round(length(x_range) / 2)), k, w(k)));
        xlabel('Time (t)');
        ylabel(sprintf('y(x=%.3f, t)', x_range(round(length(x_range) / 2))));
        grid on;
    end
end

% Define the ODE function for impulse response
function dq = odefun_impulse(t, q, a, b, c)
    A = impulse_approx(t);
    dq = [q(2); A * c - a * q(2) - b * q(1)];
end

% Define the w_i function
function w = w_i(x, beta, L)
    sin_betaL = sin(beta * L);
    sinh_betaL = sinh(beta * L);
    cos_betaL = cos(beta * L);
    cosh_betaL = cosh(beta * L);
    w = (1 / (sin_betaL - sinh_betaL)) * ...
        ((sin_betaL - sinh_betaL) * (sin(beta * x) - sinh(beta * x)) + ...
        (cos_betaL + cosh_betaL) * (cos(beta * x) - cosh(beta * x)));
end

% Approximate an impulse function
function A = impulse_approx(t, duration)
    if nargin < 2
        duration = 0.01;
    end
    A = 1 / duration * (t >= 0 & t < duration);
end
