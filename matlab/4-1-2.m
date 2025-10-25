function r_response
    % Constants
    beta_values = [7.5, 18.776, 31.419, 43.982, 56.549];
    L = 0.25;
    N = 5;
    w = [16.2338, 101.72344, 284.82312, 558.14128, 922.64504];
    a = [3.248, 20.355, 56.994, 111.685, 184.623];
    b = [263.804, 10358.158, 81206.525, 311837.785, 852137.646];
    c = [0.1136, 0.4131, 0.9280, 1.3417, 1.512];
    t_span = [0, 10];
    t_values = linspace(t_span(1), t_span(2), 1000);
    q_values = zeros(N, length(t_values));

    for i = 1:N
        initial_conditions = [0, 0];
        [t, sol] = ode45(@(t, q) odefun_double_impulse(t, q, a(i), b(i), c(i), w), t_span, initial_conditions);
        q_values(i, :) = interp1(t, sol(:, 1), t_values);
    end

    x_range = linspace(0, L, 100);
    y_values_2d = zeros(length(x_range), length(t_values));

    for j = 1:length(x_range)
        for i = 1:N
            w_i_value = w_i(x_range(j), beta_values(i), L);
            y_values_2d(j, :) = y_values_2d(j, :) + w_i_value .* q_values(i, :);
        end
    end

    plot_results(t_values, q_values, x_range, y_values_2d, N);
end

function dydt = odefun_double_impulse(t, q, a, b, c, w)
    % Approximation of double impulse function (r = delta(t) * delta(t))
    if abs(t) < 1e-6 % Even narrower than the single impulse
        A = 1e12; % Higher amplitude than the single impulse
    else
        A = 0;
    end
    dydt = [q(2); A * c - a * q(2) - b * q(1)];
end

function val = w_i(x, beta, L)
    sin_betaL = sin(beta * L);
    sinh_betaL = sinh(beta * L);
    cos_betaL = cos(beta * L);
    cosh_betaL = cosh(beta * L);
    val = (1 ./ (sin_betaL - sinh_betaL)) .* ...
        ((sin_betaL - sinh_betaL) .* (sin(beta * x) - sinh(beta * x)) + ...
         (cos_betaL + cosh_betaL) .* (cos(beta * x) - cosh(beta * x)));
end

function plot_results(t_values, q_values, x_range, y_values_2d, N)
    figure;
    for i = 1:N
        subplot(3, 3, i);
        plot(t_values, q_values(i, :));
        title(sprintf('q_%d(t)', i));
        xlabel('Time (t)');
        ylabel(sprintf('q_%d(t)', i));
        grid on;
    end

    subplot(3, 3, N + 1);
    selected_time_points = [10, 20, 30, 40, 50];
    for time_point = selected_time_points
        [~, idx] = min(abs(t_values - time_point));
        plot(x_range, y_values_2d(:, idx), 'DisplayName', sprintf('t = %.1f s', time_point));
        hold on;
    end
    title('y(x, t) for Different Time Points');
    xlabel('x');
    ylabel('y(x, t)');
    legend;
    grid on;
    hold off;
end
