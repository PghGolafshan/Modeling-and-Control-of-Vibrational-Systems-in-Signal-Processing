% Constants
natural_frequencies = [2.585, 16.198, 45.354, 88.876, 146.918]; % w1 to w5 in Hz
wi_values = [0.1136, 0.4131, 0.9280, 1.3417, 1.512]; % Given Wi values
zeta_values = linspace(0, 0.1, 10); % 10 zeta values from 0 to 0.1
t = linspace(0, 100, 2000); % Time array, change as needed
A = 1; % Amplitude of the forcing function A(t) = sin(t)

% Frequency range for the frequency response
frequency_range = linspace(0.1, 200, 500); % Adjust as needed

% Iterate over each zeta value
for zeta = zeta_values
    
    % Initialize the plot for time response
    figure;
    hold on;
    
    % Solve the differential equation for time response for each mode
    for idx = 1:length(natural_frequencies)
        wn = natural_frequencies(idx);
        wi = wi_values(idx);

        % Differential equation solver
        [T, Y] = ode45(@(t,y) [y(2); -2*zeta*wn*y(2) - wn^2*y(1) + wi*A*sin(t)], t, [0 0]);
        
        % Plotting time response
        plot(T, Y(:,1), 'DisplayName', sprintf('Mode %d (w=%.3f Hz)', idx, wn));
    end
    
    title(sprintf('Time Response for ζ=%.3f', zeta));
    xlabel('Time (seconds)');
    ylabel('Displacement');
    legend;
    grid on;
    hold off;
    
    % Initialize the plot for frequency response
    figure;
    hold on;
    
    % Calculate and plot frequency response for each mode
    for idx = 1:length(natural_frequencies)
        wn = natural_frequencies(idx);
        wi = wi_values(idx);
        response_amplitudes = zeros(size(frequency_range));
        
        for k = 1:length(frequency_range)
            frequency = frequency_range(k);
            r = frequency / wn;
            denominator = sqrt((1 - r^2)^2 + (2*zeta*r)^2);
            amplitude = abs(wi / denominator);
            response_amplitudes(k) = amplitude;
        end
        
        % Plotting frequency response
        plot(frequency_range, response_amplitudes, 'DisplayName', sprintf('Mode %d (w=%.3f Hz)', idx, wn));
    end
    
    title(sprintf('Frequency Response for ζ=%.3f', zeta));
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    legend;
    grid on;
    hold off;
end
