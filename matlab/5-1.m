% Defining the system parameters
% Assuming a linear system with 5 degrees of freedom
M = eye(5);  % Mass matrix
C = diag([3.248, 20.355, 56.994, 111.685, 184.623]);  % Damping matrix
K = diag([263.804, 10358.158, 81206.525, 311837.785, 852137.646]);  % Stiffness matrix

% Frequency range for analysis
freq_range = linspace(0.1, 100, 500);  % 0.1 to 100 Hz, 500 points

% Preallocating the frequency response matrix
H = zeros(5, 5, length(freq_range));

% Calculating the frequency response
for i = 1:length(freq_range)
    % Frequency in rad/s
    omega = 2 * pi * freq_range(i);

    % Dynamic stiffness matrix at the frequency
    K_dyn = -omega^2 * M + 1i * omega * C + K;

    % Calculating the inverse of the dynamic stiffness matrix
    H(:, :, i) = inv(K_dyn);
end

% Plotting the magnitude and phase of the frequency response
figure;

% Magnitude plot
subplot(2, 1, 1);
for i = 1:5
    plot(freq_range, 20 * log10(abs(squeeze(H(i, i, :)))), 'DisplayName', sprintf('Mode %d', i));
    hold on;
end
title('Magnitude of Frequency Response');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
legend;
grid on;

% Phase plot
subplot(2, 1, 2);
for i = 1:5
    plot(freq_range, angle(squeeze(H(i, i, :))), 'DisplayName', sprintf('Mode %d', i));
    hold on;
end
title('Phase of Frequency Response');
xlabel('Frequency [Hz]');
ylabel('Phase [Radians]');
legend;
grid on;

% Adjust layout
sgtitle('Frequency Response of the System');

