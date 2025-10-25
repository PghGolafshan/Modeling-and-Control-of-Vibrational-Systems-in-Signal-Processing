% Constants
L = 0.25; % Length in meters (m)
betas = [7.500, 18.776, 31.419, 43.982, 56.549]; % Beta values for each mode
zeta = 0.1; % Damping ratio for all modes
natural_frequencies = [2.585, 16.198, 45.354, 88.876, 146.918]; % w1 to w5
positions = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]; % 10 positions
% Calculate Wi(x) for each position and beta
wi_values = zeros(length(positions), length(betas));
for i = 1:length(positions)
   x = positions(i);
   for j = 1:length(betas)
       beta = betas(j);
       wi_values(i, j) = (1 / (sin(beta * L) - sinh(beta * L))) * ...
                        ((sin(beta * L) - sinh(beta * L)) * (sin(beta * x) - sinh(beta * x)) + ...
                        (cos(beta * L) + cosh(beta * L)) * (cos(beta * x) - cosh(beta * x)));
   end
end
% Frequency range to evaluate
frequency_range = linspace(0.1, 150, 1000); % 0.1 Hz to 150 Hz
% Calculate and plot the frequency response for each position and mode
for i = 1:length(positions)
   figure;
   for j = 1:length(natural_frequencies)
       wn = natural_frequencies(j);
       wi = wi_values(i, j);
       response_amplitudes = zeros(size(frequency_range));
      
       for k = 1:length(frequency_range)
           frequency = frequency_range(k);
           % Calculate the amplitude of the steady-state response
           r = frequency / wn;
           denominator = sqrt((1 - r^2)^2 + (2*zeta*r)^2);
           amplitude = abs(wi / denominator);
           response_amplitudes(k) = amplitude;
       end
      
       plot(frequency_range, response_amplitudes, 'DisplayName', sprintf('Mode %d (w=%f Hz)', j, wn));
       hold on;
   end
   title(sprintf('Frequency Response at Position %f meters', positions(i)));
   xlabel('Frequency (Hz)');
   ylabel('Amplitude');
   legend;
   hold off;
end

