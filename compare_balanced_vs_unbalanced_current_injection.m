close all; clear all; clc

% Run both simulations
fprintf('Running first simulation...\n');
out = sim('three_bus_current_inj_a.slx');
fprintf('Running second simulation...\n');
out_2 = sim('three_bus_current_inj_bal.slx'); 

% Simulation parameters
fs = 20000; % sampling frequency
duration = 3; % test duration (s)
t = 0:1/fs:duration;

% Fundamental frequency
f1 = 60; % Hz
omega1 = 2*pi*f1; % rad/s

% Oscillation parameters
fos = 3; % Oscillation frequency in Hz (typical for power system oscillations: 0.1-3 Hz)
omega_os = 2*pi*fos;
f_upper = f1 + fos;

% Extract phase A currents from both simulations - bus 3
i3a_sim1 = out.i3_data(:,1);
i3a_sim2 = out_2.i3_data(:,1);

% Ensure both signals have the same length
min_length = min(length(i3a_sim1), length(i3a_sim2));
i3a_sim1 = i3a_sim1(1:min_length);
i3a_sim2 = i3a_sim2(1:min_length);
t_plot = t(1:min_length);

% Calculate difference signal
i3a_diff = i3a_sim1 - i3a_sim2;

% Calculate envelopes for both signals
[envelope_upper_1, envelope_lower_1] = envelope(i3a_sim1);
[envelope_upper_2, envelope_lower_2] = envelope(i3a_sim2);

% Calculate phasor magnitudes (what PMU would measure)
window_samples = round(fs/f1); % One fundamental cycle
phasor_mag_1 = zeros(size(t_plot));
phasor_mag_2 = zeros(size(t_plot));

for i = window_samples:length(t_plot)
    window_data_1 = i3a_sim1(i-window_samples+1:i);
    window_data_2 = i3a_sim2(i-window_samples+1:i);
    phasor_mag_1(i) = sqrt(mean(window_data_1.^2));
    phasor_mag_2(i) = sqrt(mean(window_data_2.^2));
end

%% Statistical comparison
fprintf('\n=== SIMULATION COMPARISON RESULTS ===\n');
fprintf('Maximum absolute difference: %.6f\n', max(abs(i3a_diff)));
fprintf('RMS of difference signal: %.6f\n', sqrt(mean(i3a_diff.^2)));
fprintf('Mean of difference signal: %.6f\n', mean(i3a_diff));
fprintf('Standard deviation of difference: %.6f\n', std(i3a_diff));

% Correlation coefficient
correlation = corrcoef(i3a_sim1, i3a_sim2);
fprintf('Correlation coefficient: %.8f\n', correlation(1,2));

% Relative error metrics
relative_error_rms = sqrt(mean(i3a_diff.^2)) / sqrt(mean(i3a_sim1.^2)) * 100;
fprintf('Relative RMS error: %.4f%%\n', relative_error_rms);

%% Plotting
% Figure 1: Waveform comparison
figure('Position', [100 100 1400 900]);

% Subplot 1: Both complete waveforms
subplot(3,1,1);
plot(t_plot, i3a_sim1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Simulation 1');
hold on;
plot(t_plot, i3a_sim2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Simulation 2');
grid on;
xlabel('Time (s)');
ylabel('Current (pu)');
title('Waveform Comparison: Both Simulations');
legend('Location', 'northeast');
xlim([0 duration]);

% Subplot 2: Difference signal
subplot(3,1,2);
plot(t_plot, i3a_diff, 'g-', 'LineWidth', 1);
grid on;
xlabel('Time (s)');
ylabel('Current Difference (pu)');
title('Difference Signal (Simulation 1 - Simulation 2)');
xlim([0 duration]);

% Subplot 3: Phasor magnitude comparison
subplot(3,1,3);
plot(t_plot, phasor_mag_1, 'b-', 'LineWidth', 2, 'DisplayName', 'Simulation 1 RMS');
hold on;
plot(t_plot, phasor_mag_2, 'r--', 'LineWidth', 2, 'DisplayName', 'Simulation 2 RMS');
grid on;
xlabel('Time (s)');
ylabel('RMS Current (pu)');
title('Phasor Magnitude Comparison (PMU Measurements)');
legend('Location', 'northeast');
xlim([0 duration]);

% Figure 2: Detailed waveform with envelopes
figure('Position', [200 200 1400 800]);

% Subplot 1: Simulation 1 with envelope
subplot(2,1,1);
plot(t_plot, i3a_sim1, 'b-', 'LineWidth', 1);
hold on;
plot(t_plot, envelope_upper_1, 'r--', 'LineWidth', 2);
plot(t_plot, envelope_lower_1, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Current (pu)');
title('Simulation 1: Current Waveform with Envelope');
legend('Waveform', 'Envelope', 'Location', 'northeast');
xlim([0 duration]);

% Subplot 2: Simulation 2 with envelope
subplot(2,1,2);
plot(t_plot, i3a_sim2, 'b-', 'LineWidth', 1);
hold on;
plot(t_plot, envelope_upper_2, 'r--', 'LineWidth', 2);
plot(t_plot, envelope_lower_2, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Current (pu)');
title('Simulation 2: Current Waveform with Envelope');
legend('Waveform', 'Envelope', 'Location', 'northeast');
xlim([0 duration]);

%% Frequency domain comparison
figure('Position', [300 300 1200 800]);

N = length(i3a_sim1);
f_axis = (0:N-1) * fs / N;

% Calculate FFTs
V_fft_1 = abs(fft(i3a_sim1)) / N * 2; % Two-sided to one-sided conversion
V_fft_2 = abs(fft(i3a_sim2)) / N * 2;
V_fft_diff = abs(fft(i3a_diff)) / N * 2;

% Subplot 1: Frequency spectra comparison
subplot(2,1,1);
semilogy(f_axis(1:N/2), V_fft_1(1:N/2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Simulation 1');
hold on;
semilogy(f_axis(1:N/2), V_fft_2(1:N/2), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Simulation 2');

% Mark key frequencies
plot(f1, V_fft_1(round(f1*N/fs)), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot(f1, V_fft_2(round(f1*N/fs)), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(f_upper, V_fft_1(round(f_upper*N/fs)), 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot(f_upper, V_fft_2(round(f_upper*N/fs)), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum Comparison');
legend('Location', 'northeast');
xlim([20 70]);
ylim([10e-6 10]);

% Subplot 2: Difference spectrum
subplot(2,1,2);
semilogy(f_axis(1:N/2), V_fft_diff(1:N/2), 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Difference Signal');
xlim([20 70]);

%% Zoomed comparison (first 0.5 seconds for detail)
figure('Position', [400 400 1200 600]);
zoom_end = round(0.5 * fs);
t_zoom = t_plot(1:zoom_end);

subplot(2,1,1);
plot(t_zoom, i3a_sim1(1:zoom_end), 'b-', 'LineWidth', 2, 'DisplayName', 'Simulation 1');
hold on;
plot(t_zoom, i3a_sim2(1:zoom_end), 'r--', 'LineWidth', 2, 'DisplayName', 'Simulation 2');
grid on;
xlabel('Time (s)');
ylabel('Current (pu)');
title('Detailed Waveform Comparison (First 0.5 seconds)');
legend('Location', 'northeast');

subplot(2,1,2);
plot(t_zoom, i3a_diff(1:zoom_end), 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Current Difference (pu)');
title('Difference Signal (Zoomed)');

fprintf('\n=== Analysis Complete ===\n');
if max(abs(i3a_diff)) < 1e-6
    fprintf('RESULT: Simulations appear to be identical (difference < 1μ units)\n');
elseif relative_error_rms < 0.1
    fprintf('RESULT: Simulations are very similar (relative error < 0.1%%)\n');
elseif relative_error_rms < 1.0
    fprintf('RESULT: Simulations have minor differences (relative error < 1%%)\n');
else
    fprintf('RESULT: Simulations have significant differences (relative error >= 1%%)\n');
end