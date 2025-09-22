close all; clear all; clc 

out = sim('three_bus_current_inj_a.slx');

fs = 20000; % sampling frequency 
duration = 3; % test duration (s) 
t = 0:1/fs:duration; 

% Fundamental frequency
f1 = 60;           % Hz
omega1 = 2*pi*f1;  % rad/s

% Oscillation parameters
fos = 3;           % Oscillation frequency in Hz (typical for power system oscillations: 0.1-3 Hz)
omega_os = 2*pi*fos;
f_upper = f1 + fos; 

% focussing on bus 3 phase a 
i3a = out.i3_data(:,1);

% calculate envelope of bus 3 phase a
[envelope_upper, envelope_lower] = envelope(i3a);

% calculate phasor magnitude (what PMU would measure)
% Using sliding window RMS calculation
window_samples = round(fs/f1);  % One fundamental cycle
phasor_mag = zeros(size(t));
for i = window_samples:length(t)
    window_data = i3a(i-window_samples+1:i);
    phasor_mag(i) = sqrt(mean(window_data.^2));
end

figure('Position', [100 100 1200 800]);
% Subplot 1: Complete waveform with envelope
subplot(2,1,1);
plot(t, i3a, 'b-', 'LineWidth', 1);
hold on; 
plot(t, envelope_upper, 'r--', 'LineWidth', 2);
plot(t, envelope_lower, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Current (pu)');
title('Complete Current Waveform');
legend('Total Waveform', 'Envelope (Beating Pattern)', 'Location', 'northeast');

% Subplot 2: Individual components
subplot(2,1,2);
plot(t, phasor_mag, 'k-', 'LineWidth', 2);
hold on;
grid on;
xlabel('Time (s)');
ylabel('RMS Voltage (V)');
title('Phasor Magnitude (PMU Measurement)');
legend('Measured Phasor RMS', 'Location', 'northeast');
xlim([0 duration]);

%% DFT
figure('Position', [100 100 800 600]);
N = length(i3a);
f_axis = (0:N-1) * fs / N;
V_fft = abs(fft(i3a))/ N * 2;  % Two-sided to one-sided conversion

% Plot spectrum
semilogy(f_axis(1:N/2), V_fft(1:N/2), 'b-', 'LineWidth', 1.5);
hold on;
% Mark the key frequencies
plot(f1, V_fft(round(f1*N/fs)), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(f_upper, V_fft(round(f_upper*N/fs)), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
%plot(f_lower, V_fft(round(f_lower*N/fs)), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');

grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('Frequency Spectrum - Shows Fundamental and Interharmonics');
legend('Spectrum', '60 Hz', '63 Hz', '57 Hz', 'Location', 'northeast');
xlim([20 70]);
ylim([10e-5 10]);