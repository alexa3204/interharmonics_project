%% Interharmonic Voltage Analysis and Phasor Extrqaction 
% Analyzes voltage waveform: V = 1.05*cos(2*pi*60*t + 10°) + 0.01*cos(2*pi*65*t + 10°)

close all; clear all; clc 

% Parameters
fs = 1000;              % Sampling frequency (Hz)
T = 1;                  % Analysis duration (seconds)
t = 0:1/fs:T-1/fs;      % Time vector
f1 = 60;                % Fundamental frequency (Hz)
f2 = 65;                % Interharmonic frequency (Hz)

% Original signal parameters
A1 = 1.05;              % Amplitude of 60 Hz component
A2 = 0.1;              % Amplitude of 65 Hz component
phase1 = 10 * pi/180;   % Phase of 60 Hz component (radians)
phase2 = 5 * pi/180;   % Phase of 65 Hz component (radians)

% Generate composite voltage waveform
V = A1 * cos(2*pi*f1*t + phase1) + A2 * cos(2*pi*f2*t + phase2);

% FFT-based frequency separation and phasor extraction 
N = length(V); 
f_axis = (0:N-1) * fs/N;  % Frequency axis 

% Compute FFT 
V_fft = fft(V); 
V_magnitude = abs(V_fft) * 2 / N; % Two-sided to single-sided conversion 
V_phase = angle(V_fft); 

% Find peaks at target frequencies 
[~, idx1] = min(abs(f_axis - f1));  % Index for 60 Hz
[~, idx2] = min(abs(f_axis - f2));  % Index for 65 Hz

% Extract phasors from FFT
phasor1_fft = V_fft(idx1) * 2 / N;  % 60 Hz phasor
phasor2_fft = V_fft(idx2) * 2 / N;  % 65 Hz phasor

% Reconstruct individual frequency components using inverse FFT
V_fft_60Hz = zeros(size(V_fft));
V_fft_60Hz(idx1) = V_fft(idx1);
V_fft_60Hz(N-idx1+2) = V_fft(N-idx1+2);  % Add conjugate for real signal
V1_reconstructed = real(ifft(V_fft_60Hz));

V_fft_65Hz = zeros(size(V_fft));
V_fft_65Hz(idx2) = V_fft(idx2);
V_fft_65Hz(N-idx2+2) = V_fft(N-idx2+2);  % Add conjugate for real signal
V2_reconstructed = real(ifft(V_fft_65Hz));

fprintf('=== VOLTAGE WAVEFORM ANALYSIS ===\n\n');

fprintf('Original Signal Parameters:\n');
fprintf('60 Hz: Amplitude = %.3f V, Phase = %.1f°\n', A1, phase1*180/pi);
fprintf('65 Hz: Amplitude = %.3f V, Phase = %.1f°\n\n', A2, phase2*180/pi);

fprintf('=== FFT PHASOR EXTRACTION RESULTS ===\n\n');

fprintf('60 Hz Phasor: %.3f∠%.1f° V (Magnitude: %.3f V)\n', ...
    abs(phasor1_fft), angle(phasor1_fft)*180/pi, abs(phasor1_fft));
fprintf('65 Hz Phasor: %.3f∠%.1f° V (Magnitude: %.3f V)\n\n', ...
    abs(phasor2_fft), angle(phasor2_fft)*180/pi, abs(phasor2_fft));

% Calculate extraction errors
mag_error_60Hz = abs(abs(phasor1_fft) - A1);
phase_error_60Hz = abs(angle(phasor1_fft) - phase1) * 180/pi;
mag_error_65Hz = abs(abs(phasor2_fft) - A2);
phase_error_65Hz = abs(angle(phasor2_fft) - phase2) * 180/pi;

fprintf('Extraction Errors:\n');
fprintf('60 Hz: Magnitude Error = %.6f V, Phase Error = %.3f°\n', ...
    mag_error_60Hz, phase_error_60Hz);
fprintf('65 Hz: Magnitude Error = %.6f V, Phase Error = %.3f°\n\n', ...
    mag_error_65Hz, phase_error_65Hz);

% Plotting
figure('Position', [100, 100, 1200, 800]);

% Time domain plots
subplot(3,2,1);
plot(t(1:500), V(1:500), 'b-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Voltage (V)');
title('Original Composite Voltage Waveform');
grid on;

subplot(3,2,2);
plot(t(1:500), V1_reconstructed(1:500), 'r-', 'LineWidth', 1.5);
hold on;
plot(t(1:500), V2_reconstructed(1:500), 'g-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Voltage (V)');
title('FFT-Reconstructed Frequency Components');
legend('60 Hz Component', '65 Hz Component');
grid on;

% Frequency domain plot
subplot(3,2,3);
plot(f_axis(1:N/2), V_magnitude(1:N/2), 'b-', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Magnitude (V)');
title('FFT Magnitude Spectrum');
xlim([50, 75]); grid on;

% Phasor diagram
subplot(3,2,4);
compass(real(phasor1_fft), imag(phasor1_fft), 'r');
hold on;
compass(real(phasor2_fft), imag(phasor2_fft), 'g');
title('Phasor Diagram (FFT Method)');
legend('60 Hz', '65 Hz', 'Location', 'best');

% Error analysis
subplot(3,2,5);
error_data = [mag_error_60Hz*1000, phase_error_60Hz; ...  % Scale magnitude error to mV
              mag_error_65Hz*1000000, phase_error_65Hz];   % Scale 65Hz mag error to µV
bar(error_data);
set(gca, 'XTickLabel', {'60 Hz', '65 Hz'});
ylabel('Error'); 
title('FFT Extraction Errors');
legend('Magnitude Error (mV/µV)', 'Phase Error (degrees)', 'Location', 'best');
grid on;

% Reconstructed vs Original
subplot(3,2,6);
V_reconstructed = V1_reconstructed + V2_reconstructed;
plot(t(1:200), V(1:200), 'b-', 'LineWidth', 2);
hold on;
plot(t(1:200), V_reconstructed(1:200), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Voltage (V)');
title('Original vs FFT Reconstructed');
legend('Original', 'FFT Reconstructed', 'Location', 'best');
grid on;

% Calculate reconstruction error
reconstruction_error = rms(V - V_reconstructed);
text(0.02, max(V(1:200))*0.8, sprintf('RMS Error: %.2e V', reconstruction_error), ...
     'FontSize', 10, 'BackgroundColor', 'white');

sgtitle('Interharmonic Voltage Analysis - FFT Method Only', 'FontSize', 14, 'FontWeight', 'bold');

%% Test with actual simulation waveform 
% Interharmonic Current Analysis and Phasor Extraction
% Analyzes simulated current waveform from Simulink

clear; clc; close all;

% Parameters
fs = 20000;             % Sampling frequency (Hz)
T = 1;                  % Analysis duration (seconds)
f1 = 60;                % Fundamental frequency (Hz)
f2 = 63;                % Interharmonic frequency (Hz) - updated to match your simulation

% Load simulated current waveform
% Run Simulink simulation
out = sim('three_bus_current_inj_bal.slx');
I = out.i3_data(:,1);  % Extract current data

% Create time vector based on actual data length
N = length(I);
t = (0:N-1) / fs;

fprintf('Successfully loaded simulation data:\n');
fprintf('Data points: %d\n', N);
fprintf('Duration: %.3f seconds\n', t(end));
fprintf('Sampling frequency: %d Hz\n', fs);

% FFT-based frequency separation and phasor extraction
f_axis = (0:N-1) * fs / N;  % Frequency axis

% Compute FFT
I_fft = fft(I);
I_magnitude = abs(I_fft) * 2 / N;  % Two-sided to single-sided conversion
I_phase = angle(I_fft);

% Find peaks at target frequencies
[~, idx1] = min(abs(f_axis - f1));  % Index for 60 Hz
[~, idx2] = min(abs(f_axis - f2));  % Index for 63 Hz

% Extract phasors from FFT
phasor1_fft = I_fft(idx1) * 2 / N;  % 60 Hz phasor
phasor2_fft = I_fft(idx2) * 2 / N;  % 63 Hz phasor

% Reconstruct individual frequency components using inverse FFT
I_fft_60Hz = zeros(size(I_fft));
I_fft_60Hz(idx1) = I_fft(idx1);
% Handle conjugate symmetry for real signals
if idx1 > 1 && idx1 <= N/2
    I_fft_60Hz(N-idx1+2) = I_fft(N-idx1+2);  % Add conjugate for real signal
end
I1_reconstructed = real(ifft(I_fft_60Hz));

I_fft_63Hz = zeros(size(I_fft));
I_fft_63Hz(idx2) = I_fft(idx2);
% Handle conjugate symmetry for real signals  
if idx2 > 1 && idx2 <= N/2
    I_fft_63Hz(N-idx2+2) = I_fft(N-idx2+2);  % Add conjugate for real signal
end
I2_reconstructed = real(ifft(I_fft_63Hz));

% Display results
fprintf('\n=== CURRENT WAVEFORM ANALYSIS ===\n\n');

fprintf('Signal Statistics:\n');
fprintf('RMS Value: %.3f A\n', rms(I));
fprintf('Peak Value: %.3f A\n', max(abs(I)));
fprintf('Data Points: %d\n', N);
fprintf('Frequency Resolution: %.3f Hz\n\n', fs/N);

fprintf('=== FFT PHASOR EXTRACTION RESULTS ===\n\n');

fprintf('60 Hz Phasor: %.3f∠%.1f° A (Magnitude: %.3f A)\n', ...
    abs(phasor1_fft), angle(phasor1_fft)*180/pi, abs(phasor1_fft));
fprintf('63 Hz Phasor: %.3f∠%.1f° A (Magnitude: %.3f A)\n\n', ...
    abs(phasor2_fft), angle(phasor2_fft)*180/pi, abs(phasor2_fft));

% Find actual frequency bin values
actual_f1 = f_axis(idx1);
actual_f2 = f_axis(idx2);
fprintf('Frequency Bins Used:\n');
fprintf('Target 60 Hz -> Actual %.2f Hz (bin %d)\n', actual_f1, idx1);
fprintf('Target 63 Hz -> Actual %.2f Hz (bin %d)\n\n', actual_f2, idx2);

% Plotting
figure('Position', [100, 100, 1200, 800]);

% Determine plot range based on data length
plot_samples = min(10000, N);  % Plot first 1000 samples or all data if less

% Time domain plots
subplot(3,2,1);
plot(t(1:plot_samples), I(1:plot_samples), 'b-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Current (A)');
title('Simulated Current Waveform');
grid on;

subplot(3,2,2);
plot(t(1:plot_samples), I1_reconstructed(1:plot_samples), 'r-', 'LineWidth', 1.5);
hold on;
plot(t(1:plot_samples), I2_reconstructed(1:plot_samples), 'g-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Current (A)');
title('FFT-Reconstructed Frequency Components');
legend('60 Hz Component', '63 Hz Component');
grid on;

% Frequency domain plot
subplot(3,2,3);
plot(f_axis(1:N/2), I_magnitude(1:N/2), 'b-', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Magnitude (A)');
title('FFT Magnitude Spectrum');
xlim([55, 70]); grid on;
% Mark the detected frequencies
hold on;
plot(actual_f1, I_magnitude(idx1), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(actual_f2, I_magnitude(idx2), 'go', 'MarkerSize', 8, 'LineWidth', 2);
legend('Spectrum', '60 Hz', '63 Hz', 'Location', 'best');

% Phasor diagram
subplot(3,2,4);
compass(real(phasor1_fft), imag(phasor1_fft), 'r');
hold on;
compass(real(phasor2_fft), imag(phasor2_fft), 'g');
title('Phasor Diagram (FFT Method)');
legend('60 Hz', '63 Hz', 'Location', 'best');

% Spectral analysis around target frequencies
subplot(3,2,5);
freq_range = 55:0.1:70;
spectrum_interp = interp1(f_axis(1:N/2), I_magnitude(1:N/2), freq_range, 'linear', 0);
semilogy(freq_range, spectrum_interp, 'b-', 'LineWidth', 1.5);
hold on;
semilogy(actual_f1, I_magnitude(idx1), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
semilogy(actual_f2, I_magnitude(idx2), 'go', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Magnitude (A) - Log Scale');
title('Detailed Spectrum (55-70 Hz)');
legend('Spectrum', '60 Hz Peak', '63 Hz Peak', 'Location', 'best');
grid on;

% Reconstructed vs Original
subplot(3,2,6);
I_reconstructed = I1_reconstructed + I2_reconstructed;
plot_range = min(10000, N);  % Plot first 400 samples for clarity
plot(t(1:plot_range), I(1:plot_range), 'b-', 'LineWidth', 2);
hold on;
plot(t(1:plot_range), I_reconstructed(1:plot_range), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Current (A)');
title('Original vs FFT Reconstructed');
legend('Original', 'FFT Reconstructed', 'Location', 'best');
grid on;

% Calculate reconstruction error
reconstruction_error = rms(I - I_reconstructed);
reconstruction_percent = 100 * reconstruction_error / rms(I);
text(0.02*max(t(1:plot_range)), max(I(1:plot_range))*0.8, ...
     sprintf('RMS Error: %.2e A (%.2f%%)', reconstruction_error, reconstruction_percent), ...
     'FontSize', 9, 'BackgroundColor', 'white');

sgtitle('Simulated Current Interharmonic Analysis - FFT Method', 'FontSize', 14, 'FontWeight', 'bold');


%% Test phasor extraction on system without harmonics 
% test to see what the angle output would be from the three-bus system
% running at 60hz 

close all; clear all; clc 

out = sim('three_bus_fundamental.slx');

fs = 20000; % sampling frequency 
duration = 1; % test duration (s) 
t = 0:1/fs:duration; 

% Fundamental frequency
f1 = 60;           % Hz
omega1 = 2*pi*f1;  % rad/s

% get current vectors from phase a at all buses 
I = zeros(3, length(t)); 
I(1,:) = out.i1_data(:,1); 
I(2,:) = out.i2_data(:,1); 
I(3,:) = out.i3_data(:,1); 

% FFT-based frequency phasor extraction 
N = length(I(1,:));
f_axis = (0:N-1) * fs / N; % Frequency axis 

for i = 1:3 
    % Compute FFT 
    I_fft(i,:) = fft(I(i,:));

    [~, idx1] = min(abs(f_axis - f1)); % Index for 60Hz

    % Extract phasors from FFT 
    current_phasor_fft(i,1) = I_fft(i,idx1) * 2/N; 

end 

fprintf('Current Phasors:\n');
[abs(current_phasor_fft) rad2deg(angle(current_phasor_fft))]

% Voltage phasors
V = zeros(3, length(t)); 
V(1,:) = out.v1_data(:,1); 
V(2,:) = out.v2_data(:,1); 
V(3,:) = out.v3_data(:,1); 

for i = 1:3 
    V_fft(i,:) = fft(V(i,:));
    [~, idx1] = min(abs(f_axis - f1));
    voltage_phasor_fft(i,1) = V_fft(i,idx1) * 2/N; 
end 

fprintf('Voltage Phasors:\n');
[abs(voltage_phasor_fft) rad2deg(angle(voltage_phasor_fft))]

%% Test phasor extraction on system with harmonics 
% 60hz with current injection at 63hz 

close all; clear all; clc 

out = sim('three_bus_current_inj_bal.slx');

fs = 20000; % sampling frequency 
duration = 1; % test duration (s) 
t = 0:1/fs:duration; 

% Fundamental frequency
f1 = 63;           % Hz
omega1 = 2*pi*f1;  % rad/s

% get current vectors from phase a at all buses 
I = zeros(3, length(t)); 
I(1,:) = out.i1_data(:,1); 
I(2,:) = out.i2_data(:,1); 
I(3,:) = out.i3_data(:,1); 

% FFT-based frequency phasor extraction 
N = length(I(1,:));
f_axis = (0:N-1) * fs / N; % Frequency axis 

for i = 1:3 
    % Compute FFT 
    I_fft(i,:) = fft(I(i,:));

    [~, idx1] = min(abs(f_axis - f1)); % Index for 60Hz

    % Extract phasors from FFT 
    current_phasor_fft(i,1) = I_fft(i,idx1) * 2/N; 

end 

fprintf('Current Phasors:\n');
[abs(current_phasor_fft) rad2deg(angle(current_phasor_fft))]

% Voltage phasors
V = zeros(3, length(t)); 
V(1,:) = out.v1_data(:,1); 
V(2,:) = out.v2_data(:,1); 
V(3,:) = out.v3_data(:,1); 

for i = 1:3 
    V_fft(i,:) = fft(V(i,:));
    [~, idx1] = min(abs(f_axis - f1));
    voltage_phasor_fft(i,1) = V_fft(i,idx1) * 2/N; 
end 

fprintf('Voltage Phasors:\n');
[abs(voltage_phasor_fft) rad2deg(angle(voltage_phasor_fft))]

%% Test PMU Ybus validity  

close all; clear all; clc 

out = sim('PMU_three_bus_fundamental.slx');

fs = 20000; % sampling frequency 
duration = 1; % test duration (s) 
t = 0:1/fs:duration; 

% Fundamental frequency
f1 = 60;           % Hz
omega1 = 2*pi*f1;  % rad/s

% get current vectors from phase a at all buses 
I(1,:) = out.i1_data(end,:); 
I(2,:) = out.i2_data(end,:); 
I(3,:) = out.i3_data(end,:); 

% get current vectors from phase a at all buses 
V(1,:) = out.v1_data(end,:); 
V(2,:) = out.v2_data(end,:); 
V(3,:) = out.v3_data(end,:); 


% convert from phasor to rectangular 
I_c = I(:,1).*exp(j.*(deg2rad(I(:,2))));
V_c = V(:,1).*exp(j.*(deg2rad(V(:,2))));


%% Line admittances
% Y-bus

% Define base values
f = 60;
w = 2 * pi * f;

v_base = 10e3;
s_base = 100e6;
z_base = v_base^2 / s_base;

length = 1; % line length (km)
% Initialize 3x3 Ybus
Ybus = complex(zeros(3, 3));

% Define line data (buses and R, L, C)
% Format: [from to R L C] (value per km)
pi_lines = [
    1 2 0.05 0.0013528 2.8011e-05;     % line 1-2
    1 3 0.1 0.0029285 2.5146e-05;    % line 1-3
    2 3 0.2 0.0051248 4.8701e-05;   % line 2-3
    ];

for k = 1:size(pi_lines, 1)
    i = pi_lines(k, 1);
    g = pi_lines(k, 2);
    R = pi_lines(k, 3) * length;
    L = pi_lines(k, 4) * length;
    C = pi_lines(k, 5) * length;
 
    Z = (R + 1j * w * L)/z_base;
    Yseries = 1 / Z;

    Yshunt = 1j * w * C * z_base;

    % Add to Ybus
    Ybus(i,i) = Ybus(i,i) + Yseries + Yshunt/2;
    Ybus(g,g) = Ybus(g,g) + Yseries + Yshunt/2;
    Ybus(i,g) = Ybus(i,g) - Yseries;
    Ybus(g,i) = Ybus(g,i) - Yseries;
end

Y_r = Ybus;