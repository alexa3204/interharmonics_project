%% Interharmonic Power-Based Oscillation Source Location
% Implementation of Section V methodology for a 3-machine system

clear all; close all; clc;

%% System Parameters (Based on Fig. 24)
f1 = 60;                    % Fundamental frequency (Hz)
omega1 = 2*pi*f1;           % Angular frequency (rad/s)

% Grid impedance (Bus 1)
R_grid = 0.01;
L_grid = 0.001;

% Series capacitor compensation
C_series = 100e-6;          % Series capacitor

% Line impedances
R_line12 = 0.05;
L_line12 = 0.01;
R_line23 = 0.1;
L_line23 = 0.02;
R_line24 = 0.1;
L_line24 = 0.02;

% Induction Generator Parameters (simplified model)
% All three IGs are identical
R_s = 0.01;                 % Stator resistance
L_s = 0.05;                 % Stator inductance
R_r = 0.015;                % Rotor resistance
L_r = 0.05;                 % Rotor inductance
L_m = 1.5;                  % Magnetizing inductance

%% Step 1: 2D Frequency Scan to Find Poles

% Define scan range
sigma_range = linspace(-5, 5, 100);     % Real part of s
freq_range = linspace(30, 80, 200);     % Frequency range (Hz)

% Initialize result matrix
Z_min = zeros(length(sigma_range), length(freq_range));

for i = 1:length(sigma_range)
    for j = 1:length(freq_range)
        sigma = sigma_range(i);
        f = freq_range(j);
        s = sigma + 1j*2*pi*f;
        
        % Construct Y matrix at this s value
        Y = construct_Y_matrix(s, R_grid, L_grid, C_series, ...
                               R_line12, L_line12, R_line23, L_line23, ...
                               R_line24, L_line24, R_s, L_s, R_r, L_r, L_m);
        
        % Calculate minimum eigenvalue magnitude
        eigenvalues = eig(Y);
        Z_min(i,j) = min(abs(eigenvalues));
    end
end

% Find the pole (minimum point)
[min_val, idx] = min(Z_min(:));
[idx_sigma, idx_freq] = ind2sub(size(Z_min), idx);
sigma_pole = sigma_range(idx_sigma);
f_pole = freq_range(idx_freq);
s_pole = sigma_pole + 1j*2*pi*f_pole;

fprintf('\n=== Pole Location Results ===\n');
fprintf('Modal frequency (f_pole): %.2f Hz\n', f_pole);
fprintf('Damping (sigma_pole): %.4f\n', sigma_pole);
fprintf('Interharmonic frequency: %.2f Hz\n', f_pole);

% Calculate oscillation frequency using Eq. 3.3
h = round(f_pole/f1);  % Nearest harmonic number
f_os = abs(f_pole - h*f1);
fprintf('Phasor oscillation frequency (f_os): %.2f Hz\n', f_os);

if sigma_pole > 0
    fprintf('System is UNSTABLE\n');
else
    fprintf('System is STABLE\n');
end

%% Plot 2D Frequency Scan Results

figure('Position', [100 100 1200 400]);

% Surface plot
subplot(1,2,1);
surf(freq_range, sigma_range, Z_min);
xlabel('Frequency (Hz)');
ylabel('\sigma (damping)');
zlabel('|min(\lambda)|');
title('2D Frequency Scan - Surface Plot');
colorbar;
view(45, 30);
hold on;
plot3(f_pole, sigma_pole, min_val, 'ro', 'MarkerSize', 15, 'LineWidth', 3);
legend('System Response', 'Pole Location');

% Contour plot
subplot(1,2,2);
contourf(freq_range, sigma_range, log10(Z_min+1e-10), 50);
xlabel('Frequency (Hz)');
ylabel('\sigma (damping)');
title('2D Frequency Scan - Contour Plot (log scale)');
colorbar;
hold on;
plot(f_pole, sigma_pole, 'r*', 'MarkerSize', 20, 'LineWidth', 3);
plot([f1 f1], [min(sigma_range) max(sigma_range)], 'w--', 'LineWidth', 2);
legend('log_{10}|min(\lambda)|', 'Pole', 'f_1 = 60Hz');
grid on;

%% Step 2: Calculate Generator Participation Factors (Eq. 5.5)

% Construct Y matrix at the pole
Y_pole = construct_Y_matrix(s_pole, R_grid, L_grid, C_series, ...
                            R_line12, L_line12, R_line23, L_line23, ...
                            R_line24, L_line24, R_s, L_s, R_r, L_r, L_m);

% Eigenvalue decomposition (Eq. 5.3)
[T, Lambda] = eig(Y_pole);
L = inv(T);

% Find the eigenvalue closest to zero
[~, idx_critical] = min(abs(diag(Lambda)));
lambda_1 = Lambda(idx_critical, idx_critical);

% Extract corresponding eigenvectors
L_1 = L(idx_critical, :).';  % Left eigenvector (column)
T_1 = T(:, idx_critical);     % Right eigenvector (column)

fprintf('\n=== Critical Eigenvalue ===\n');
fprintf('lambda_1 = %.6f + j%.6f\n', real(lambda_1), imag(lambda_1));

% Calculate generator impedances and conductances at s_pole
Z_gen = @(s) impedance_IG(s, R_s, L_s, R_r, L_r, L_m);
Z_IG1 = Z_gen(s_pole);
Z_IG2 = Z_gen(s_pole);
Z_IG3 = Z_gen(s_pole);

G_IG1 = real(1/Z_IG1);  % Conductance
G_IG2 = real(1/Z_IG2);
G_IG3 = real(1/Z_IG3);

% Generator buses are 2, 3, 4
bus_IG1 = 2;
bus_IG2 = 3;
bus_IG3 = 4;

% Calculate generator participation factors (Eq. 5.5)
P_gen_IG1 = -G_IG1 * abs(L_1(bus_IG1))^2;
P_gen_IG2 = -G_IG2 * abs(L_1(bus_IG2))^2;
P_gen_IG3 = -G_IG3 * abs(L_1(bus_IG3))^2;

% Normalize
P_gen_normalized = [P_gen_IG1, P_gen_IG2, P_gen_IG3];
P_gen_normalized = P_gen_normalized / max(abs(P_gen_normalized));

fprintf('\n=== Generator Participation Factors ===\n');
fprintf('Generator Conductances at f_pole = %.2f Hz:\n', f_pole);
fprintf('  G_IG1 = %.6f S\n', G_IG1);
fprintf('  G_IG2 = %.6f S\n', G_IG2);
fprintf('  G_IG3 = %.6f S\n', G_IG3);

fprintf('\nNormalized Generator Participation Factors:\n');
fprintf('  IG1 (Bus 2): %.4f %s\n', P_gen_normalized(1), ...
        get_interpretation(P_gen_normalized(1)));
fprintf('  IG2 (Bus 3): %.4f %s\n', P_gen_normalized(2), ...
        get_interpretation(P_gen_normalized(2)));
fprintf('  IG3 (Bus 4): %.4f %s\n', P_gen_normalized(3), ...
        get_interpretation(P_gen_normalized(3)));

%% Step 3: Calculate Resonance Participation Factors (Eq. 5.6 & 5.7)

% Series capacitor (between buses 1 and 2)
Z_cap = 1/(s_pole * C_series);
B_cap = imag(1/Z_cap);
P_res_cap = B_cap * abs(L_1(1) - L_1(2))^2;

% Line 1-2 (between buses 1 and 2, excluding capacitor)
Z_line12 = R_line12 + s_pole*L_line12;
B_line12 = imag(1/Z_line12);
P_res_line12 = B_line12 * abs(L_1(1) - L_1(2))^2;

% Line 2-3
Z_line23 = R_line23 + s_pole*L_line23;
B_line23 = imag(1/Z_line23);
P_res_line23 = B_line23 * abs(L_1(2) - L_1(3))^2;

% Line 2-4
Z_line24 = R_line24 + s_pole*L_line24;
B_line24 = imag(1/Z_line24);
P_res_line24 = B_line24 * abs(L_1(2) - L_1(4))^2;

% Induction generators (shunt components at buses 2, 3, 4)
B_IG1 = imag(1/Z_IG1);
B_IG2 = imag(1/Z_IG2);
B_IG3 = imag(1/Z_IG3);
P_res_IG1 = B_IG1 * abs(L_1(bus_IG1))^2;
P_res_IG2 = B_IG2 * abs(L_1(bus_IG2))^2;
P_res_IG3 = B_IG3 * abs(L_1(bus_IG3))^2;

% Normalize
P_res_all = [P_res_cap, P_res_line12, P_res_line23, P_res_line24, ...
             P_res_IG1, P_res_IG2, P_res_IG3];
P_res_normalized = P_res_all / max(abs(P_res_all));

fprintf('\n=== Resonance Participation Factors ===\n');
fprintf('Component               P_res (normalized)  Type\n');
fprintf('Series Capacitor:       %8.4f           %s\n', P_res_normalized(1), ...
        get_resonance_type(P_res_normalized(1)));
fprintf('Line 1-2 (inductor):    %8.4f           %s\n', P_res_normalized(2), ...
        get_resonance_type(P_res_normalized(2)));
fprintf('Line 2-3:               %8.4f           %s\n', P_res_normalized(3), ...
        get_resonance_type(P_res_normalized(3)));
fprintf('Line 2-4:               %8.4f           %s\n', P_res_normalized(4), ...
        get_resonance_type(P_res_normalized(4)));
fprintf('IG1:                    %8.4f           %s\n', P_res_normalized(5), ...
        get_resonance_type(P_res_normalized(5)));
fprintf('IG2:                    %8.4f           %s\n', P_res_normalized(6), ...
        get_resonance_type(P_res_normalized(6)));
fprintf('IG3:                    %8.4f           %s\n', P_res_normalized(7), ...
        get_resonance_type(P_res_normalized(7)));

%% Visualization of Participation Factors

figure('Position', [100 600 1200 400]);

subplot(1,2,1);
bar([1 2 3], P_gen_normalized);
xlabel('Generator Number');
ylabel('Normalized Participation Factor');
title('Generator Participation Factors');
grid on;
set(gca, 'XTickLabel', {'IG1', 'IG2', 'IG3'});
yline(0, 'r--', 'LineWidth', 2);

subplot(1,2,2);
components = {'Cap', 'L12', 'L23', 'L24', 'IG1', 'IG2', 'IG3'};
bar(P_res_normalized);
xlabel('Component');
ylabel('Normalized Participation Factor');
title('Resonance Participation Factors');
grid on;
set(gca, 'XTickLabel', components);
yline(0, 'r--', 'LineWidth', 2);
legend('P_{res}', 'Zero line');

%% Additional Analysis: Bus Participation Factors

P_bus = abs(L_1).^2;
P_bus_normalized = P_bus / max(P_bus);

fprintf('\n=== Bus Participation Factors ===\n');
for i = 1:length(P_bus)
    fprintf('Bus %d: %.4f\n', i, P_bus_normalized(i));
end

%% Helper Functions

function Y = construct_Y_matrix(s, R_grid, L_grid, C_series, ...
                                R_line12, L_line12, R_line23, L_line23, ...
                                R_line24, L_line24, R_s, L_s, R_r, L_r, L_m)
    % Construct 4-bus Y matrix at complex frequency s
    
    % Component impedances
    Z_grid = R_grid + s*L_grid;
    Z_cap = 1/(s*C_series);
    Z_line12 = R_line12 + s*L_line12 + Z_cap;  % Including series cap
    Z_line23 = R_line23 + s*L_line23;
    Z_line24 = R_line24 + s*L_line24;
    
    % Generator impedances (simplified IG model)
    Z_gen = impedance_IG(s, R_s, L_s, R_r, L_r, L_m);
    
    % Initialize Y matrix (4x4)
    Y = zeros(4,4);
    
    % Add admittances
    % Grid connection (bus 1 to ground)
    Y(1,1) = Y(1,1) + 1/Z_grid;
    
    % Line 1-2
    Y(1,1) = Y(1,1) + 1/Z_line12;
    Y(2,2) = Y(2,2) + 1/Z_line12;
    Y(1,2) = Y(1,2) - 1/Z_line12;
    Y(2,1) = Y(2,1) - 1/Z_line12;
    
    % Line 2-3
    Y(2,2) = Y(2,2) + 1/Z_line23;
    Y(3,3) = Y(3,3) + 1/Z_line23;
    Y(2,3) = Y(2,3) - 1/Z_line23;
    Y(3,2) = Y(3,2) - 1/Z_line23;
    
    % Line 2-4
    Y(2,2) = Y(2,2) + 1/Z_line24;
    Y(4,4) = Y(4,4) + 1/Z_line24;
    Y(2,4) = Y(2,4) - 1/Z_line24;
    Y(4,2) = Y(4,2) - 1/Z_line24;
    
    % Generators (shunt at buses 2, 3, 4)
    Y(2,2) = Y(2,2) + 1/Z_gen;
    Y(3,3) = Y(3,3) + 1/Z_gen;
    Y(4,4) = Y(4,4) + 1/Z_gen;
end

function Z = impedance_IG(s, R_s, L_s, R_r, L_r, L_m)
    % Simplified induction generator impedance model
    % This is a basic model; more sophisticated models can be used
    
    Z_s = R_s + s*L_s;
    Z_r = R_r + s*L_r;
    Z_m = s*L_m;
    
    % Equivalent impedance
    Z = Z_s + (Z_m * Z_r)/(Z_m + Z_r);
end

function str = get_interpretation(P_gen)
    if P_gen > 0
        str = '(Generates interharmonic power - SOURCE)';
    else
        str = '(Consumes interharmonic power - DAMPING)';
    end
end

function str = get_resonance_type(P_res)
    if P_res > 0
        str = 'Capacitive resonator';
    else
        str = 'Inductive resonator';
    end
end