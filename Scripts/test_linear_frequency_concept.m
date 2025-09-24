%% Test phasor extraction on system with harmonics 
% 60hz with current injection at 63hz 

close all; clear all; clc 

out = sim('three_bus_current_inj_bal.slx');

fs = 20000; % sampling frequency 
duration = 1; % test duration (s) 
t = 0:1/fs:duration; 

% Frequencies (fundamental & interharmonic)
f = [60 63];      % Hz
omega1 = 2*pi*f;  % rad/s

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

    [~, idx1] = min(abs(f_axis - f(1))); % Index for 60Hz
    [~, idx2] = min(abs(f_axis - f(2))); % Index for 63Hz

    % Extract phasors from FFT 
    I_r(i,1) = I_fft(i,idx1) * 2/N; 
    I_r(i,2) = I_fft(i,idx2) * 2/N; 
end 


fprintf('Current Phasors:\n');
[abs(I_r) rad2deg(angle(I_r))]

% Voltage phasors
V = zeros(3, length(t)); 
V(1,:) = out.v1_data(:,1); 
V(2,:) = out.v2_data(:,1); 
V(3,:) = out.v3_data(:,1); 

for i = 1:3 
    V_fft(i,:) = fft(V(i,:));

    [~, idx1] = min(abs(f_axis - f(1)));
    [~, idx2] = min(abs(f_axis - f(2)));

    V_r(i,1) = V_fft(i,idx1) * 2/N; 
    V_r(i,2) = V_fft(i,idx2) * 2/N; 
end 

fprintf('Voltage Phasors:\n');
[abs(V_r) rad2deg(angle(V_r))]

% Y-matrix for 60Hz 
w = 2*pi*f(1); 

% Define base values
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

Y60 = Ybus;

% Y-matrix for 63Hz 
w = 2*pi*f(2); 

% Define base values
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

Y63 = Ybus;

% Compare I = YV
fprintf('60 Hz: \n')
Y60
fprintf('Y_60*V_60:')
Y60*V_r(:,1)
fprintf('I_60:')
I_r(:,1)

fprintf('63 Hz: \n')
Y63
fprintf('Y_63*V_63:')
Y63*V_r(:,2)
fprintf('I_63:')
I_r(:,2)