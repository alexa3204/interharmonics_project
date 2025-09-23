close all; clear all; clc 

% out = sim('three_bus_fundamental.slx');

% % voltage at bus 3 for phases abc 
% V3a = out.v3_data(:,1);
% V3b = out.v3_data(:,2);
% V3c = out.v3_data(:,3);

fs = 20000; % sampling frequency 
duration = 1; % test duration (s) 
t = 0:1/fs:duration; 
t = t';

% Fundamental frequency
f1 = 60;           % Hz (North American grid frequency)
omega1 = 2*pi*f1;  % rad/s

% voltage at bus 3 for phases abc 
V3a = cos(omega1 * t);
V3b = cos(omega1 * t - 2*pi/3);
V3c = cos(omega1 * t - 4*pi/3);

% % transpose to column vectors 
% V3a = V3a'; 
% V3b = V3b';
% V3c = V3c'; 

figure
plot(t, V3a); 
hold on; 
plot(t, V3b); 
plot(t, V3c);
legend('Va','Vb', 'Vc');

v3_ps = positive_sequence_waveform(V3a, V3b, V3c, 60, t); % positive sequence phasor 

figure 
plot(t,[v3_ps, V3a]);

%% test 3-bus system using just phase A 

%% Line admittances
y12 = 1/(0.01 + 1j*0.085); 
y23 = 1/(0.02 + 1j*0.161); 
y13 = 1/(0.01 + 1j*0.092); 


% Shunt admittances
y12_sh = 1j*0.088; 
y23_sh = 1j*0.153; 
y13_sh = 1j*0.079; 

% Add each shunt admittance to the corresponding bus diagonals
Y_r = [y12 + y13 + y12_sh + y13_sh, -y12, -y13; ...
       -y12, y12 + y23 + y12_sh + y23_sh, -y23; ...
       -y13, -y23, y13 + y23 + y13_sh + y23_sh];
 
