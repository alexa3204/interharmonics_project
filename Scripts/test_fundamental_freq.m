%% test fundamental frequency I=YV from phase a 
% logic is that in a balanced system, phase a = pos-sequence component 
close all; clear all; clc 

out = sim('three_bus_fundamental.slx');

fs = 20000; % sampling frequency 
duration = 1; % test duration (s) 
t = 0:1/fs:duration; 

% Fundamental frequency
f1 = 60;           % Hz
omega1 = 2*pi*f1;  % rad/s

% get current vectors from phase a at all buses 
i = zeros(3, length(t)); 
i(1,:) = out.i1_data(:,1); 
i(2,:) = out.i2_data(:,1); 
i(3,:) = out.i3_data(:,1); 

% get voltage vectors from phase a at all buses 
v = zeros(3, length(t)); 
v(1,:) = out.v1_data(:,1); 
v(2,:) = out.v2_data(:,1); 
v(3,:) = out.v3_data(:,1); 

% voltage phasor 
V_ph = ones*length()

% Line admittances
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

%% Test phasor extraction 
