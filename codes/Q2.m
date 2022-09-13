clear; close all; clc;
dt = 0.01; % Simulation time step
Duration = 100; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
%% Q1.3.2
a = 0.02; b = 0.2; c = -65; d = 2; h = 15; vRest = -70;
I = zeros(1, T);
I(t >= 10) = h;
v = simulate_Izh(a, b, c, d, dt, I, T, vRest);
figure;
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('Voltage(mv)', 'Interpreter','latex')
title('Tonic Spiking', 'Interpreter','latex')
%% Q2.3.2/4.2
a = 0.02; b = 0.25; c = -65; d = 6; h = 1; vRest = -54.3;
I = zeros(1, T);
I(t >= 10) = h;
[v, u] = simulate_Izh(a, b, c, d, dt, I, T, vRest);
figure;
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('Voltage(mv)', 'Interpreter','latex')
title('Phasic Spiking', 'Interpreter','latex')
figure;
plot(v, u)
xlabel('v', 'Interpreter','latex')
ylabel('u', 'Interpreter','latex')
title('v-u figure', 'Interpreter','latex')
%% Q3.3.2
a = 0.02; b = 0.2; c = -50; d = 2; h = 15; vRest = -70;
I = zeros(1, T);
I(t >= 10) = h;
v = simulate_Izh(a, b, c, d, dt, I, T, vRest);
figure;
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('Voltage(mv)', 'Interpreter','latex')
title('Tonic Bursting', 'Interpreter','latex')
%% Q4.3.2
a = 0.02; b = 0.25; c = -55; d = 0.05; h = 0.6; vRest = -64.4;
I = zeros(1, T);
I(t >= 10) = h;
v = simulate_Izh(a, b, c, d, dt, I, T, vRest);
figure;
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('Voltage(mv)', 'Interpreter','latex')
title('Phasic Bursting', 'Interpreter','latex')
%% Q5.3.2
a = 0.02; b = 0.2; c = -55; d = 4; h = 10; vRest = -70;
I = zeros(1, T);
I(t >= 10) = h;
v = simulate_Izh(a, b, c, d, dt, I, T, vRest);
figure;
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('Voltage(mv)', 'Interpreter','latex')
title('Phasic Bursting', 'Interpreter','latex')


























