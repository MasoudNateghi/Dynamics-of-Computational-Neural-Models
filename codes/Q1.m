%% Q1.1
clear; close all; clc;
vRest = -60;
v = -80:0.01:10;
u = vRest - v;
alpha_n = (.1 * u + 1)./(exp(1 + .1 * u) - 1) / 10;
beta_n = .125 * exp(u/80);
alpha_m = (u+25) ./ (exp(2.5+.1*u)-1)/10;
beta_m = 4*exp(u/18);
alpha_h = .07 * exp(u/20);
beta_h = 1 ./ (1+exp(3 + .1*u));
% time constants
tau_n = 1 ./ (alpha_n + beta_n);
tau_m = 1 ./ (alpha_m + beta_m);
tau_h = 1 ./ (alpha_h + beta_h);

figure;
hold on; grid minor;
plot(v, tau_h, 'LineWidth', 2)
plot(v, tau_m, 'LineWidth', 2)
plot(v, tau_n, 'LineWidth', 2)
title('time constants', 'Interpreter','latex')
legend('\tau_h', '\tau_m', '\tau_n')

% steady state values
n_inf = alpha_n ./ (alpha_n + beta_n);
m_inf = alpha_m ./ (alpha_m + beta_m);
h_inf = alpha_h ./ (alpha_h + beta_h);

figure;
hold on; grid minor;
plot(v, h_inf, 'LineWidth', 2)
plot(v, m_inf, 'LineWidth', 2)
plot(v, n_inf, 'LineWidth', 2)
title('steady state values', 'Interpreter','latex')
legend('$$h_{\infty}$$', '$$m_{\infty}$$', '$$n_{\infty}$$', 'Interpreter', 'latex')
%% Q2.1
dt = 0.01; % Simulation time step
Duration = 200; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
I = zeros(1, T); % in uA, external stimulus (external current)
I(1:T) = 200; % an input current pulse
[v, ~, ~, ~] = simulate_HH(dt, I, T);
figure;
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('v(t)', 'Interpreter','latex')
%% Q5.1
% !!!!!!!!!!!!!!!!! LONG RUNTIME !!!!!!!!!!!!!!!!!
dt = 0.01; % Simulation time step
Duration = 2000; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; %#ok<NASGU> % Simulation time points in ms
external_input = linspace(6, 112, 100);
amp = zeros(size(external_input)); % amplitude
FR  = zeros(size(external_input)); % firing rate

for count = 1:length(external_input)
    count %#ok<NOPTS> 
    I = external_input(count) * ones(1, T); % in uA, external stimulus (external current)
    [v, ~, ~, ~] = simulate_HH(dt, I, T);
    [pks,locs] = findpeaks(v);
    amp(count) = mean(pks(pks>=-20));
    FR(count) = sum(pks>-20) / (Duration*1e-3);
end

figure;
plot(external_input, amp)
xlabel('External Current($$\mu$$A)', 'Interpreter','latex')
ylabel('Amplitude(mv)', 'Interpreter','latex')
title('Amplitude vs. External Current', 'Interpreter','latex')
figure;
plot(external_input, FR)
xlabel('External Current($$\mu$$A)', 'Interpreter','latex')
ylabel('Firing Rate(Hz)', 'Interpreter','latex')
title('Firing Rate vs. External Current', 'Interpreter','latex')
%% Q7.1
dt = 0.01; % Simulation time step
Duration = 2000; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
I = linspace(1, 200, T);
[v, ~, ~, ~] = simulate_HH(dt, I, T);
figure;
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('v(t)', 'Interpreter','latex')
% figure;
% plot(t, I)
% xlabel('Time(ms)', 'Interpreter','latex')
% ylabel('I(t)', 'Interpreter','latex')
%% Q8.1
% very low input current
dt = 0.01; % Simulation time step
Duration = 2000; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; %#ok<NASGU> % Simulation time points in ms
I = 2*ones(1, T);
[v, n, m, h] = simulate_HH(dt, I, T);

figure;
subplot(1, 3, 1)
plot(v, n); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('n', 'Interpreter','latex')
subplot(1, 3, 2)
plot(v, m); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('m', 'Interpreter','latex')
subplot(1, 3, 3)
plot(v, h); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('h', 'Interpreter','latex')
sgtitle('$$I = 2\mu A$$', 'Interpreter', 'latex')

% normal input current
I = 8*ones(1, T);
[v, n, m, h] = simulate_HH(dt, I, T);

figure;
subplot(1, 3, 1)
plot(v, n); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('n', 'Interpreter','latex')
subplot(1, 3, 2)
plot(v, m); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('m', 'Interpreter','latex')
subplot(1, 3, 3)
plot(v, h); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('h', 'Interpreter','latex')
sgtitle('$$I = 8\mu A$$', 'Interpreter', 'latex')

% very high input current
I = 200*ones(1, T);
[v, n, m, h] = simulate_HH(dt, I, T);

figure;
subplot(1, 3, 1)
plot(v, n); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('n', 'Interpreter','latex')
subplot(1, 3, 2)
plot(v, m); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('m', 'Interpreter','latex')
subplot(1, 3, 3)
plot(v, h); grid minor;
xlabel('voltage(mv)', 'Interpreter','latex')
ylabel('h', 'Interpreter','latex')
sgtitle('$$I = 200\mu A$$', 'Interpreter', 'latex')
%% Q9.1
dt = 0.001; % Simulation time step
Duration = 2000; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
% triangular input
f = 10;
I = 10*sawtooth(2*pi*f*t*1e-3,1/2);
[v, ~, ~, ~] = simulate_HH(dt, I, T);
figure;
subplot(2, 1, 1)
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('v(t)', 'Interpreter','latex')
subplot(2, 1, 2)
plot(t, I)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('External Current($$\mu$$A)', 'Interpreter','latex')
%% sinusoidal input
f = 25;
I = 10*sin(2*pi*f*t*1e-3);
% I = 10*ones(1, T);
[v, ~, ~, ~] = simulate_HH(dt, I, T);
figure;
subplot(2, 1, 1)
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('v(t)', 'Interpreter','latex')
subplot(2, 1, 2)
plot(t, I)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('External Current($$\mu$$A)', 'Interpreter','latex')
str = strcat('f = ', num2str(f), '(Hz)');
sgtitle(str, 'Interpreter', 'latex')
%% square input
f = 10;
I = 10*square(2*pi*f*t*1e-3,50);
[v, ~, ~, ~] = simulate_HH(dt, I, T);
figure;
subplot(2, 1, 1)
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('v(t)', 'Interpreter','latex')
subplot(2, 1, 2)
plot(t, I)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('External Current($$\mu$$A)', 'Interpreter','latex')
%% chirp input
f0 = 10; f1 = 40; t1 = 1;
I = 10*chirp(t*1e-3, f0, t1, f1);
[v, ~, ~, ~] = simulate_HH(dt, I, T);
figure;
subplot(2, 1, 1)
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('v(t)', 'Interpreter','latex')
subplot(2, 1, 2)
plot(t, I)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('External Current($$\mu$$A)', 'Interpreter','latex')




































