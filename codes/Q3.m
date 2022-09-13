clear; close all; clc;
%% Q1.3
gamma = 1;
pd = makedist("Logistic", "mu", 0, "sigma", 1/(2*gamma));
x = -15:0.01:15;
y_pdf = pdf(pd, x);
y_cdf = cdf(pd, x);
figure;
plot(x, y_cdf); grid minor;
xlabel('x', 'Interpreter','latex')
ylabel('$$F_X(x)$$', 'Interpreter','latex')
str = "cdf of logistic distribution with $$\gamma = $$" + num2str(gamma);
title(str, 'Interpreter','latex')
figure;
plot(x, y_pdf); grid minor;
xlabel('x', 'Interpreter','latex')
ylabel('$$f_X(x)$$', 'Interpreter','latex')
str = "pdf of logistic distribution with $$\gamma = $$" + num2str(gamma);
title(str, 'Interpreter','latex')
%% Q2.3
dt = 0.01; % Simulation time step
Duration = 20; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
RI = 120*ones(1,T); %mV
v = simulate_NOM(dt, T, RI, pd, false);
figure;
plot(t, v)
xlabel('Time(ms)', 'Interpreter','latex')
ylabel('Voltage(mv)', 'Interpreter','latex')
str = "Noisy Output Model with $$\gamma= $$" + num2str(gamma) + " and RI = " + num2str(RI(1));  
title(str, 'Interpreter','latex')
%% Q3.3
dt = 0.01; % Simulation time step
Duration = 100; % Simulation length
T = ceil(Duration/dt);
% t = (1:T) * dt; % Simulation time points in ms
RI = 20*ones(1,T); %mV
simulate_NOM(dt, T, RI, pd, true);
%% Q4.3
% !!!!!!!!!!!!!!!!!!!!! LONG RUNTIME !!!!!!!!!!!!!!!!!!!!!
dt = 0.01; % Simulation time step
Duration = 100; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
R = 1; %mOhm
I = linspace(-70, 100, 100);
F = zeros(size(I));
for i =  1:length(I)
    i %#ok<NOPTS> 
    RI = R*I(i)*ones(1, T);
    v = simulate_NOM(dt, T, RI, pd, false);
    [pks, ~] = findpeaks(v);
    F(i) = numel(pks) / (Duration*1e-3);
end
plot(I, F)
xlabel('external current($$\mu A$$)', "Interpreter","latex")
ylabel("Firing Rate(Hz)", "Interpreter","latex")
title("F-I plot", "Interpreter","latex")





















