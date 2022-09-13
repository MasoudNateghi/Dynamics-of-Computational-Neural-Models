clear; close all; clc;
EL = -70e-3;
Es_exc = 0e-3;
Es_inh = -80e-3;
tau_m = 20e-3;
IeRm = 25e-3;
rm = 100e-3; 
tau_peak = 10e-3;
vRest = -80e-3;
vTh = -54e-3;
vPeak = 0e-3;
%% single neuron
dt = 0.01e-3; % Simulation time step
Duration = 1; % Simulation length
T = ceil(Duration/dt);
t = (0:T) * dt; % Simulation time points in ms
rho = zeros(size(t)); % construct train impulse signal
Tau = 10e-3; % period of pulse train
% create pulse train
idx = int32((Tau/dt:Tau/dt:T)+1);
rho(idx) = 1; 
K = t/tau_peak .* exp(1-(t/tau_peak)); % K(t)
g = conv(rho, K); % g(t)
% excitatory
v_exc = zeros(size(t));
v_exc(1) = -60e-3;
i = 1;
while i < length(t)
    v_exc(i+1) = v_exc(i) + dt/tau_m*(-(v_exc(i)-EL) -g(i)*(v_exc(i)-Es_exc)*rm + IeRm);
    if v_exc(i+1) > vTh
        v_exc(i+1) = vPeak;
        v_exc(i+2) = vRest;
        i = i + 1;
    end
    i = i + 1;
end
% inhibitory
v_inh = zeros(size(t));
v_inh(1) = -60e-3;
i = 1;
while i < length(t)
    v_inh(i+1) = v_inh(i) + dt/tau_m*(-(v_inh(i)-EL) -g(i)*(v_inh(i)-Es_inh)*rm + IeRm);
    if v_inh(i+1) > vTh
        v_inh(i+1) = vPeak;
        v_inh(i+2) = vRest;
        i = i + 1;
    end
    i = i + 1;
end

figure;
subplot(3, 1, 1)
plot(t*1000, v_exc*1000)
grid minor
title('excitatory (single neuron)', "Interpreter","latex")
xlabel("Time(ms)", "Interpreter","latex")
ylabel("Voltage(mV)", "Interpreter","latex")
subplot(3, 1, 2)
plot(t*1000, v_inh*1000)
grid minor
title('inhibitory (single neuron)', "Interpreter","latex")
xlabel("Time(ms)", "Interpreter","latex")
ylabel("Voltage(mV)", "Interpreter","latex")
subplot(3, 1, 3)
plot(t*1000, rho)
grid minor
title('$$\rho(t)$$', "Interpreter","latex")
% str = "period of palse train = " + num2str(Tau*1000) + "ms";
str = "$$V_0 = $$" + num2str(v_exc(1)*1000) + "mV";
sgtitle(str, "interpreter", "latex")
%% connected neurons
dt = 0.01e-3; % Simulation time step
Duration = 0.7; % Simulation length
T = ceil(Duration/dt);
t = (0:T) * dt; % Simulation time points in ms
K = t/tau_peak .* exp(1-(t/tau_peak)); % K(t)
% excitatory
rho1 = zeros(size(t));
rho2 = zeros(size(t));
v1 = zeros(size(t));
v1(1) = -60e-3;
v2 = zeros(size(t));
v2(1) = -80e-3;
g1 = zeros(size(t));
g2 = zeros(size(t));
spikeFlag1 = 0; 
spikeFlag2 = 0;
for i = 1:length(t)-1
    if ~spikeFlag1
        v1(i+1) = v1(i) + dt/tau_m*(-(v1(i)-EL) -g1(i)*(v1(i)-Es_exc)*rm + IeRm);
    else
        v1(i+1) = vRest;
        spikeFlag1 = 0; 
    end
    if ~spikeFlag2
        v2(i+1) = v2(i) + dt/tau_m*(-(v2(i)-EL) -g2(i)*(v2(i)-Es_exc)*rm + IeRm);
    else
        v2(i+1) = vRest;
        spikeFlag2 = 0; 
    end
    if v1(i+1) > vTh
        v1(i+1) = 0;
        spikeFlag1 = 1;
        rho2(i+1) = 1;
        g2 = conv(K, rho2);
    end
    if v2(i+1) > vTh
        v2(i+1) = 0;
        spikeFlag2 = 1;
        rho1(i+1) = 1;
        g1 = conv(K, rho1);
    end
end
v1_exc = v1;
v2_exc = v2;
% inhibitory
rho1 = zeros(size(t));
rho2 = zeros(size(t));
v1 = zeros(size(t));
v1(1) = -60e-3;
v2 = zeros(size(t));
v2(1) = -80e-3;
g1 = zeros(size(t));
g2 = zeros(size(t));
spikeFlag1 = 0; 
spikeFlag2 = 0;
for i = 1:length(t)-1
    if ~spikeFlag1
        v1(i+1) = v1(i) + dt/tau_m*(-(v1(i)-EL) -g1(i)*(v1(i)-Es_inh)*rm + IeRm);
    else
        v1(i+1) = vRest;
        spikeFlag1 = 0; 
    end
    if ~spikeFlag2
        v2(i+1) = v2(i) + dt/tau_m*(-(v2(i)-EL) -g2(i)*(v2(i)-Es_inh)*rm + IeRm);
    else
        v2(i+1) = vRest;
        spikeFlag2 = 0; 
    end
    if v1(i+1) > vTh
        v1(i+1) = 0;
        spikeFlag1 = 1;
        rho2(i+1) = 1;
        g2 = conv(K, rho2);
    end
    if v2(i+1) > vTh
        v2(i+1) = 0;
        spikeFlag2 = 1;
        rho1(i+1) = 1;
        g1 = conv(K, rho1);
    end
end
v1_inh = v1;
v2_inh = v2;
figure;
subplot(2, 1, 1)
plot(t*1e3, v1_exc*1e3)
hold on
plot(t*1e3, v2_exc*1e3)
grid minor
title('excitatory (connected neurons)', "Interpreter","latex")
xlabel("Time(ms)", "Interpreter","latex")
ylabel("Voltage(mV)", "Interpreter","latex")
subplot(2, 1, 2)
plot(t*1e3, v1_inh*1e3)
hold on
plot(t*1e3, v2_inh*1e3)
grid minor
title('inhibitory (connected neurons)', "Interpreter","latex")
xlabel("Time(ms)", "Interpreter","latex")
ylabel("Voltage(mV)", "Interpreter","latex")


















