function v = simulate_NOM(dt, T, RI, pd, plotHist)
vRest = -70;
vPeak = 30;
vTh = -45;
tau_m = 5;
% initialization
v = zeros(1, T);
v(1) = vRest;
thr = zeros(1, T);
% 2nd order Runge-kutta
i = 1;
while i < T
    thr(i+1) = random(pd) + vTh; % determine the threshold
    k1 = dt * 1/tau_m * (-v(i)+RI(i));
    k2 = dt * 1/tau_m * (-(v(i)+k1)+RI(i));
    v(i+1) = v(i) + (k1+k2)/2;
    i = i + 1;
    if v(i) > thr(i)
        v(i) = vPeak;
        v(i+1) = vRest;
        i = i + 1;
    end
end
if plotHist
    histogram(thr(v==vPeak))
end





























