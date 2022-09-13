function [v, u] = simulate_Izh(a, b, c, d, dt, I, T, vRest)
v = zeros(1, T); % output voltage
v(1) = vRest;
u = zeros(1, T);
u(1) = b*v(1);
% 2nd order Runge-Kutta
k1 = zeros(2, 1);
k2 = zeros(2, 1);
for i = 1:T-1
    k1(1) = dt * (0.04*v(i)^2 + 5*v(i) + 140 - u(i) + I(i));
    k1(2) = dt * a * (b*v(i) - u(i));
    k2(1) = dt * (0.04*(v(i)+k1(1))^2 + 5*(v(i)+k1(1)) + 140 - (u(i)+k1(2)) + I(i));
    k2(2) = dt * a * (b*(v(i)+k1(1)) - (u(i)+k1(2)));

    v(i+1) = v(i) + (k1(1)+k2(1))/2;
    u(i+1) = u(i) + (k1(2)+k2(2))/2;

    if v(i+1) > 30
        v(i+1) = c;
        u(i+1) = u(i+1) + d;
    end
end






























