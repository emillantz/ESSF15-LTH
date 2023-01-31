%% Plot simulink values
out = sim('rl.slx');
i1 = out.i1;
u1 = out.u1;
t = out.tout;

figure(1);
plot(t, i1)
hold on
plot(t, u1)

xlabel('time (t)')
legend('Current', 'Voltage')

% Constants
R = 40; % Ohm
L = 0.095; % Henry
f = 50; % Hz
V = 220*sqrt(2); % V
V_RMS = V/sqrt(2); % V

%% Voltage (1b)
v_top = max(u1(3:end))
v_cpx = v_top % Voltage is leading phase => vtop*exp(0) = vtop
%% Current (1c)
phase_diff = acos(dot(u1, i1) / (norm(u1) * norm(i1))) % rad
i_top = max(i1(3:end))
i_cpx = i_top*exp(-i*phase_diff)

%% Impedance + Resistance + Inductance (1d+1e)
xr = v_cpx / i_cpx;
fprintf('Z = %f, X_L = %f, R = %f \n', abs(xr), imag(xr), real(xr))

%% Power (2a+2b)
P = real(v_cpx * conj(i_cpx));
Q = imag(v_cpx * conj(i_cpx));
S = abs(v_cpx * conj(i_cpx));
EF = P/S;

fprintf(['P = %f W, Q = %f VAr, S = %f VA, effect factor = %f' ...
    '\n\n'], P, Q, S, EF);

figure(2);
plot(t, u1.*i1)
fprintf('Plot top values ~873 W (after first period), mean value = %f \n\n', mean(u1.*i1))

%% Triple-phase
out_3 = sim('trefas.slx');
% 3a => Change the phase for L2 and L3 by 120 and 240 degrees respectively

voltages = [out_3.u1, out_3.u2, out_3.u3]';
currents = [out_3.i1, out_3.i2, out_3.i3]';
t = out_3.tout;
% 3b
figure(3)
hold on
for i = 1:3 
    plot(t, voltages(i,:))
end
xlabel('time (t)')
legend('L1 Voltage', 'L2 Voltage', 'L3 Voltage')

figure(4)
hold on
for i = 1:3 
    plot(t, currents(i,:))
end
xlabel('time (t)')
legend('L1 Current', 'L2 Current', 'L3 Current')

% Plot sum of currents at each index
figure(5)
plot(t, sum(currents))
xlabel('time (t)')
legend('Sum of currents')

% 3c and 3d on paper

% 3e
figure(6)
plot(t, sum(voltages(1:2, :)))
hold on
plot(t, voltages(1, :))
xlabel('time(t)')
legend('Sum of voltages in L1 and L2', 'L1 Voltage')

% 3f
powers = voltages .* currents;
p1 = powers(1, :);
p12 = powers(1, :) + powers(2, :);
p123 = powers(1, :) + powers(2, :) + powers(3, :);

figure(7)
plot(t, p1)
hold on
plot(t, p12)
plot(t, p123)
xlabel('time(t)')
legend('P1', 'P1+P2', 'P1+P2+P3')

%% 4
out_gen = sim("infasning.slx");

figure(8)
hold on
plot(out_gen, uelnat)
plot(out_gen, usg)
plot(out_gen, u_over_brytare)
xlabel('time(t)')
legend('Nätspänning', 'Generator', 'Skillnadsspänning')
