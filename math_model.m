clear; clc; close all

%% PARAMETERS 
Rm   = 2.6;
Lm   = 0.18e-3;
Km   = 7.68e-3;
Kt   = 7.68e-3;
N    = 70;
Jm   = 3.90e-7;
Jeq  = 2.087e-3;
beta = 0.015;
Jtot = Jm + Jeq;
conv = 360/(2*pi); % conversion to degrees

%% TRANSFER FUNCTIONS
s = tf('s');

% Full model (with Lm)
G_elec       = 1/(Lm*s + Rm);
G_mech       = 1/(Jtot*s + beta);
G_omega_full = feedback(N*Kt*G_elec*G_mech, N*Km);
G_theta_full = conv * G_omega_full/s;                 % integrator adds pole at 0

% Simplified model (Lm neglected)
G_simple       = N*Kt/(Rm*(Jtot*s + beta));
G_omega_simple = feedback(G_simple, N*Km);
G_theta_simple = conv * G_omega_simple/s;

%% CALCULATIONS
% --- (Lm neglected) ---
c1 = Rm*Jtot;
c0 = beta*Rm + Kt*Km*N^2;

w_open   = beta/Jtot;               % open loop, no back-EMF   [rad/s]
w_closed = c0/c1;                   % closed loop, den = 0     [rad/s]
w_unity  = (Kt*N - c0)/c1;          % G(s) = 1                 [rad/s]
w_elec   = Rm/Lm;                   % electrical pole          [rad/s]

f_open   = w_open/(2*pi);
f_closed = w_closed/(2*pi);
f_unity  = w_unity/(2*pi);
f_elec   = w_elec/(2*pi);

tau_open   = 1/w_open;
tau_closed = 1/w_closed;
tau_elec   = 1/w_elec;

% --- Poles from the transfer functions, sorted slow to fast ---
p_theta_simple = pole(G_theta_simple);
p_theta_full   = pole(G_theta_full);
[~, i1] = sort(abs(p_theta_simple));  p_theta_simple = p_theta_simple(i1);
[~, i2] = sort(abs(p_theta_full));    p_theta_full   = p_theta_full(i2);

f_theta_simple = abs(p_theta_simple)/(2*pi);
f_theta_full   = abs(p_theta_full)/(2*pi);

%% PRINT
fprintf('\n--- MECHANICAL (Lm neglected) ---\n')
fprintf('Open-loop pole  : %8.2f rad/s = %6.2f Hz | tau = %.4f s\n', w_open,   f_open,   tau_open)
fprintf('Closed-loop pole: %8.2f rad/s = %6.2f Hz | tau = %.4f s\n', w_closed, f_closed, tau_closed)
fprintf('G(s) = 1        : %8.2f rad/s = %6.2f Hz\n',                w_unity,  f_unity)

fprintf('\n--- ELECTRICAL ---\n')
fprintf('Electrical pole : %8.0f rad/s = %6.1f Hz | tau = %.6f s\n', w_elec, f_elec, tau_elec)

fprintf('\n--- SIMPLIFIED (2 poles), theta/Va POLES ---\n')
for k = 1:numel(p_theta_simple)
    fprintf('p%d = %10.2f rad/s = %8.2f Hz\n', k, abs(p_theta_simple(k)), f_theta_simple(k))
end

fprintf('\n--- FULL (3 poles), theta/Va POLES ---\n')
for k = 1:numel(p_theta_full)
    fprintf('p%d = %10.2f rad/s = %8.2f Hz\n', k, abs(p_theta_full(k)), f_theta_full(k))
end

%% BODE theta/Va: SIMPLIFIED vs FULL
figure
bode(G_theta_simple, 'b', G_theta_full, 'r--', {0.1, 1e5})
grid on
legend('Simplified (Lm neglected)', 'Full (Lm included)', 'Location', 'southwest')
title('Bode: \theta/V_a')

%% CUTOFF FREQUENCIES
% -3 dB bandwidth of omega/Va (theta/Va has an integrator, so no finite DC gain)
wc_simple = bandwidth(G_omega_simple);      % [rad/s]
wc_full   = bandwidth(G_omega_full);        % [rad/s]
fc_simple = wc_simple/(2*pi);               % [Hz]
fc_full   = wc_full/(2*pi);                 % [Hz]

fprintf('\n--- CUTOFF (-3 dB) of omega/Va ---\n')
fprintf('Simplified: %8.2f rad/s = %6.2f Hz\n', wc_simple, fc_simple)
fprintf('Full      : %8.2f rad/s = %6.2f Hz\n', wc_full,   fc_full)


