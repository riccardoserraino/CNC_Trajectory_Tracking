clear; clc; close all

%% INUPT

load('.mat')    % modify the name of data obj to load for continuing the analysis
Ak = 0.5;
rng(29);

f = [5 8 12 17 20 23 27 32 35 40];
phi = 2*pi*rand(1,10);

%% MODEL IDENTIFICATION
% -- PREPARE EXPERIMENTAL DATA
% Convert data to N x 3 if necessary

if size(data,1) == 3 && size(data,2) > 3
    data = data.';          % was 3 x N -> N x 3
end

t_exp = data(:,1);
u_exp = data(:,2);
y_exp = data(:,3);

Ts = 0.0020;

% -- CREATE IDENTIFICATION DATA

data_idd = iddata(y_exp, u_exp, Ts);
N = size(data_idd, 1);      % number of samples
fprintf('Campioni totali: %d\n', N);
assert(N > 10, ...
    'too few samples: check data variable');

% -- TRAINING / VALIDATION SPLIT

Ntrain = floor(0.8*N);
data_train = data_idd(1:Ntrain);
data_val   = data_idd(Ntrain+1:end);

% -- SYSTEM IDENTIFICATION

sys = tfest(data_train, 3, 0);

% -- MODEL VALIDATION

[yhat, fit] = compare(data_val, sys);
fprintf('\nQualità fitting: %.2f%%\n', fit);
sys
bode(sys)


%% MECHANICAL PARAMETER ESTIMATION
% Closing the back-EMF loop and integrating gives:
%   theta(s)/V(s) = (Kt*N*conv/(Lm*Jeq)) / (s^3 + a2*s^2 + a1*s)
%   a2 = (Lm*beta + Rm*Jeq)/(Lm*Jeq)
%   a1 = (Rm*beta + Kt^2*N^2)/(Lm*Jeq)
% This is why tfest(data_train,3,0) is used. Lm, Rm, Kt, N are KNOWN
% (electrical part not estimated); we solve the two equations above for
% the two mechanical unknowns Jeq, beta.

Lm   = Lm;              % known electrical inductance [H]
conv = 360/(2*pi);      % 1 if y_exp is in rad, 360/(2*pi) if in degrees

[num, den] = tfdata(sys, 'v');   % den = [1 a2 a1 a0], num = [b0]
a2 = den(2);
a1 = den(3);

Jeq_est  = Lm*Kt^2*N^2 / (Rm^2 + a1*Lm^2 - a2*Lm*Rm);
beta_est = (a1*Lm*Jeq_est - Kt^2*N^2) / Rm;

%% ================= RESULTS =================
fprintf('\n--- ESTIMATED MECHANICAL PARAMETERS ---\n');
fprintf('Estimated beta : %.6f Nms/rad\n', beta_est);
fprintf('Estimated Jeq  : %.6e kg m^2\n', Jeq_est);

fprintf('\n--- COMPARISON ---\n');
fprintf('J original     : %.6e kg m^2\n', Jtot);
fprintf('J estimated    : %.6e kg m^2\n', Jeq_est);
fprintf('beta original  : %.6f Nms/rad\n', beta);
fprintf('beta estimated : %.6f Nms/rad\n', beta_est);
