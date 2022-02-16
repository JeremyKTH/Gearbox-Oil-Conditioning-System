clc
clear all

%The ARMAX model is based on the paper "Heat Exchanger Model (Data-driven modelling)"

num = [0.0001701, 0.0001295];          % B
den = [1, -0.2492, -0.7469, 0, 0, 0];  % A
Ts = 1;
sys = tf(num, den, Ts)

%% Pole-Zero Cancellation
%Can use this to simplify the model if there exists pole-zero cancellation
sys = minreal(sys, 1e-3)   % in this case, it does not

%% Checking  the system polse, zeros and DC gain
[ze, p, k] = tf2zp(num, den) %[zero, pole, gain]

%% Finding dominate pole (reducing order)
z = tf('z', Ts);
sys_dominant = k*(z+0.7613)/((z-0.9978)*(z+0.7486))  % remember to include the DC gain (k)

%% Proving reduced-order system is still a good approximation of original
figure(1)
step(sys, sys_dominant)
legend({'original system', 'dominant pole'})

S = stepinfo(sys)
display('--------------------------')
S = stepinfo(sys_dominant)
% From the output below it is clear that sys and sys_dominant are almost identical

%% Converting the poles to continuous time. 
%The rest is just to double check if the converted poles are correct.
sys_cont = d2c(sys_dominant, 'tustin');

pole_cont = pole(sys_cont)
zero_cont = zero(sys_cont)

s = tf('s');
%sys_cont_proof = ((s-0.5)*(s+3.6894))/((s+0.0006)*(s+3.477))    % Ts = 4    , P = -0.0006, -3.477  , Z = 0.5   , -3.6894    
%sys_cont_proof = ((s-2)*(s+14.7574))/((s+0.0022)*(s+13.9109))   % Ts = 1    , P = -0.0022, -13.9109, Z = 2     , -14.7574
sys_cont_proof = ((s+1475.7)*(s+200))/((s+1391.1)*(s+0.2))       % Ts = 0.01 , P = -0.2   , -1391.1 , Z = -2000 , -1475.7
%sys_cont_proof = ((s+14757)*(s+20000))/((s+13911)*(s+2))        % Ts = 0.001, P = -2     , -13911  , Z = -20000, -14757

% check dc gain
% k = DC gain of original system
k_d = dcgain(sys_cont) % k_d : sys_dominant = sys_cont
k_c = dcgain(sys_cont_proof)
sys_cont_proof = sys_cont_proof/k_c

figure(2)
step((sys/k_d), (sys_cont/k_d), sys_cont_proof)
legend({'original system', 'dominant pole', 'proof_cont'})

%% Pole Placement - PI Controller
 
% Retreive coefficients of the plant...
Gp = sys_dominant;
[num, den] = tfdata(Gp);
b1 = num{1}(2);  b0 = num{1}(3);
a1 = den{1}(2);  a0 = den{1}(3);


%--- A_m ------
w_m = 35;  % 10 times faster, origial is at -3.4777
zeta_m = 0.7;

%--- A_o ------
w_o = 80;

%===============================
% Converting to discrete:
%--- A_m ------
p2 = 2*exp(-zeta_m*w_m*Ts)*cos(w_m*Ts*sqrt(1-(zeta_m)^2)); % z1 + z2
p1 = exp(-2*zeta_m*w_m*Ts); % z1*z2

%--- A_o ------
p0 = exp(-w_o*Ts);


% Diophantine Eqn A_cl = AR + BS:
syms P I z
A_cl = (z^2+a1*z+a0)*(z-1) + (b1*z+b0)*(P*z + I*Ts - P);
A_d = (z^2-p2*z+p1)*(z-p0);
A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
A_d_c = fliplr(coeffs(A_d, z));     % retreive coeficients

% Solve Eqn
equ1 = A_cl_c(2) == A_d_c(2); % z^2
equ2 = A_cl_c(3) == A_d_c(3); % z^1
equ3 = A_cl_c(4) == A_d_c(4); % z^0

sol = solve([equ1, equ2, equ3], [P, I]);

P = double(sol.P)
I = double(sol.I)

%% Pole Placement - PID Controller
%(PI was unsolvable so we have to add a D term w/ LP)

%--- A_m ------
w_m = 5;
zeta_m = 1;
%--- A_o ------
w_o = 5;
zeta_o = zeta_m;

p3 = 2*exp(-zeta_m*w_m*Ts)*cos(w_m*Ts*sqrt(1-(zeta_m)^2));
p2 = exp(-2*zeta_m*w_m*Ts);
p1 = 2*exp(-zeta_o*w_o*Ts)*cos(w_o*Ts*sqrt(1-(zeta_o)^2));
p0 = exp(-2*zeta_o*w_o*Ts);

fprintf('The original cont. poles are: -0.2 and -1391.1 \n');
fprintf('The chosen cont. poles are: %d and %d \n', -w_m, -w_o);

% =====Convert to associating discrete poles==========
w_m_d = exp(-zeta_m*w_m*Ts);
w_o_d = exp(-zeta_o*w_o*Ts);

chosen_disc_poles = ['The chosen disc. poles are: ', num2str(w_m_d), ' and ', num2str(w_o_d)];
disp(chosen_disc_poles)

% =========Diophantine Eqn A_cl = AR + BS==========
syms P I D N z
A_cl = (z^2+a1*z+a0)*(z-1)*(z-1+N*Ts) + (b1*z+b0)*((z-1)*(z-1+N*Ts)*P+(z-1+N*Ts)*I*Ts+(z-1)^2*D*N);
A_d = (z^2-p3*z+p2)*(z^2-p1*z+p0);
A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
A_d_c = fliplr(coeffs(A_d, z));   % retreive coeficients

% =========Solve Eqn=========
equ1 = A_cl_c(2) == A_d_c(2); % z^3
equ2 = A_cl_c(3) == A_d_c(3); % z^2
equ3 = A_cl_c(4) == A_d_c(4); % z^1
equ4 = A_cl_c(5) == A_d_c(5); % z^0

sol = solve([equ1, equ2, equ3, equ4], [P, I, D, N]);

P = double(sol.P)
I = double(sol.I)
D = double(sol.D)
N = double(sol.N)

z = tf('z', Ts);
Gc = P + I*Ts/(z-1) + D*(N)/(1+(N*Ts)/(z-1))

% =======Calculating T(z)===========
t_o = (1-p3+p2)/(b1+b0);
A_o = [1, -p1, p0];
T = t_o*A_o;
R = [1, N*Ts-2, 1-N*Ts];

fprintf('##################################### \n');
fprintf('Plant Model Gp: \n');
Gp

fprintf('##################################### \n');
fprintf('FB Controller Gc: \n');
Gc

fprintf('##################################### \n');
fprintf('FF Controller Gff: \n');
Gff = tf(T, R, Ts)

fprintf('##################################### \n');
fprintf('Full state FB Controller Gyr: \n');
Gyr = minreal(Gff*Gp/(1+Gc*Gp), 1e-3)

disp('Recall that:')
disp(chosen_disc_poles)

[num_yr, den_yr] = tfdata(Gyr);

b1_yr = num_yr{1}(1);  b0_yr = num_yr{1}(2);
a2_yr = den_yr{1}(1);  a1_yr = den_yr{1}(2);  a0_yr = den_yr{1}(3);

B_yr = [b1_yr, b0_yr];
A_yr = [a2_yr, a1_yr, a0_yr];
[z_yr, p_yr, k_yr] = tf2zp(B_yr, A_yr)

figure(3)
pzmap(Gyr)
grid on;

dc_gain_sys = dcgain(sys)
dc_gain_dominant = dcgain(sys_dominant)
dc_gain_Gyr = dcgain(Gyr)

figure(4)
step(sys/dc_gain_sys, sys_dominant/dc_gain_dominant, Gyr)
legend({'original system', 'dominant pole', 'New'})

%% Solution with balred()

b2_rd = 5.0280e-05; b1_rd = -2.0388e-04; b0_rd = 2.4500e-04;
a2_rd = 1; a1_rd = -1.4644; a0_rd = 0.4655;
B_rd = [b2_rd, b1_rd, b0_rd];
A_rd = [a2_rd, a1_rd, a0_rd];

Gp_rd = tf(B_rd, A_rd, Ts)

[z_rd, p_rd, k_rd] = tf2zp(B_rd, A_rd)

figure(5)
step(sys, Gp, Gp_rd)
legend('original system', 'dominant pole', 'model reduction');























