clc
clear all

%The ARMAX model is based on the paper "Heat Exchanger Model (Data-driven modelling)"

num = [0.0001701, 0.0001295];          % B
den = [1, -0.2492, -0.7469, 0, 0, 0];  % A
Ts = 4; %"2000s with a sampling time of 4s" .01
sys = tf(num, den, Ts)

%% Pole-Zero Cancellation
%Can use this to simplify the model if there exists pole-zero cancellation
sys = minreal(sys, 1e-3)   % in this case, it does not
figure(1)
pzmap(sys)
%% Checking  the system polse, zeros and DC gain
[ze, p, k] = tf2zp(num, den) %[zero, pole, gain]

%% Finding dominate pole - reducing order
z = tf('z', Ts);

sys_dominant = k*(z-ze)/((z-p(4))*(z-p(5)))  % remember to include the DC gain (k)

%Proving reduced-order system is still a good approximation of original
figure(2)
step(sys, sys_dominant)
legend({'original system', 'dominant pole'})

S = stepinfo(sys)
display('--------------------------')
S = stepinfo(sys_dominant)
% REDUCED IS VIABLE BECAUSE PLOT SHOWS THE SAME

%% Continuous Time Conversion
%The rest is just to double check if the converted poles are correct.
sys_cont = d2c(sys_dominant, 'tustin') %"Tustin response is faster by half a period"
% zoh was used in the past DMC courses for plant (Tustin for controller)
%NOTE: ZOH produces a order too large; must use tustin
sys_cont_exact = d2c(sys_dominant, 'zoh')

figure(3)
pzmap(sys_cont)
%zero_cont = zero(sys_cont);
%pole_cont = pole(sys_cont);

[num, den] = tfdata(sys_cont);
b1_c = (num{1}(2))/(num{1}(1));  b0_c = (num{1}(3))/(num{1}(1));
a1_c = (den{1}(2))/(den{1}(1));  a0_c = (den{1}(3))/(den{1}(1));

s = tf('s');
%sys_cont_proof = ((s-zero_cont(1))*(s-zero_cont(2)))/((s-pole_cont(1))*(s-pole_cont(2)))% Ts = 4    , P = -0.0006, -3.477  , Z = 0.5   , -3.6894 
sys_cont_proof = (s^2 + b1_c*s + b0_c)/(s^2 + a1_c*s + a0_c) %Just another way of solving

%sys_cont_proof = ((s-2)*(s+14.7574))/((s+0.0022)*(s+13.9109))   % Ts = 1    , P = -0.0022, -13.9109, Z = 2     , -14.7574
%sys_cont_proof = ((s+1475.7)*(s+200))/((s+1391.1)*(s+0.2))       % Ts = 0.01 , P = -0.2   , -1391.1 , Z = -2000 , -1475.7
%sys_cont_proof = ((s+14757)*(s+20000))/((s+13911)*(s+2))        % Ts = 0.001, P = -2     , -13911  , Z = -20000, -14757

% check dc gain
% k = DC gain of original system
k_d = dcgain(sys_cont) % DC gain for the tustin converted 
k_c = dcgain(sys_cont_proof)

figure(4)
step((sys_dominant/k_d), (sys_cont/k_d), (sys_cont_proof/k_c))
legend({'original system (discrete)', 'dominant pole (cont.)', 'convert proof (cont.)'})

%% Pole Placement - PID Controller (Cont. --> Discrete)

%--- A_m ------
w_m = 10;
zeta_m = 1;
%--- A_o ------
w_o = 5;
zeta_o = zeta_m;

s_0 = (w_m*w_o)^2/b0_c;
s_2 = 0;

%A = [1, 1;a1_c, b1_c];

%B = [2*(zeta_m*w_m + zeta_o*w_o) - a1_c - b1_c*s_2 ; 
%    (w_m^2) + (4*zeta_m*zeta_o*w_m*w_o) + (w_o^2) - s_0 - a0_c - b0_c*s_2];

%[r_0c, s_1c] =  A\B;

syms r_0c s_1c
equ1_c = r_0c + s_1c == 2*(zeta_m*w_m + zeta_o*w_o) - b1_c*s_2 - a1_c;
equ2_c = (a1_c*r_0c) + (b1_c*s_1c) == (w_m^2) + (4*zeta_m*zeta_o*w_m*w_o) + (w_o^2) - a0_c - s_0 - b0_c*s_2;
% equ3_c_1 = (a0_c*r_0) + (b0_c*s_1) == 2*w_m*w_o*(zeta_m*w_o + zeta_o*w_m) - b1_c*s_0;

sol_cont = solve([equ1_c, equ2_c], [r_0c, s_1c]);
r_0 = double(sol_cont.r_0c)
s_1 = double(sol_cont.s_1c)

%Checking to make sure solution still holds in 3rd equation:
equ3_c_1 = (a0_c*r_0) + (b0_c*s_1)
equ3_c_2 = 2*w_m*w_o*(zeta_m*w_o + zeta_o*w_m) - b1_c*s_0
Q = isAlways(equ3_c_1 == equ3_c_2)
%***NOT FEASIBLE! MUST HAVE MORE POLES THAN ZEROS***


%% Pole Placement - PI Controller (Direct - Discrete)
 
% Retreive coefficients of the plant...
Gp = sys_dominant; %discrete
[num, den] = tfdata(Gp);
b1_d = num{1}(2);  b0_d = num{1}(3);
a1_d = den{1}(2);  a0_d = den{1}(3);

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
A_cl = (z^2+a1_d*z+a0_d)*(z-1) + (b1_d*z+b0_d)*(P*z + I*Ts - P);
A_d = (z^2-p2*z+p1)*(z-p0);
A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
A_d_c = fliplr(coeffs(A_d, z));     % retreive coeficients

% Solve Eqn
equ1 = A_cl_c(2) == A_d_c(2); % z^2
equ2 = A_cl_c(3) == A_d_c(3); % z^1
equ3 = A_cl_c(4) == A_d_c(4); % z^0

sol = solve([equ1, equ2, equ3], [P, I]);

p = double(sol.P)
I = double(sol.I)


%% Pole Placement - PID Controller (Direct - Discrete)
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
A_cl = (z^2+a1_d*z+a0_d)*(z-1)*(z-1+N*Ts) + (b1_d*z+b0_d)*((z-1)*(z-1+N*Ts)*P+(z-1+N*Ts)*I*Ts+(z-1)^2*D*N);
A_d = (z^2-p3*z+p2)*(z^2-p1*z+p0);
A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
A_d_c = fliplr(coeffs(A_d, z));   % retreive coeficients

% =========Solve Eqn=========
equ1 = A_cl_c(2) == A_d_c(2); % z^3
equ2 = A_cl_c(3) == A_d_c(3); % z^2
equ3 = A_cl_c(4) == A_d_c(4); % z^1
equ4 = A_cl_c(5) == A_d_c(5); % z^0

sol = solve([equ1, equ2, equ3, equ4], [, I, D, N]);

P = double(sol.P)
I = double(sol.I)
D = double(sol.D)
N = double(sol.N)

z = tf('z', Ts);
Gc = P + I*Ts/(z-1) + D*(N)/(1+(N*Ts)/(z-1))

% =======Calculating T(z)===========
t_o = (1-p3+p2)/(b1_d+b0_d);
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

figure(5)
pzmap(Gyr)
grid on;

dc_gain_sys = dcgain(sys)
dc_gain_dominant = dcgain(sys_dominant)
dc_gain_Gyr = dcgain(Gyr)

figure(6)
step(sys/dc_gain_sys, sys_dominant/dc_gain_dominant, Gyr)
legend({'original system', 'dominant pole', 'New'})



%% Solution with balred()

b2_rd = 5.0280e-05; b1_rd = -2.0388e-04; b0_rd = 2.4500e-04;
a2_rd = 1; a1_rd = -1.4644; a0_rd = 0.4655;
B_rd = [b2_rd, b1_rd, b0_rd];
A_rd = [a2_rd, a1_rd, a0_rd];

Gp_rd = tf(B_rd, A_rd, Ts)

[z_rd, p_rd, k_rd] = tf2zp(B_rd, A_rd)

figure(7)
step(sys, Gp, Gp_rd)
legend('original system', 'dominant pole', 'model reduction');























