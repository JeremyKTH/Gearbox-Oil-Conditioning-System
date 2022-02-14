%% Original System
clc
clear all

num = [0.0001701, 0.0001295];          % B
den = [1, -0.2492, -0.7469, 0, 0, 0];  % A
Ts = 0.01;
sys = tf(num, den, Ts);
%% Reduce LTI model order using balanced truncation
System = sys; % Define System to reduce
Order = 2;
 
% Create option set for balred command
Options = balredOptions();
 
% Compute reduced order approximation
Gp = balred(System,Order,Options)

[num, den] = tfdata(Gp);
b2 = num{1}(1); b1 = num{1}(2); b0 = num{1}(3);
a2 = 1; a1 = den{1}(2); a0 = den{1}(3);

B = [b2, b1, b0];
A = [a2, a1, a0];

%% Choose Poles
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

% Convert to associating discrete poles
w_m_d = exp(-zeta_m*w_m*Ts);
w_o_d = exp(-zeta_o*w_o*Ts);

chosen_disc_poles = ['The chosen disc. poles are: ', num2str(w_m_d), ' and ', num2str(w_o_d)];
disp(chosen_disc_poles)

%% Diophantine Eqn A_cl = AR + BS
syms P I D N z
A_cl = (z^2+a1*z+a0)*(z-1)*(z-1+N*Ts) + (b1*z+b0)*((z-1)*(z-1+N*Ts)*P+(z-1+N*Ts)*I*Ts+(z-1)^2*D*N);
A_d = (z^2-p3*z+p2)*(z^2-p1*z+p0);
A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
A_d_c = fliplr(coeffs(A_d, z));   % retreive coeficients

% Solve Eqn
equ1 = A_cl_c(2) == A_d_c(2); % z^3
equ2 = A_cl_c(3) == A_d_c(3); % z^2
equ3 = A_cl_c(4) == A_d_c(4); % z^1
equ4 = A_cl_c(5) == A_d_c(5); % z^0

sol = solve([equ1, equ2, equ3, equ4], [P, I, D, N]);

P = double(sol.P)
I = double(sol.I)
D = double(sol.D)
N = double(sol.N)

%% Gc
z = tf('z', Ts);
Gc = P + I*Ts/(z-1) + D*(N)/(1+(N*Ts)/(z-1))

%%
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

%%
disp('Recall that:')
disp(chosen_disc_poles)

[num_yr, den_yr] = tfdata(Gyr);

b1_yr = num_yr{1}(1);  b0_yr = num_yr{1}(2);
a2_yr = den_yr{1}(1);  a1_yr = den_yr{1}(2);  a0_yr = den_yr{1}(3);

B_yr = [b1_yr, b0_yr];
A_yr = [a2_yr, a1_yr, a0_yr];
[z_yr, p_yr, k_yr] = tf2zp(B_yr, A_yr)

pzmap(Gyr)
grid on;