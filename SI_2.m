%% Original System
clc

num = [0.0001701, 0.0001295];          % B
den = [1, -0.2492, -0.7469, 0, 0, 0];  % A
Ts = 0.01;
sys = tf(num, den, Ts);

% Check System Zeros, Poles, Gain
[ze, p, k] = tf2zp(num, den);

% Dominant Poles
z = tf('z', Ts);
sys_dominant = k*(z+0.7613)/((z-0.9978)*(z+0.7486));
Gp = sys_dominant;
%% Pole Placement
% 1. Retreive coefficients of the plant
[num, den] = tfdata(sys_dominant);
b1 = num{1}(2); b0 = num{1}(3);
a1 = den{1}(2); a0 = den{1}(3);

%% Choose Poles (w_m, zeta_m) (w_o, zeta_o)
%--- A_m ------
w_m = 5;  % 10 times faster, origial is at -3.4777
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
Gc = P + I*Ts/(z-1) + D*(N)/(1+(N*Ts)/(z-1));
%% Gff and Gyr
t_o = (1-p3+p2)/(b1+b0);
A_o = [1, -p1, p0];
T = t_o*A_o;
R = [1, N*Ts-2, 1-N*Ts];
Gff = tf(T, R, Ts)
%Gyr = Gff*Gp/(1+Gc*Gp)
Gyr = minreal(Gff*Gp/(1+Gc*Gp), 1e-2)
pzmap(Gyr)
grid on;
%% For Simulink solve b and c
syms b c x
Gff_L = (x-1)*(x-1+N*Ts)*b*P+c*N*D*(x-1)^2+I*Ts*(x-1+N*Ts);
Gff_R = T(1)*x^2 + T(2)*x^1 + T(3);
Gff_L_c = fliplr(coeffs(Gff_L, x)); % retreive coeficients
Gff_R_c = fliplr(coeffs(Gff_R, x)); % retreive coeficients

% Solve Eqn
equ1 = Gff_L_c(1) == Gff_R_c(1); % z^2
equ2 = Gff_L_c(2) == Gff_R_c(2); % z^1
equ3 = Gff_L_c(3) == Gff_R_c(3); % z^0

sol = solve([equ1, equ2], [b, c]); % Use any 2 of the above 3 equations reach same answers

b = double(sol.b)
c = double(sol.c)
