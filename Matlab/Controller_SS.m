T_train = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\Test_Data\Training_Data\D_48404_train.csv');
T_train.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
T_test = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\Test_Data\Testing_Data\D_48404_test.csv');
T_test.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
u_train = table2array(T_train(:, 1));   % OTSV1
y_train = table2array(T_train(:, 2));   % TV12
u_test =  table2array(T_test(:, 1));   % OTSV1
y_test =  table2array(T_test(:, 2));   % TV12
Ts = 0.4;  % second
% Training Data
data_train = iddata(y_train, u_train, Ts);
data_train.Name = 'Data_Train';
data_train.TimeUnit = 'seconds';
data_train.InputName = 'OTSV1';   data_train.InputUnit = 'Percentage';
data_train.OutputName = 'TV12';   data_train.OutputUnit = 'Celsius';
% Testing Data
data_test = iddata(y_test, u_test, Ts);
data_test.Name = 'Data_Test';
data_test.TimeUnit = 'seconds';
data_test.InputName = 'OTSV1';   data_test.InputUnit = 'Percentage';
data_test.OutputName = 'TV12';   data_test.OutputUnit = 'Celsius';
data_train = detrend(data_train);
data_test = detrend(data_test);
%
cost_func = 'NRMSE';
%%
sysSS = n4sid(data_train, 2)
Gp = tf(sysSS)
[num, den] = tfdata(Gp);
b1 = num{1}(2);
b0 = num{1}(3);
a1 = den{1}(2);
a0 = den{1}(3);
B = [b1, b0];
A = [1, a1, a0];
%% Best Avg. Model: SS
% A_ss = [[1.037, 0.4309];
%        [-0.01128, 0.8728]];
% B_ss = [-0.003958;
%         0.0011];
% C_ss = [-3.142, -0.6329];
% D_ss = 0;
% 
% [B, A] = ss2tf(A_ss, B_ss, C_ss, D_ss);
% b1 = B(2)
% b0 = B(3)
% a1 = A(2)
% a0 = A(3)
% Ts = 0.4;
% Gp = tf(B, A, Ts)
% 
% poles = pole(Gp);  % 0.9983, 0.9115
% zeros = zero(Gp);  % 0.9923

%pzmap(Gp)
%%
Con = [sysSS.B, sysSS.A*sysSS.B];
rank(Con)
%% Choose Poles (w_m, zeta_m) (w_o, zeta_o)
%--- A_m ------
w_m = 1;
zeta_m = 1;
%--- A_o ------
w_o = 1;
zeta_o = zeta_m;

p3 = 2*exp(-zeta_m*w_m*Ts)*cos(w_m*Ts*sqrt(1-(zeta_m)^2));
p2 = exp(-2*zeta_m*w_m*Ts);
p1 = 2*exp(-zeta_o*w_o*Ts)*cos(w_o*Ts*sqrt(1-(zeta_o)^2));
p0 = exp(-2*zeta_o*w_o*Ts);

% Convert to associating discrete poles
w_m_d = exp(-zeta_m*w_m*Ts);
w_o_d = exp(-zeta_o*w_o*Ts);

fprintf('The original disc. poles are: 0.9983 and 0.9115 \n');
fprintf('The chosen cont. poles are: %d and %d \n', -w_m, -w_o);
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

P = double(sol.P)   % 2.9864
I = double(sol.I)   % 0.6575
D = double(sol.D)   % 4.6858
N = double(sol.N)   % 3.2530
%% Gc
z = tf('z', Ts);
Gc = P + I*Ts/(z-1) + D*(N)/(1+(N*Ts)/(z-1));

%% Gff and Gyr
t_o = (1-p3+p2)/(b0);
A_o = [1, -p1, p0];
T = t_o*A_o;
R = [1, N*Ts-2, 1-N*Ts];
Gff = tf(T, R, Ts);

%% Entire closed loop check
Gyr = Gff*Gp/(1+Gc*Gp)
pole(Gyr);
zero(Gyr);

figure(1)
pzmap(Gyr)

figure(2)
step(Gyr)
grid on

stepinfo(Gyr)
%% The following parts are just for Simulink Usage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve b & c for Simulink
syms b c z
G_ff = (z-1)*(z-1+N*Ts)*b*P+c*N*D*(z-1)^2+I*Ts*(z-1+N*Ts);
G_ff_c = fliplr(coeffs(G_ff, z)); % retreive coeficients

equ1 = G_ff_c(1) == T(1); % z^2
equ2 = G_ff_c(2) == T(2); % z^1
equ3 = G_ff_c(3) == T(3); % z^0

sol = solve([equ1, equ3], [b, c]);

b = double(sol.b)   % 0.4666
c = double(sol.c)   % 0.1151
%% Get S for Simulink
[S, ~] = tfdata(Gc);
%% PID Tuner
P = 8.7423;
I = 0.1305;
D = 144.4523;
N = 0.44817;
b = 0.041887;
c = 0.0010044;