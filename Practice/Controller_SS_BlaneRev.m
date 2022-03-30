clc

%T_data = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 0.4\Training_Data\D_47414_train.csv');
%T_train.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
T_train = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 1\Training_Data\D_48404_train.csv');
T_train.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11'};

%T_data = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 0.4\Training_Data\D_47414_test.csv');
%T_test.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
T_test = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 1\Testing_Data\D_48404_test.csv');
T_test.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11'};

u_train = table2array(T_train(:, 1));   % OTSV1
y_train = table2array(T_train(:, 2));   % TV12
u_test =  table2array(T_test(:, 1));   % OTSV1
y_test =  table2array(T_test(:, 2));   % TV12

Ts = 1;  
%Ts2 = .4;

%% Training Data
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

%% Detrending Data
T_train = getTrend(data_train);
T_test = getTrend(data_test);
T_train.InputOffset = 40;
T_train.OutputOffset = 48;
T_test.InputOffset = 40;
T_test.OutputOffset = 48;

data_train = detrend(data_train, T_train);
data_test = detrend(data_test, T_test);

cost_func = 'NRMSE';


%% SS Fitting

sysSS = n4sid(data_train, 2, 'Ts', Ts);
%sysSS = ssregest(data_train, 2, 'Ts', Ts);

opt = compareOptions('InitialCondition','e');
figure(1)
compare(data_test, sysSS, opt)

Gp = tf(sysSS)

[num, den] = tfdata(Gp); %2 poles, 1 zero
b2 = num{1}(1); %z^2 (0)
b1 = num{1}(2); %zero z^1
b0 = num{1}(3); %zero z^0

a2 = den{1}(1); %z^2 (1)
a1 = den{1}(2); %pole z^1
a0 = den{1}(3); %pole 

B = [b1, b0];  % B
A = [a2, a1, a0];  % A

poles_Gp_disc = pole(Gp)
zeros_Gp_disc = zero(Gp)

figure(2)
pzmap(Gp)

Gp_cont = d2c(Gp, 'zoh') %converting to cont. 
%Gp = c2d(Gp_cont, Ts, 'zoh')

%% Best +/-40% Model: ARMAX
% b0 = 3.264e-05;
% a1 = -1.98;
% a0 = 0.9805;
% B = [b0, 0];      % B
% A = [1, a1, a0];  % A
% Ts = 0.4;
% Gp_ar = tf(B, A, Ts);
% poles = pole(Gp);  % 0.99 +/- 0.02i
% zeros = zero(Gp);  % No zeros

% pzmap(Gp_ar)
%% Choose Poles (w_m, zeta_m) (w_o, zeta_o)
%--- A_m ------
w_m = .1; % .8
zeta_m = 1;

% z^2 
p1_m = 2*exp(-zeta_m*w_m*Ts)*cos(w_m*Ts*sqrt(1-(zeta_m)^2)); %
p0_m = exp(-2*zeta_m*w_m*Ts);

%--- A_o ------
w_o = .4; %.1
zeta_o = zeta_m;

p1_o = 2*exp(-zeta_o*w_o*Ts)*cos(w_o*Ts*sqrt(1-(zeta_o)^2));
p0_o = exp(-2*zeta_o*w_o*Ts);

%NOTE: 
% -Faster observer for modelling error rejection *
% -Slower observer for filtering high frequency noise

%----------------------------------

% Convert to associating discrete poles (ONLY WORK IF ZETA = 1)
w_m_d = exp(-zeta_m*w_m*Ts);
w_o_d = exp(-zeta_o*w_o*Ts);

% fprintf('The original disc. poles are: 0.9915 +/- 0.0113i \n');
% fprintf('The chosen cont. poles are: %d and %d \n', -w_m, -w_o);
chosen_disc_poles = ['***The chosen disc. poles are: ', num2str(w_m_d), '(machine) and ', num2str(w_o_d), '(observer)***'];
disp(chosen_disc_poles)

%% Diophantine Equ - JERMY ATTEMPT (Matlab Version)

% syms P I D N z
% A_cl = (z^2+a1_ar*z+a0_ar)*(z-1)*(z-1+N*Ts) + (b0_ar*z)*((z-1)*(z-1+N*Ts)*P+(z-1+N*Ts)*I*Ts+(z-1)^2*D*N);
% A_d = (z^2-p1_m*z+p0_m)*(z^2-p1_o*z+p0_o);
% A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
% A_d_c = fliplr(coeffs(A_d, z));   % retreive coeficients

% equ1 = A_cl_c(2) == A_d_c(2); % z^3
% equ2 = A_cl_c(3) == A_d_c(3); % z^2
% equ3 = A_cl_c(4) == A_d_c(4); % z^1
% equ4 = A_cl_c(5) == A_d_c(5); % z^0
% 
% sol = solve([equ1, equ2, equ3, equ4], [P, I, D, N]);
% 
% P = double(sol.P)   % 3151.3
% I = double(sol.I)   % 696.698
% D = double(sol.D)   % 4933.2
% N = double(sol.N)   % 3.2468

%% Diophantine Equ - BLANE ATTEMPT (DMC Version)
syms S2 S1 S0 r0 z
A_cl = (z^2+a1*z+a0)*(z-1)*(z+r0) + (b1*z + b0)*(S2*z^2 + S1*z + S0);
A_d = (z^2-p1_m*z+p0_m)*(z^2-p1_o*z+p0_o);
A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
A_d_c = fliplr(coeffs(A_d, z));

equ1 = A_cl_c(2) == A_d_c(2); % z^3
equ2 = A_cl_c(3) == A_d_c(3); % z^2
equ3 = A_cl_c(4) == A_d_c(4); % z^1
equ4 = A_cl_c(5) == A_d_c(5); % z^0

sol = solve([equ1, equ2, equ3, equ4], [S2, S1, S0, r0]);

S2 = double(sol.S2);   %
S1 = double(sol.S1);   % 
S0 = double(sol.S0);   % 
r0 = double(sol.r0);   % 

%% Controller TF - Gc (S/R)
%z = tf('z', Ts);
% Gc = P + I*Ts/(z-1) + D*(N)/(1+(N*Ts)/(z-1)); %FOR JEREMY'S WAY

S = [S2, S1, S0];
R = [1, r0-1, -r0];

Gc = tf(S, R, Ts)

%% Gff (T/R) and Gyr

% %------FOR JEREMY'S WAY-------
% t_o = (1-p1_m+p0_m)/(b0_ar);
% A_o = [1, -p1_o, p0_o];
% T = t_o*A_o;
% R = [1, N*Ts-2, 1-N*Ts];
% Gff = tf(T, R, Ts);

%------FOR BLANE'S WAY-------
t_o = (1-p1_m+p0_m)/(b1 + b0);
A_o = [1, -p1_o, p0_o];
T = t_o*A_o;
R = [1, r0-1, -r0];

Gff = tf(T, R, Ts);

%% MATLAB PID TUNER SUBSTITUTION

% (EULER FORWARD)
syms P I D N 
equ1 = S2 == P + (D*N);
%equ2 = S1 == (I*Ts) + (P*N*Ts) - (2*P) - (2*D*N);
equ2 = S1 == (P*N*Ts) - (2*P) + (I*Ts) - (2*D*N);
%equ3 = S0 == (P-P*N*Ts) + (I*Ts) + (I*N*(Ts^2)) + (D*N);
equ3 = S0 == (I*N*(Ts^2)) - (I*Ts) + (D*N) - (P*N*Ts) + P;
equ4 = r0 ==(N*Ts)-1;

sol = solve([equ4, equ1, equ2, equ3], [N, P, I, D]);

N = double(sol.N)
P = double(sol.P)
I = double(sol.I)
D = double(sol.D)

% ----- Getting b & c terms for Simulink ------------
%b = "set point weight for proportional term
%c = "set point weight for deriviative term

syms b c z
% EULEER FORWARD VERSION:
G_ff_PID = (P*b*(z-1)*(z-1+N*Ts)) + (I*Ts*(z-1+N*Ts)) + (D*N*c*(z-1)^2);

% EULER BACKWARD VERSION
%G_ff_PID = (P*b*(z-1)*(z-1+N*Ts*z)) + (I*Ts*(z-1+N*Ts*z)) + (D*N*c*(z-1)^2);

G_ff_c = fliplr(coeffs(G_ff_PID, z)); % retreive coeficients

equ1 = G_ff_c(1) == T(1); % z^2
equ2 = G_ff_c(2) == T(2); % z^1
equ3 = G_ff_c(3) == T(3); % z^0

sol = solve([equ1, equ3], [b, c]);

b = double(sol.b)  
c = double(sol.c)  

%% Entire closed loop check (Gyr)
Gyr = Gff*Gp/(1+Gc*Gp)
poles_Gyr = pole(Gyr)
zeros_Gyr = zero(Gyr)

Gyr = minreal(Gyr, 1e-5)
poles_Gyr_min = pole(Gyr)
zeros_Gyr_min = zero(Gyr)

figure(3)
pzmap(Gyr)

figure(4)
step(Gyr)
grid on

stepinfo(Gyr)

%% AntiWindup Changes

T_S2 = T(1);
T_S1 = T(2);
T_S0 = T(3);

% For I terms
c0_Gc = (S0 + S1 + S2)/(1+r0);
c0_Gff = (T_S2 + T_S1 + T_S0)/(1+r0);

% For PD Terms
d1_Gc = S2;
d1_Gff = T_S2;
d0_Gc = (r0*(S1+S2) - S0)/(1+r0);
d0_Gff = (r0*(T_S1+T_S2) - T_S0)/(1+r0); 

% PD terms:
PD_SR_num = [d1_Gc, d0_Gc];
PD_TR_num = [d1_Gff, d0_Gff];
PD_den = [1, r0];

%K_ant = 2.6; %antiwindup coefficient

%% Sensitivy Analysis

S_e = 1/(1 + Gp*Gc);
T_e = 1 - S_e;

figure(5)
bode(S_e, T_e)
legend('Sensitivity Ftn', 'Comp. Sensitivity Ftn')

%% PID Tuner Parameter Method

% % SLOW (279,.72) ~ 400 sec settling:
% P_tun = 1.0665;
% I_tun = .021406;
% D_tun = 0;
% N_tun = 0;
% b_tun = 1;
% c_tun = 1;
% % k_ant = 2.8;

% MED (88.12, .6) ~ 300 sec settling:
P_tun = 5.589;
I_tun = .080445;
D_tun = 69.3573;
N_tun = .88551;
b_tun = .022431;
c_tun = 8.8241e-05;
% k_ant = 2.2;

% % FAST (23.18, .72):
% P_tun = 19.8114;
% I_tun = .1883;
% D_tun = 483.883;
% N_tun = 1.6041;
% b_tun = 1;
% c_tun = 1;
% % k_ant = .015;

% % Test (23.18, .622):
% P_tun = 18.9039;
% I_tun = .23283;
% D_tun = 352.9657;
% N_tun = 3.1466;
% b_tun = .5117;
% c_tun = .0054235;
% % k_ant = ;

% Euler Forward:
z = tf('z', Ts);
Gc_tun = P_tun + I_tun*Ts/(z-1) + D_tun*((N_tun)/(1+((N_tun*Ts)/(z-1))));
Gff_tun = b_tun*P_tun + I_tun*Ts/(z-1) + c_tun*D_tun*((N_tun)/(1+((N_tun*Ts)/(z-1))));
Gyr_tun = Gff_tun*Gp/(1+Gc_tun*Gp);
Gyr_tun = minreal(Gyr_tun, 1e-5);

pole_tun = pole(Gyr_tun)
zero_tun = zero(Gyr_tun)

figure(6)
pzmap(Gyr_tun)
figure(7)
step(Gyr_tun)



