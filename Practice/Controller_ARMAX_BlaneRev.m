clc

T_train = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\2022-Feb-24\Training_Data\D_48404_train.csv');
T_train.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
T_test = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\2022-Feb-24\Testing_Data\D_48404_test.csv');
T_test.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
u_train = table2array(T_train(:, 1));   % OTSV1
y_train = table2array(T_train(:, 2));   % TV12
u_test =  table2array(T_test(:, 1));   % OTSV1
y_test =  table2array(T_test(:, 2));   % TV12

Ts_ID = .4;
Ts = .4;  % second

%% Training Data
data_train = iddata(y_train, u_train, Ts_ID);
data_train.Name = 'Data_Train';
data_train.TimeUnit = 'seconds';
data_train.InputName = 'OTSV1';   data_train.InputUnit = 'Percentage';
data_train.OutputName = 'TV12';   data_train.OutputUnit = 'Celsius';

data_train = detrend(data_train);

%% Testing Data
data_test = iddata(y_test, u_test, Ts_ID);
data_test.Name = 'Data_Test';
data_test.TimeUnit = 'seconds';
data_test.InputName = 'OTSV1';   data_test.InputUnit = 'Percentage';
data_test.OutputName = 'TV12';   data_test.OutputUnit = 'Celsius';

data_test = detrend(data_test);

cost_func = 'NRMSE';

%% ARMAX
opt = armaxOptions;
opt.Focus = 'prediction';
opt.SearchOptions.MaxIterations = 1000;
opt.SearchOptions.Tolerance = 1e-5;

%FINAL: [2, 1, 1, 1]
na = 2;
nb = 1;
nc = 1;
nk = 1;
    
sysARMAX = armax(data_train, [na nb nc nk], opt)

opt = compareOptions('InitialCondition','e');

figure(1)
compare(data_test, sysARMAX, opt)

Gp = tf(sysARMAX)

[num, den] = tfdata(Gp); %2 poles, 1 zero
b0 = num{1}(2); %zero
a1 = den{1}(2); %pole
a0 = den{1}(3); %pole
B = [b0, 0];      % B
A = [1, a1, a0];  % A

poles_Gp_disc = pole(Gp)
zeros_Gp_disc = zero(Gp)
figure(2)
pzmap(Gp)

Gp_ar_cont = d2c(Gp, 'zoh') %converting to cont. 

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
w_m = 2;
zeta_m = 1;

p1_m = 2*exp(-zeta_m*w_m*Ts)*cos(w_m*Ts*sqrt(1-(zeta_m)^2)); %
p0_m = exp(-2*zeta_m*w_m*Ts);

%--- A_o ------
w_o = 500; %min. of 14*w_m to be stable
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
% 
% fprintf('The original disc. poles are: 0.9915 +/- 0.0113i \n');
% fprintf('The chosen cont. poles are: %d and %d \n', -w_m, -w_o);
chosen_disc_poles = ['The chosen disc. poles are: ', num2str(w_m_d), ' and ', num2str(w_o_d)];
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
A_cl = (z^2+a1*z+a0)*(z-1)*(z+r0) + (b0*z)*(S2*z^2 + S1*z + S0);
A_d = (z^2-p1_m*z+p0_m)*(z^2-p1_o*z+p0_o);
A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
A_d_c = fliplr(coeffs(A_d, z));

equ1 = A_cl_c(2) == A_d_c(2); % z^3
equ2 = A_cl_c(3) == A_d_c(3); % z^2
equ3 = A_cl_c(4) == A_d_c(4); % z^1
equ4 = A_cl_c(5) == A_d_c(5); % z^0

sol = solve([equ1, equ2, equ3, equ4], [S2, S1, S0, r0]);

S2 = double(sol.S2)   %
S1 = double(sol.S1)   % 
S0 = double(sol.S0)   % 
r0 = double(sol.r0)   % 


%% MATLAB PID TUNER SUBSTITUTION
% P = 0;
% I = 0; 
% D = 0;
% N = 0;
% 
% S2 = P + D*N;
% S1 = I*Ts + P*N*Ts-2*P-2*D*N;
% S0 = P-P*N*Ts + I*Ts + I*N*Ts^2 + D*N;
% ro = N*Ts-1;
% 
% % ----- Getting b & c terms for Simulink ------------
% %b = "set point weight for proportional term
% %c = "set point weight for deriviative term
% 
% syms b c z
% G_ff = (z-1)*(z-1+N*Ts)*b*P+c*N*D*(z-1)^2+I*Ts*(z-1+N*Ts);
% G_ff_c = fliplr(coeffs(G_ff, z)); % retreive coeficients
% 
% equ1 = G_ff_c(1) == T(1); % z^2
% equ2 = G_ff_c(2) == T(2); % z^1
% equ3 = G_ff_c(3) == T(3); % z^0
% 
% sol = solve([equ1, equ3], [b, c]);
% 
% b = double(sol.b)   % 0.4684
% c = double(sol.c)   % 0.1157
% 
% % Get S for Simulink
% [S, ~] = tfdata(Gc);

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
t_o = (1-p1_m+p0_m)/(b0);
A_o = [1, -p1_o, p0_o];
T = t_o*A_o;
R = [1, r0-1, r0];

Gff = tf(T, R, Ts);

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

K_ant = 10; %antiwindup coefficitnt

%% Sensitivy Analysis

S_e = 1/(1 + Gp*Gc);
T_e = 1 - S_e;

figure(5)
bode(S_e, T_e)
legend('Sensitivity Ftn', 'Comp. Sensitivity Ftn')




